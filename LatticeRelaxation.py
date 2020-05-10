import os
import sys
import csv
import cv2 as cv
import numpy as np
from bisect import bisect

from modules.bondFunctions import findJMatrix
from modules.bondFunctions import triangulation

peaks = []
v_peaks = []
a_peaks = []
neighbors = [[]]
selected = -1
lattice_spacing = 15
stress_factor = 0.06

width = 512
height = 512

k = 0.001
kappa = 0.01
b = 0.01
timeScale = 1

def createPoint(x,y):
	global peaks, v_peaks, neighbors
	peaks.append((x,y))
	v_peaks.append((0,0))
	a_peaks.append((0,0))
	neighbors.append([])

def deletePoint(index):
	global peaks, v_peaks, neighbors
	del peaks[index]
	del v_peaks[index]
	del a_peaks[index]
	del neighbors[index]
	for i in range(len(neighbors)):
		for j in range(len(neighbors[i])):
			jRev = len(neighbors[i])-1-j
			if neighbors[i][jRev] == index:
				del neighbors[i][jRev]
			elif neighbors[i][jRev] > index:
				neighbors[i][jRev] -= 1

def click_event(event, x, y, flags, param):
	global peaks,selected
	if event == cv.EVENT_LBUTTONDOWN:
		if selected == -1:
			index = -1
			minDist = (width**2+height**2)**0.5
			for i,peak in enumerate(peaks):
				x2,y2 = peak
				dist = ((x-x2)**2+(y-y2)**2)**0.5
				if dist < minDist:
					minDist = dist
					index = i
			selected = index
		else:
			index = -1
			minDist = (width**2+height**2)**0.5
			for i,peak in enumerate(peaks):
				x2,y2 = peak
				dist = ((x-x2)**2+(y-y2)**2)**0.5
				if dist < minDist:
					minDist = dist
					index = i
			neighbors[index].remove(selected)
			neighbors[selected].remove(index)

			selected = -1

def genPhaseShift(stress_factor):
	print("Generating phase shifts...")
	phaseShift = np.zeros([2*width,2*height,3],float)
	ave = [0,0,0]
	randomWalk = stress_factor*lattice_spacing
	pDone = 0
	for i in range(2*width):
		if 10*int(10*(i+1)/(2*width)) > pDone:
			pDone = 10*int(10*(i+1)/(2*width))
			print(str(pDone)+"%")
		for j in range(2*height):
			for k in [0,2]:
				if i == 0 and j == 0:
					phaseShift[i][j][k] = 0.5
					ave[k] += 0.5
				elif j == 0:
					shift = phaseShift[i-1][j][k]+randomWalk*(2*np.random.sample()-1)
					phaseShift[i][j][k] = shift
					ave[k] += shift
				elif i == 0:
					shift = phaseShift[i][j-1][k]+randomWalk*(2*np.random.sample()-1)
					phaseShift[i][j][k] = shift
					ave[k] += shift
				else:
					shift1 = phaseShift[i][j-1][k]
					shift2 = phaseShift[i-1][j][k]
					weight = np.random.sample()
					shift = weight*shift1+(1-weight)*shift2+randomWalk*2**0.5*(2*np.random.sample()-1)
					phaseShift[i][j][k] = shift
					ave[k] += shift

	for k in [0,2]:
		ave[k] /= 4*width**2+4*height**2
		for i in range(2*width):
			for j in range(2*height):
				phaseShift[i][j][k] -= ave[k]

	phaseShift2 = np.zeros([width,height,3],float)
	for i in range(width):
		for j in range(height):
			for k in [0,2]:
				phaseShift2[i][j][k] = phaseShift[i+int(width/2)][j+int(height/2)][k]/lattice_spacing+0.5

	cv.imshow('phaseShift',phaseShift2)
	key = cv.waitKey(1000)
	cv.destroyAllWindows()

	return phaseShift

def initializePeaks(phaseShift):
	for i in range(int(round(5*width/lattice_spacing/3))):
		for j in range(int(round(5*height/(lattice_spacing*(3**0.5)/2)/3))):
			x0 = (i+(j%2)/2.0)*lattice_spacing-width/3
			y0 = j*lattice_spacing*3**0.5/2-height/3
			x_int = int(min(max(int(x0),-width/3),4*width/3)+width/3)
			y_int = int(min(max(int(y0),-height/3),4*height/3)+height/3)
			dx = phaseShift[x_int][y_int][0]
			dy = phaseShift[x_int][y_int][2]
			x = x0+dx
			y = y0+dy
			createPoint(x,y)

	img = np.zeros([width,height,3],np.uint8)
	for peak in peaks:
		x,y = peak
		if x < width and x >= 0 and y < height and y >= 0:
			img = cv.circle(img, (int(x),int(y)),3,(255,255,255),-1)
	cv.imshow('initialPosition',img)
	key = cv.waitKey(1000)
	cv.destroyAllWindows()

def runAnimation():
	global selected, peaks, v_peaks, a_peaks

	fourcc = cv.VideoWriter_fourcc('X','V','I','D')
	out = cv.VideoWriter("RelaxationAnimation.avi", fourcc, 20, (width,height))

	frame = 1

	while True:
		img = np.zeros([width,height,3],np.uint8)

		for i in range(len(peaks)):
			peak = peaks[i]
			x,y = peak
			vx,vy = v_peaks[i]
			ax, ay = (0,0)

			for neighbor in neighbors[i]:
				torque = 0
				x0,y0 = peaks[neighbor]
				dist = ((x-x0)**2+(y-y0)**2)**0.5

				theta1 = np.arctan2(y-y0,x-x0)
				theta_bond = 2*np.pi/len(neighbors[neighbor])

				for neighbor2 in neighbors[neighbor]:
					x2,y2 = peaks[neighbor2]
					theta2 = np.arctan2(y2-y0,x2-x0)
					theta = theta1-theta2
					theta_eq = theta_bond*round(theta/theta_bond)
					delTheta = theta-theta_eq
					torque += -kappa*delTheta

				ax += -torque/dist*np.sin(theta1)
				ay += torque/dist*np.cos(theta1)

				ax += k*(dist-lattice_spacing)*(x0-x)/dist
				ay += k*(dist-lattice_spacing)*(y0-y)/dist
				ax -= b*vx
				ay -= b*vy

			a_peaks[i] = (ax,ay)

		v_cm = [0,0]
		for i in range(len(peaks)):
			ax,ay = a_peaks[i]
			vx,vy = v_peaks[i]
			vx += timeScale*ax
			vy += timeScale*ay
			v_peaks[i] = (vx,vy)

			v_cm[0] += vx/len(peaks)
			v_cm[1] += vy/len(peaks)

		for i in range(len(v_peaks)):
			vx,vy = v_peaks[i]
			vx -= v_cm[0]
			vy -= v_cm[1]
			v_peaks[i] = (vx,vy)

		for i in range(len(peaks)):
			peak = peaks[i]
			x,y = peak
			vx,vy = v_peaks[i]
			x += timeScale*vx
			y += timeScale*vy
			peaks[i] = (x,y)

		for i,peak in enumerate(peaks):
			x,y = peak
			if i == selected:
				img = cv.circle(img, (int(x),int(y)),3,(0,0,255),-1)
			elif len(neighbors[i]) != 6:
				img = cv.circle(img, (int(x),int(y)),3,(0,255,255),-1)
			else:
				img = cv.circle(img, (int(x),int(y)),3,(255,255,255),-1)
			for neighbor in neighbors[i]:
				x2,y2 = peaks[neighbor]
				img = cv.line(img, (int(x),int(y)),(int(x2),int(y2)),(255,255,255),1)

		cv.imshow('frame',img)
		cv.setMouseCallback('frame',click_event)
		out.write(img)

		print(frame)
		frame += 1

		key = cv.waitKey(1) & 0xFF
		if key == ord('c'):
			selected = -1
		elif key == ord('q') or frame == 250:
			with open("SimulatedData_Peaks.csv", 'w',newline='') as f:
				wr = csv.writer(f)
				for peak in peaks:
					wr.writerow(peak)
			findJMatrix("SimulatedData",peaks,width,height)
			break

	cv.destroyAllWindows()

def main():
	global width,height, stress_factor
	if len(sys.argv) > 1:
		width = int(sys.argv[1])
		height = int(sys.argv[1])
	if len(sys.argv) > 2:
		stress_factor = float(sys.argv[2])
	os.chdir("SimulatedData")
	with open("SimulatedData.txt","w+") as file:
		for i in range(height):
			if i < height-1:
				str = '1.'+('\t1.' * (width-1))+'\n'
			else:
				str = '1'+('\t1' * (width-1))
			file.write(str)
	if not os.path.exists("SimulatedData_Peaks.csv"):
		phaseShift = genPhaseShift(stress_factor)
		initializePeaks(phaseShift)
		with open("SimulatedData_Peaks.csv", 'w',newline='') as f:
			wr = csv.writer(f)
			for peak in peaks:
				wr.writerow(peak)

		findJMatrix("SimulatedData",peaks,width,height)
	else:
		with open("SimulatedData_Peaks.csv",newline='') as file:
			reader = csv.reader(file,delimiter=',',quotechar='|')
			for row in reader:
				x,y = row
				x = float(x)
				y = float(y)
				createPoint(x,y)
		if not os.path.exists("SimulatedData_jMatrix.csv"):
			findJMatrix("SimulatedData",peaks,width,height)

	jMatrix = []
	with open("SimulatedData_jMatrix.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			jMatrix.append(row)

	bondMatrix = triangulation(jMatrix,len(peaks),width,height)

	print("Creating bond matrix...")
	peaksToDelete = []
	for i,row in enumerate(bondMatrix):
		for j,bond in enumerate(row):
			if bond == 1:
				neighbors[i].append(j)
		if len(neighbors[i]) == 0:
			peaksToDelete.append(i)

	print("Deleting extra points...")
	for i in range(len(peaksToDelete)):
		iRev = len(peaksToDelete)-1-i
		index = peaksToDelete[iRev]
		deletePoint(index)

	runAnimation()

main()