import os
import csv
import cv2 as cv
import numpy as np
from bisect import bisect

peaks = []
v_peaks = []
a_peaks = []
maxBonds = []
neighbors = [[]]
peakClosed = []
selected = -1
firstUnbonded = 0
lattice_spacing = 15
defectDensity = 0.1
period = 100

width = 512
height = 512

k1 = 0.001
k2 = 0.0003
k3 = 0.001
kappa = 0.01
b = 0.01
pulse = 0.1
pulse2 = 0.7
timeScale = 5
showImage = False

os.chdir("SimulatedData")

def createPoint(x,y,n):
	global peaks, v_peaks, neighbors
	peaks.append((x,y))
	v_peaks.append((0,0))
	a_peaks.append((0,0))
	maxBonds.append(n)
	neighbors.append([])
	peakClosed.append(False)

def click_event(event, x, y, flags, param):
	global peaks,selectedx
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
			showImage = True
		else:
			index = -1
			minDist = (width**2+height**2)**0.5
			for i,peak in enumerate(peaks):
				x2,y2 = peak
				dist = ((x-x2)**2+(y-y2)**2)**0.5
				if dist < minDist:
					minDist = dist
					index = i
			if not (selected in neighbors[index]):
				neighbors[index].append(selected)
				if len(neighbors[index]) > maxBonds[index]:
					maxBonds[index] = len(neighbors[index])
			if not (index in neighbors[selected]):
				neighbors[selected].append(index)
				if len(neighbors[selected]) > maxBonds[selected]:
					maxBonds[selected] = len(neighbors[selected])
			selected = -1
			showImage = True
	elif event == cv.EVENT_RBUTTONDOWN:
		key = cv.waitKey(0) & 0xFF
		if key == ord('5'):
			createPoint(x,y,5)
		if key == ord('6'):
			createPoint(x,y,6)
		elif key == ord('7'):
			createPoint(x,y,7)

# createPoint(256,256,6)

phaseShift = np.zeros([width,height,3],float)
ave = [0,0,0]
randomWalk = 0.04*lattice_spacing
for j in range(height):
	for i in range(width):
		iRev = width-1-i
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
	ave[k] /= 512**2
	for i in range(width):
		for j in range(height):
			phaseShift[i][j][k] -= ave[k]

phaseShift2 = np.zeros([width,height,3],float)
for i in range(width):
	for j in range(height):
		for k in [0,2]:
			phaseShift2[i][j][k] = phaseShift[i][j][k]/lattice_spacing+0.5

cv.imshow('phaseShift',phaseShift2)
key = cv.waitKey(0)
cv.destroyAllWindows()

for i in range(int(round(width/lattice_spacing))):
	for j in range(int(round(height/(lattice_spacing*(3**0.5)/2)))):
		x0 = (i+(j%2)/2.0)*lattice_spacing
		y0 = j*lattice_spacing*3**0.5/2
		x_int = min(int(x0),width)
		y_int = min(int(y0),height)
		dx = phaseShift[x_int][y_int][0]
		dy = phaseShift[x_int][y_int][2]
		x = x0+dx
		y = y0+dy
		createPoint(x,y,6)

img = np.zeros([width,height,3],np.uint8)
for peak in peaks:
	x,y = peak
	img = cv.circle(img, (int(x),int(y)),3,(255,255,255),-1)
cv.imshow('initialPosition',img)
key = cv.waitKey(0)
cv.destroyAllWindows()

with open("SimulatedData_Peaks.csv", 'w',newline='') as f:
	wr = csv.writer(f)
	for peak in peaks:
		wr.writerow(peak)

t = 0

img = np.zeros([width,height,3],np.uint8)

while True:
	showImage = True
	if (t % period == 0):
		img = np.zeros([width,height,3],np.uint8)
		# for i in range(max(0,len(peaks)-20),len(peaks)):
		# 	peak = peaks[i]
		# 	x,y = peak
		# 	if len(neighbors[i]) == 0:
		# 		theta_exp = 2*np.pi/maxBonds[i]
		# 		for j in range(maxBonds[i]):
		# 			theta = j*theta_exp
		# 			n = 6
		# 			defect = np.random.sample()
		# 			if (defect < defectDensity):
		# 				n += int(2*round(np.random.sample())-1)
		# 			createPoint(x+lattice_spacing*np.cos(theta),y+lattice_spacing*np.sin(theta),n)
		# 		firstUnbonded = i
		# 		break
		# 	elif not peakClosed[i]:
		# 		angles = []
		# 		for neighbor in neighbors[i]:
		# 			x2,y2 = peaks[neighbor]
		# 			theta = np.arctan2(y2-y,x2-x)
		# 			angles.append(theta)
		# 		angles.sort()
		# 		dAngles = []
		# 		for j in range(len(angles)):
		# 			if j == len(angles)-1:
		# 				dAngles.append(angles[0]+2*np.pi-angles[j])
		# 			else:
		# 				dAngles.append(angles[j+1]-angles[j])
		# 		theta_exp = 2*np.pi/maxBonds[i]
		# 		for j,dAngle in enumerate(dAngles):
		# 			nBonds = int(round(dAngle/theta_exp))-1
		# 			midPoint = angles[j]+dAngle/2
		# 			bondsCreated = 0
		# 			for k in range(nBonds):
		# 				theta = midPoint+(k-(nBonds-1)/2)*theta_exp
		# 				n = 6
		# 				defect = np.random.sample()
		# 				if (defect < defectDensity):
		# 					n += int(2*round(np.random.sample())-1)
		# 				xp = x+lattice_spacing*np.cos(theta)
		# 				yp = y+lattice_spacing*np.sin(theta)
		# 				overlay = -1
		# 				for l,peak2 in enumerate(peaks):
		# 					if i == l:
		# 						continue
		# 					x3,y3 = peak2
		# 					dist = ((xp-x3)**2+(yp-y3)**2)**0.5
		# 					if dist < 0.5*lattice_spacing:
		# 						overlay = l
		# 				if overlay == -1:
		# 					createPoint(xp,yp,n)
		# 					bondsCreated += 1
		# 				else:
		# 					if not (overlay in neighbors[i]):
		# 						neighbors[i].append(overlay)
		# 						if len(neighbors[i]) > maxBonds[i]:
		# 							maxBonds[i] = len(neighbors[i])
		# 					if not (i in neighbors[overlay]):
		# 						neighbors[overlay].append(i)
		# 						if len(neighbors[overlay]) > maxBonds[overlay]:
		# 							maxBonds[overlay] = len(neighbors[overlay])

		# 			if bondsCreated == 0:
		# 				peakClosed[i] = True
		# 				# maxBonds[i] = len(neighbors[i])
		# 		firstUnbonded = i
		# 		break
		# # theta = 2*np.pi*np.random.sample()
		# # x = 256*(1+np.cos(theta))
		# # y = 256*(1+np.sin(theta))
		# # n = 6
		# # if (t % (period*defect_period) == 0):
		# # 	n += int(2*round(np.random.sample()-1))
		# # peak = createPoint(x,y,n)
		# # dist = ((256-x)**2+(256-y)**2)**0.5
		# # v_peaks[-1] = (pulse*(256-x)/dist,pulse*(256-y)/dist)
	elif showImage:
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
			theta_bond = 2*np.pi/maxBonds[neighbor]

			for neighbor2 in neighbors[neighbor]:
				x2,y2 = peaks[neighbor2]
				theta2 = np.arctan2(y2-y0,x2-x0)
				theta = theta1-theta2
				theta_eq = theta_bond*round(theta/theta_bond)
				delTheta = theta-theta_eq
				torque += -kappa*delTheta

			ax += -torque/dist*np.sin(theta1)
			ay += torque/dist*np.cos(theta1)

		for j,peak2 in enumerate(peaks):
			if j == i:
				continue
			x2,y2 = peak2
			dist = ((x-x2)**2+(y-y2)**2)**0.5
			if dist > 4*lattice_spacing:
				continue
			if (j in neighbors[i]):

				ax += k1*(dist-lattice_spacing)*(x2-x)/dist
				ay += k1*(dist-lattice_spacing)*(y2-y)/dist
				ax -= b*vx
				ay -= b*vy
			# elif len(neighbors[j]) < maxBonds[j] and len(neighbors[i]) < maxBonds[i]:
			# 	dist = ((x-x2)**2+(y-y2)**2)**0.5
			# 	ax += -k2*(x2-x)/(dist**2)
			# 	ay += -k2*(y2-y)/(dist**2)
		a_peaks[i] = (ax,ay)

	for i in range(len(peaks)):
		ax,ay = a_peaks[i]
		vx,vy = v_peaks[i]
		vx += timeScale*ax
		vy += timeScale*ay
		v_peaks[i] = (vx,vy)

	for i in range(len(v_peaks)):
		vx,vy = v_peaks[i]
		vx0,vy0 = v_peaks[0]
		vx -= vx0
		vy -= vy0
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
		if (t % period == 0) or showImage:
			if i == selected:
				img = cv.circle(img, (int(x),int(y)),3,(0,0,255),-1)
			elif i == firstUnbonded:
				img = cv.circle(img, (int(x),int(y)),3,(0,255,0),-1)
			elif maxBonds[i] != 6:
				img = cv.circle(img, (int(x),int(y)),3,(0,255,255),-1)
			else:
				img = cv.circle(img, (int(x),int(y)),3,(255,255,255),-1)

		for j, peak2 in enumerate(peaks):
			if j == i:
				continue
			x2,y2 = peak2
			if (j in neighbors[i]):
				if (t % period == 0) or showImage:
					img = cv.line(img, (int(x),int(y)),(int(x2),int(y2)),(255,255,255),1)				
				continue
			dist = ((x-x2)**2+(y-y2)**2)**0.5
			commonNeighbors = False
			equilibriumLatticeSpacings = []
			for neighbor in neighbors[i]:
				if j in neighbors[neighbor]:
					theta = 2*np.pi/maxBonds[neighbor]
					commonNeighbors = True
					equilibrium = lattice_spacing*(2-2*np.cos(theta))**0.5
					equilibriumLatticeSpacings.append(equilibrium)
			adjusted_lattice_spacing = lattice_spacing
			if commonNeighbors:
				adjusted_lattice_spacing = 1
				for equilibrium in equilibriumLatticeSpacings:
					adjusted_lattice_spacing *= equilibrium
				adjusted_lattice_spacing = adjusted_lattice_spacing**(1.0/len(equilibriumLatticeSpacings))
			if dist < adjusted_lattice_spacing*1.3 and len(neighbors[i]) < maxBonds[i] and len(neighbors[j]) < maxBonds[j]:
				possibleBond = True
				theta0 = np.arctan2(y2-y,x2-x)
				angles = []
				for neighbor in neighbors[i]:
					x3,y3 = peaks[neighbor]
					theta = np.arctan2(y3-y,x3-x)
					angles.append((theta,neighbor))
				if len(angles) > 0:
					angles.sort()
					index = bisect(angles,(theta,j))
					if index == 0 or index == len(angles):
						if angles[len(angles)-1][1] in neighbors[angles[0][1]]:
							possibleBond = False
					else:
						if angles[index-1][1] in neighbors[angles[0][1]]:
							possibleBond = False
					if not possibleBond:
						continue
				neighbors[i].append(j)
				neighbors[j].append(i)
				if len(neighbors[i]) == maxBonds[i]:
					peakClosed[i] = True
				if len(neighbors[j]) == maxBonds[j]:
					peakClosed[j] = True

	cv.imshow('frame',img)
	cv.setMouseCallback('frame',click_event)

	showImage = False

	key = cv.waitKey(1) & 0xFF
	if key == ord('p'):
		for i,peak in enumerate(peaks):
			if i == 0:
				continue
			x,y = peak
			dist = ((width/2-x)**2+(height/2-y)**2)**0.5
			vx,vy = v_peaks[i]
			vx += pulse*(width/2-x)/dist
			vy += pulse*(height/2-y)/dist
			v_peaks[i] = vx,vy
	elif key == ord('c'):
		selected = -1
	elif key == ord('i'):
		showImage = True
	elif key == ord(' '):
		if selected >= 0:
			for neighbor in neighbors[selected]:
				neighbors[neighbor].remove(selected)
			neighbors[selected] = []
			x,y = peaks[selected]
			for i,peak in enumerate(peaks):
				if i == selected:
					continue
				x2,y2 = peak
				dist = ((x2-x)**2+(y2-y)**2)**0.5
				vx,vy = v_peaks[i]
				vx += pulse2*(x2-x)/dist
				vy += pulse2*(y2-y)/dist
				v_peaks[i] = vx,vy
			selected = -1
	elif key == ord('q'):
		with open("SimulatedData_Peaks.csv", 'w',newline='') as f:
			wr = csv.writer(f)
			for peak in peaks:
				wr.writerow(peak)
		jMatrix = np.zeros([len(peaks),len(peaks)],int)
		for i in range(len(peaks)):
			for j in neighbors[i]:
				jMatrix[i][j] = 1
				jMatrix[j][i] = 1
		with open("SimulatedData_jMatrix.csv", 'w',newline='') as f:
			wr = csv.writer(f)
			for row in jMatrix:
				wr.writerow(row)
		break

	t += 1

cv.destroyAllWindows()