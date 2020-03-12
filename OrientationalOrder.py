import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import math
import random
import csv

from graphics import *
from PIL import Image as NewImage

# Declaring global variables_____________________
filename = "CDW_Data"
# Cutoff radius for orientational order correlation function G6(r)
cutoff = 100
# Window and size of window
win = 0
size_x = 512
size_y = 512
scale = 1

def findJMatrix(peaks):
	noJMatrix = not os.path.exists(filename+"_jMatrix.csv")

	if noJMatrix:
		num_peaks = len(peaks)
		jMatrix = [[0 for i in range(size_x)] for j in range(size_y)]
		for y in range(size_y):
			for x in range(size_x):
				dmin = math.hypot(size_x-1, size_y-1)
				j = -1
				for i in range(num_peaks):
					xPeak,yPeak = peaks[i]
					d = math.hypot(xPeak-x, yPeak-y)
					if d < dmin:
						dmin = d
						j = i
				jMatrix[y][x] = j

		with open(filename+"_jMatrix.csv", 'w',newline='') as f:
			wr = csv.writer(f)
			for row in jMatrix:
				wr.writerow(row)
	else:
		jMatrix = []
		with open(filename+"_jMatrix.csv",newline='') as file:
			reader = csv.reader(file,delimiter=',',quotechar='|')
			for row in reader:
				for i,x in enumerate(row):
					row[i] = float(x)
				jMatrix.append(row)

	return jMatrix

def voronoi(jMatrix,num_peaks):
	img = NewImage.new("RGB", (size_x, size_y))
	putpixel = img.putpixel
	colors = []
	for i in range(num_peaks):
		colors.append((0,random.randrange(255),random.randrange(255)))
	for y in range(size_y):
		for x in range(size_x):
			j = int(jMatrix[y][x])
			putpixel((x,y),colors[j])

	return img

def background():
	noImage = not os.path.exists(filename+"_Scaled.gif")

	if noImage:
		img = NewImage.new('RGB',(size_x,size_y),(255,255,255))
		img.save(filename+"_Scaled.gif","gif")
		img = Image(Point(int(scale*size_x/2),int(scale*size_y/2)), filename+"_Scaled.gif")
		return img

	img = NewImage.open(filename+".gif")
	img = img.resize((int(scale*size_x),int(scale*size_y)),NewImage.ANTIALIAS)
	img.save(filename+"_Scaled.gif","gif")
	img = Image(Point(int(scale*size_x/2),int(scale*size_y/2)), filename+"_Scaled.gif")
	return img

def triangulation(jMatrix,num_peaks):
	bondMatrix = [[0 for i in range(num_peaks)] for j in range(num_peaks)]
	for y in range(size_y):
		for x in range(size_x):
			j = int(jMatrix[y][x])
			for dy in [-1,0,1]:
				for dx in [-1,0,1]:
					if x+dx < 0 or x+dx >= size_x or y+dy < 0 or y+dy >= size_y:
						continue
					i = int(jMatrix[y+dy][x+dx])
					if i != j:
						bondMatrix[j][i] = 1
						bondMatrix[i][j] = 1
	return bondMatrix

def psi6(peak,neighbors):
	x,y = peak
	re = 0
	im = 0
	N = len(neighbors)
	for neighbor in neighbors:
		x2,y2 = neighbor
		angle = np.arctan2(y2-y,x2-x)
		re += np.cos(6*angle)
		im += np.sin(6*angle)
	re /= N
	im /= N
	return (re,im)

def findNeighbors(peak_index,peaks,bondMatrix):
	neighbors = []
	for j in range(len(bondMatrix[peak_index])):
		if bondMatrix[peak_index][j] == 1:
			neighbors.append(peaks[j])
	return neighbors

def orientationalOrder(peaks,edges,bondMatrix,latticeSpacing):
	G6_r_re = [0 for i in range(cutoff)]
	G6_r_im = [0 for i in range(cutoff)]
	N_r = [0 for i in range(cutoff)]
	percentDone = 0
	for i,peak in enumerate(peaks):
		pDone = 10*int(10*float(i)/len(peaks))
		if pDone > percentDone:
			percentDone = pDone
			print(str(pDone)+"% Done")

		isEdge = False
		for j in edges:
			if i == j:
				isEdge = True
		if isEdge:
			continue

		neighbors = findNeighbors(i,peaks,bondMatrix)
		re,im = psi6(peak,neighbors)

		(x,y) = peak
		for j,peak2 in enumerate(peaks):
			if peak == peak2:
				continue

			isEdge = False
			for k in edges:
				if j == k:
					isEdge = True
			if isEdge:
				continue

			(x2,y2) = peak2
			distance = np.sqrt((x2-x)**2+(y2-y)**2)/latticeSpacing
			if distance >= cutoff or distance < 0.5:
				continue
			neighbors2 = findNeighbors(j,peaks,bondMatrix)
			re2,im2 = psi6(peak2,neighbors2)
			g6_r = (re*re2+im*im2, re2*im-re*im2)

			rBin = int(distance-0.5)
			G6_r_re[rBin] += g6_r[0]
			G6_r_im[rBin] += g6_r[1]
			N_r[rBin] += 1

	for rBin in range(cutoff):
		if N_r[rBin] > 0:
			G6_r_re[rBin] = G6_r_re[rBin]/N_r[rBin]
			G6_r_im[rBin] = G6_r_im[rBin]/N_r[rBin]

	return (G6_r_re,G6_r_im)

def storeG6(G6_r):
	# Save G6(r) plot as an image
	maxLen = len(G6_r)
	plt.plot([i+1 for i in range(maxLen)],G6_r,'bo')
	plt.ylim((0, 1))
	plt.title("Orientational Order Correlation Function")
	plt.xlabel("r")
	plt.ylabel("G6(r)")
	plt.savefig(filename+"_G6.png")
	# Store G6(r) as a csv file
	data = [i+1 for i in range(maxLen)],G6_r
	with open(filename+"_G6.csv", 'w',newline='') as f:
		wr = csv.writer(f, quoting=csv.QUOTE_ALL)
		wr.writerow(G6_r)

	img = NewImage.open(filename+"_G6.png")
	width,height = img.size

	G6Win = GraphWin('G6(r)', width, height)
	G6Plot = Image(Point(width/2,height/2), filename+"_G6.png")
	G6Plot.draw(G6Win)
	G6Win.getMouse()
	G6Win.close()

def main():
	global filename, scale, drawVoronoi
	if len(sys.argv) > 1:
		filename = sys.argv[1]
	if len(sys.argv) > 2:
		scale = float(sys.argv[2])
	if not os.path.isdir(filename):
		filename = "TestCDW_512px"
	os.chdir(filename)

	data = np.loadtxt(filename+".txt")
	global size_x, size_y
	size_x = len(data[0])
	size_y = len(data)

	peaks = []
	with open(filename+"_Peaks.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			x,y = row
			x = float(x)
			y = float(y)
			peaks.append([x,y])
	num_peaks = len(peaks)

	global win
	jMatrix = findJMatrix(peaks)
	win = GraphWin('CDW Data', int(size_x*scale), int(size_y*scale))
	img = background()
	img.draw(win)

	for peak in peaks:
		x,y = peak
		pt = Circle(Point(scale*x,scale*y),2)
		pt.setFill(color_rgb(255,0,0))
		pt.setOutline(color_rgb(255,0,0))
		pt.draw(win)

	defects = []
	bondMatrix = triangulation(jMatrix,len(peaks))
	for j,peak in enumerate(peaks):
		x,y = peak
		num_bonds = sum(bondMatrix[j])
		if num_bonds != 6:
			defects.append(peak)
		# for i,bond in enumerate(bondMatrix[j]):
		# 	if i > j and bond == 1:
		# 		x2,y2 = peaks[i]
		# 		ln = Line(Point(scale*x,scale*y),Point(scale*x2,scale*y2))
		# 		ln.setOutline(color_rgb(255,255,255))
		# 		ln.draw(win)

	latticeSpacing = 0
	total_bonds = 0
	edges = []
	for i,peak in enumerate(peaks):
		# bondAngles = []
		x,y = peak
		for j,bonded in enumerate(bondMatrix[i]):
			if bonded == 1:
				x2,y2 = peaks[j]
		# 		angle = np.arctan2(y2-y,x2-x)
		# 		bondAngles.append(angle)
				total_bonds += 1
				latticeSpacing += np.sqrt((x2-x)**2 + (y2-y)**2)
		# bondAngles = sorted(bondAngles)
		# maxAdjAngle = 0
		# for j in range(len(bondAngles)):
		# 	angle1 = bondAngles[j]
		# 	if j < len(bondAngles)-1:
		# 		angle2 = bondAngles[j+1]
		# 	else:
		# 		angle2 = bondAngles[0]+2*np.pi
		# 	adjAngle = angle2-angle1
		# 	if adjAngle > maxAdjAngle:
		# 		maxAdjAngle = adjAngle
		# if maxAdjAngle > 0.6*np.pi:
		# 	edges.append(i)
		if y < 0.05*size_y or y > 0.95*size_y:
			edges.append(i)
		elif x < 0.05*size_x or x > 0.95*size_x:
			edges.append(i)

	latticeSpacing /= total_bonds

	for defect in defects:
		x,y = defect
		pt = Circle(Point(scale*x,scale*y),5)
		pt.setFill(color_rgb(255,255,0))
		pt.setOutline(color_rgb(255,255,0))
		pt.draw(win)		

	for j in range(len(edges)):
		i = edges[j]
		x,y = peaks[i]
		pt = Circle(Point(scale*x,scale*y),5)
		pt.setFill(color_rgb(0,255,0))
		pt.setOutline(color_rgb(0,255,0))
		pt.draw(win)


	G6_r = orientationalOrder(peaks,edges,bondMatrix,latticeSpacing)
	win.getMouse()
	win.close()

	# G6_r_Clipped = []
	# G6_0 = 0
	# for i in range(20):
	# 	G6_i = G6_r[0][i]
	# 	if i == 0:
	# 		G6_2 = G6_r[0][1]
	# 		G6_0 = 2*G6_i-G6_2
	# 	if G6_i != 0:
	# 		G6_r_Clipped.append(G6_i/G6_0)
	G6_r_Clipped = G6_r[0][0:20]
	storeG6(G6_r_Clipped)

main()