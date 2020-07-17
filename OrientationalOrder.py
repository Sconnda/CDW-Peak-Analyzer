import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import math
import random
import csv

from graphics import *
from PIL import Image as NewImage

from modules.bondFunctions import findJMatrix
from modules.bondFunctions import triangulation

# Declaring global variables_____________________
filename = "CDW_Data"
# Cutoff radius for orientational order correlation function G6(r)
cutoff = 100
# Borders of image to select edge peaks
border_x = 0.1
border_y = 0.07
# Window and size of window
win = 0
size_x = 512
size_y = 512
scale = 1

sin = np.sin
cos = np.cos

def retJMatrix(peaks,size_x,size_y):
	noJMatrix = not os.path.exists(filename+"_jMatrix.csv")

	if noJMatrix:
		return findJMatrix(filename,peaks,size_x,size_y)

	jMatrix = []
	with open(filename+"_jMatrix.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			jMatrix.append(row)

	return jMatrix

def retTriangulation(jMatrix,n_peaks,size_x,size_y):
	noTriangulation = not os.path.exists(filename+"_Triangulation.csv")

	if noTriangulation:
		return triangulation(jMatrix,n_peaks,size_x,size_y)

	bondMatrix = []
	with open(filename+"_Triangulation.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			bondMatrix.append(row)

	return bondMatrix

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
	noImage = not os.path.exists(filename+".gif")

	if noImage:
		img = NewImage.new('RGB',(int(size_x*scale),int(size_y*scale)),(255,255,255))
		img.save(filename+"_Scaled.gif","gif")
		img = Image(Point(int(scale*size_x/2),int(scale*size_y/2)), filename+"_Scaled.gif")
		return img

	img = NewImage.open(filename+".gif")
	img = img.resize((int(scale*size_x),int(scale*size_y)),NewImage.ANTIALIAS)
	img.save(filename+"_Scaled.gif","gif")
	img = Image(Point(int(scale*size_x/2),int(scale*size_y/2)), filename+"_Scaled.gif")
	return img

def psi6(peak,neighbors):
	x,y = peak
	re = 0
	im = 0
	re_err = 0
	im_err = 0
	N = len(neighbors)
	for neighbor in neighbors:
		x2,y2 = neighbor
		X = x2-x
		Y = y2-y
		dX = 4
		dY = 4
		theta = np.arctan2(Y,X)
		dtheta = ((X*dY)**2+(Y*dX)**2)**0.5/(X**2+Y**2)
		re += cos(6*theta)
		im += sin(6*theta)
		re_err += (6*sin(6*theta)*dtheta)**2
		im_err += (6*cos(6*theta)*dtheta)**2
	re_err = re_err**0.5
	im_err = im_err**0.5
	re /= N
	im /= N
	re_err /= N
	im_err /= N
	return [(re,im),(re_err,im_err)]

def findNeighbors(peak_index,peaks,bondMatrix):
	neighbors = []
	for j in range(len(bondMatrix[peak_index])):
		if bondMatrix[peak_index][j] == 1:
			neighbors.append(peaks[j])
	return neighbors

def orientationalOrder(peaks,edges,bondMatrix,latticeSpacing):
	G6_r_re = [0 for i in range(cutoff)]
	G6_r_im = [0 for i in range(cutoff)]
	G6_r_re_err = [0 for i in range(cutoff)]
	G6_r_im_err = [0 for i in range(cutoff)]
	N_r = [0 for i in range(cutoff)]
	percentDone = 0
	for i,peak in enumerate(peaks):
		pDone = int(100*float(i)/len(peaks))
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
		[(re,im),(re_err,im_err)] = psi6(peak,neighbors)

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
			[(re2,im2),(re2_err,im2_err)] = psi6(peak2,neighbors2)
			g6_r = (re*re2+im*im2, re2*im-re*im2)
			g6_r_err = (((re*re2_err)**2+(re2*re_err)**2+(im*im2_err)**2+(im2*im_err)**2)**0.5, ((re2*im_err)**2+(im*re2_err)**2+(re*im2_err)**2+(im2*re_err)**2)**0.5)

			rBin = int(distance-0.5)
			G6_r_re[rBin] += g6_r[0]
			G6_r_im[rBin] += g6_r[1]
			G6_r_re_err[rBin] += g6_r_err[0]**2
			G6_r_im_err[rBin] += g6_r_err[1]**2
			N_r[rBin] += 1
	for rBin in range(cutoff):
		if N_r[rBin] > 0:
			G6_r_re_err[rBin] = G6_r_re_err[rBin]**0.5
			G6_r_im_err[rBin] = G6_r_im_err[rBin]**0.5
			G6_r_re[rBin] = G6_r_re[rBin]/N_r[rBin]
			G6_r_im[rBin] = G6_r_im[rBin]/N_r[rBin]
			G6_r_re_err[rBin] = G6_r_re_err[rBin]/N_r[rBin]
			G6_r_im_err[rBin] = G6_r_im_err[rBin]/N_r[rBin]

	return [(G6_r_re,G6_r_re_err),(G6_r_im,G6_r_im_err)]

def storeG6(G6_r,G6_r_err):
	maxLen = len(G6_r)
	G6_r0 = 2*G6_r[0]-G6_r[1]
	for i in range(maxLen):
		G6_r[i] /= G6_r0
		G6_r_err[i] /= G6_r0
	# Save G6(r) plot as an image
	plt.plot([i+1 for i in range(maxLen)],G6_r,'bo')
	plt.ylim((0, 1))
	plt.title("Orientational Order Correlation Function")
	plt.xlabel("r")
	plt.ylabel("G6(r)")
	plt.savefig(filename+"_G6.png")
	# Store G6(r) as a csv file
	with open(filename+"_G6.csv", 'w',newline='') as f:
		wr = csv.writer(f, quoting=csv.QUOTE_ALL)
		for i in range(maxLen):
			wr.writerow([i+1,G6_r[i],G6_r_err[i]])

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
	jMatrix = retJMatrix(peaks,size_x,size_y)
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
	bondMatrix = retTriangulation(jMatrix,len(peaks),size_x,size_y)
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
				total_bonds += 1
				latticeSpacing += np.sqrt((x2-x)**2 + (y2-y)**2)
		if y < border_y*size_y or y > (1-border_y)*size_y:
			edges.append(i)
		elif x < border_x*size_x or x > (1-border_x)*size_x:
			edges.append(i)

	latticeSpacing /= total_bonds
	print(latticeSpacing)

	latticeSpacing = 39.0428

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
	# G6_r_Clipped = G6_r[0][0:20]
	G6_r_Clipped = G6_r[0]
	storeG6(G6_r_Clipped[0],G6_r_Clipped[1])

main()