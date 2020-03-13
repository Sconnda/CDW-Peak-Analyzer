import os
import sys

import numpy as np
import math
import random
import csv

from graphics import *
from PIL import Image as NewImage

from modules.bondFunctions import findJMatrix
from modules.bondFunctions import triangulation

# Declaring global variables_____________________
filename = "CDW_Data"
# Window and size of window
win = 0
size_x = 512
size_y = 512
scale = 1
# Mode: draw Voronoi diagram or draw CDW image
drawVoronoi = True

def retJMatrix(peaks,size_x,size_y):
	noJMatrix = not os.path.exists(filename+"_jMatrix.csv")

	if noJMatrix:
		findJMatrix(filename,peaks,size_x,size_y)

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

def voronoiImage(jMatrix,num_peaks):
	noImage = not os.path.exists(filename+"_Voronoi.gif")

	if noImage:
		img = voronoi(jMatrix,num_peaks)
		img.save(filename+"_Voronoi.gif",'gif')	

# Change back to Voronoi later
	if drawVoronoi:
		img = NewImage.open(filename+"_Voronoi.gif")
		img = img.resize((int(scale*size_x),int(scale*size_y)),NewImage.ANTIALIAS)
		img.save(filename+"_Voronoi_Scaled.gif","gif")
		img = Image(Point(int(scale*size_x/2),int(scale*size_y/2)), filename+"_Voronoi_Scaled.gif")
		return img

	img = NewImage.open(filename+".gif")
	img = img.resize((int(scale*size_x),int(scale*size_y)),NewImage.ANTIALIAS)
	img.save(filename+"_Scaled.gif","gif")
	img = Image(Point(int(scale*size_x/2),int(scale*size_y/2)), filename+"_Scaled.gif")
	return img

def main():
	global filename, scale, drawVoronoi
	if len(sys.argv) > 1:
		filename = sys.argv[1]
	if len(sys.argv) > 2:
		scale = float(sys.argv[2])
	if len(sys.argv) > 3:
		if sys.argv[3] == "background":
			drawVoronoi = False
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
	img = voronoiImage(jMatrix, num_peaks)
	img.draw(win)

	for peak in peaks:
		x,y = peak
		pt = Circle(Point(scale*x,scale*y),2)
		pt.setFill(color_rgb(255,0,0))
		pt.setOutline(color_rgb(255,0,0))
		pt.draw(win)

	defects = []
	bondMatrix = triangulation(jMatrix,len(peaks),size_x,size_y)
	for j,peak in enumerate(peaks):
		x,y = peak
		num_bonds = sum(bondMatrix[j])
		if num_bonds != 6:
			defects.append(peak)
		for i,bond in enumerate(bondMatrix[j]):
			if i > j and bond == 1:
				x2,y2 = peaks[i]
				ln = Line(Point(scale*x,scale*y),Point(scale*x2,scale*y2))
				ln.setOutline(color_rgb(255,255,255))
				ln.draw(win)

	for defect in defects:
		x,y = defect
		pt = Circle(Point(scale*x,scale*y),5)
		pt.setFill(color_rgb(255,255,0))
		pt.setOutline(color_rgb(255,255,0))
		pt.draw(win)		

	win.getMouse()
	win.close()

main()