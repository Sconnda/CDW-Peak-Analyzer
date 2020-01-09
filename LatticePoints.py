import os
from os import path

import numpy as np

from graphics import *
from PIL import Image as NewImage

threshold = 0
size_x = 512
size_y = 512

def createDataImage(filename,data):

	# Determine data value corresponding to white
	max_point = max(max(line) for line in data)
	print(max_point)

	# Create window for drawing
	size_x = len(data[0])
	size_y = len(data)
	winGen = GraphWin('Generating CDW Image...', size_x, size_y)
	winGen.setBackground('black')

	# Draw data
	for y in range(size_y):
		for x in range(size_x):
			pt = Point(x,y)
			color = int(data[y][x]/max_point*255)
			color = max(color,0)
			pt.setOutline(color_rgb(color,color,color))
			pt.draw(winGen)

	winGen.postscript(file=filename+".eps",colormode="gray")
	winGen.close()

	img = NewImage.open(filename+".eps")
	img.save(filename+".png","png")

def setBackground(filename, data):
	# True if image does not exist
	noImage = not path.exists(filename+".png")

	if noImage:
		createDataImage(filename, data)

	img = Image(Point(size_x/2,size_y/2), filename+".png")
	return img

def findNeighbors(data,win,x,y):
	peak = [[data[y][x],x,y]]
	data[y][x] = 0
	pt = Point(x,y)
	pt.setOutline(color_rgb(255,0,0))
	pt.draw(win)

	for dy in [-1,0,1]:
		for dx in [-1,0,1]:
			if data[y+dy][x+dx] > threshold:
				peak2,data = findNeighbors(data,win,x+dx,y+dy)
				peak += peak2
				print(data[y][x])
				time.sleep(0.1)

	return peak,data

def findPeaks(data, win):
	data2 = data
	peaks = []

	# !!Change back to size_x, size_y
	for y in range(50):
		for x in range(50):
			if data2[y][x] > threshold:
				peak,data2 = findNeighbors(data2, win, x,y)
				peaks.append(peak)

	print(peaks)
	return 0

def main():
	filename = "CDW_GonTaS2"
	data = np.loadtxt(filename+".txt")

	size_x = len(data[0])
	size_y = len(data)

	win = GraphWin('CDW Data', size_x, size_y)
	win.setBackground('black')

	# Show CDW images
	bg = setBackground(filename, data)
	bg.draw(win)

	# findPeaks(data, win)

	win.getMouse()
	win.close()

main()