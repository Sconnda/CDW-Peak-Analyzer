import os
from os import path

import numpy as np

from graphics import *
from PIL import Image

def createDataImage(filename,data):

	# Determine data value corresponding to white
	max_point = max(max(line) for line in data)
	print(max_point)

	# Create window for drawing
	size_x = len(data[0])
	size_y = len(data)
	win = GraphWin('CDW Data', size_x, size_y)
	win.setBackground('black')

	# Draw data
	for y in range(size_y):
		for x in range(size_x):
			pt = Point(x,y)
			color = int(data[y][x]/max_point*255)
			color = max(color,0)
			pt.setOutline(color_rgb(color,color,color))
			pt.draw(win)

	win.postscript(file=filename+".eps",colormode="gray")
	win.close()

	img = Image.open(filename+".eps")
	img.save(filename+".png","png")

def setBackground(filename, data):
	# True if image exists
	existingImage = path.exists(filename+".png")

	if existingImage:
		size_x = len(data[0])
		size_y = len(data)

		win = GraphWin('CDW Data', size_x, size_y)
		win.setBackground('black')

		img = Image(Point(0,0), filename+".png")
		print("Imagine all the people")
	else:
		createDataImage(filename, data)

def main():
	filename = "CDW_GonTaS2"
	data = np.loadtxt(filename+".txt")

	# Show CDW images
	setBackground(filename, data)

main()