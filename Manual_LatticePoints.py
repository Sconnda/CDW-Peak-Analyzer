import os
import sys
import keyboard

import numpy as np
import csv

from graphics import *
from PIL import Image as NewImage

# Declaring global variables_____________________
filename = "CDW_Data"
# Peak searching variables
searching = True
peaks = []
# Window and size of window
win = 0
size_x = 512
size_y = 512
scale = 1

# Generate an image file for the data
def createDataImage(data):

	# Create window for drawing
	size_x = len(data[0])
	size_y = len(data)
	winGen = GraphWin('Generating CDW Image...', size_x, size_y)
	winGen.setBackground('black')

	# Draw data
	for y in range(size_y):
		for x in range(size_x):
			pt = Point(x,y)
			color = int((data[y][x]-min_point)/(max_point-min_point)*255)
			pt.setOutline(color_rgb(color,color,color))
			pt.draw(winGen)

	winGen.getMouse()

	# Save drawing as image file for data
	winGen.postscript(file=filename+".eps",colormode="gray")
	winGen.close()
	img = NewImage.open(filename+".eps")
	img = img.resize((size_x,size_y),NewImage.ANTIALIAS)
	img.save(filename+".gif","gif")

# Show raw data as image backdrop
def setBackground(data):

	# True if image does not exist
	noImage = not os.path.exists(filename+".gif")

	if noImage:
		createDataImage(data)

	# Find image file
	img = NewImage.open(filename+".gif")
	img = img.resize((int(scale*size_x),int(scale*size_y)),NewImage.ANTIALIAS)
	img.save(filename+"_Scaled.gif","gif")
	img = Image(Point(int(scale*size_x/2),int(scale*size_y/2)), filename+"_Scaled.gif")
	return img

def undo():
	if len(peaks) > 0:
		x,y = peaks.pop()
	print(peaks)

def save():
	global searching
	searching = False

def main():
	global filename, scale
	if len(sys.argv) > 1:
		filename = sys.argv[1]
	if len(sys.argv) > 2:
		scale = float(sys.argv[2])
	if not os.path.isdir(filename):
		filename = "TestCDW_512px"
	os.chdir(filename)
	data = np.loadtxt(filename+".txt")

	# Identify width and height of data array
	global size_x, size_y
	size_x = len(data[0])
	size_y = len(data)

	# Determine data value corresponding to white
	global max_point, min_point, threshold
	max_point = max(max(line) for line in data)
	min_point = min(min(line) for line in data)
	threshold = 0.6*(max_point-min_point)+min_point

	global win
	win = GraphWin('CDW Data', int(size_x*scale), int(size_y*scale))
	win.setBackground('black')

	# Show CDW images
	bg = setBackground(data)
	bg.draw(win)

	keyboard.add_hotkey('s',lambda: save())
	keyboard.add_hotkey('z',lambda: undo())
	while searching:
		pt = win.getMouse()

		mark = Circle(pt,1)
		mark.setFill(color_rgb(255,0,0))
		mark.setOutline(color_rgb(255,0,0))
		mark.draw(win)

		if searching == True:
			peaks.append((pt.getX()/scale,pt.getY()/scale))
			print(peaks)

	with open(filename+"_Peaks.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for peak in peaks:
			wr.writerow(peak)

	win.close()

main()