import os
from os import path

# For ease of debugging
import time

import numpy as np

from graphics import *
from PIL import Image as NewImage

# Declaring global variables_____________________
# Maximum value in data
max_point = 0
# Threshold value for considering a data point part of a cluster around any lattice point (rather than the space between clusters)
threshold = 0
# Window and size of window
win = 0
size_x = 512
size_y = 512

# Generate an image file for the data
def createDataImage(filename,data):

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

	# Save drawing as image file for data
	winGen.postscript(file=filename+".eps",colormode="gray")
	winGen.close()
	img = NewImage.open(filename+".eps")
	img.save(filename+".png","png")

# Show raw data as image backdrop
def setBackground(filename, data):

	# True if image does not exist
	noImage = not path.exists(filename+".png")

	if noImage:
		createDataImage(filename, data)

	# Find image file
	img = Image(Point(size_x/2,size_y/2), filename+".png")
	return img

# Identify the points corresponding to a single cluster by finding neighboring points that exceed the threshold value
def findNeighbors(data,x0,y0):

	# List of points identified in the cluster and the associated data value
	cluster = []
	# List of points to investigate, updated with neighbors from each identified point that meet the threshold
	pointList = [(x0,y0)]

	while len(pointList) > 0:
		(x,y) = pointList.pop()
		# Remove redundancies from previous updates to pointList that have already been investigated
		if data[y][x] < threshold:
			continue
		# Add point to cluster and make sure the point is not added back to the point list
		cluster.append((data[y][x],x,y))
		data[y][x] = threshold-1

		# Mark cluster
		# pt = Point(x,y)
		# pt.setOutline(color_rgb(255,0,0))
		# pt.draw(win)

		# Look through neighbors to find points to investigate
		for dy in [-1,0,1]:
			for dx in [-1,0,1]:
				xp = x+dx
				yp = y+dy
				if xp < 0 or yp < 0 or xp >= size_x or yp >= size_y:
					break
				if data[yp][xp] > threshold:
					pointList.append((xp,yp))

	return cluster,data

# Find and return all clusters of bright spots where peaks are located
def findClusters(data):
	clusters = []

	# Search through all points to start building clusters
	for y in range(size_y):
		for x in range(size_x):
			if data[y][x] > threshold:
				# Data array for the search is updated to avoid searching through points that have already been associated to clusters
				cluster,data = findNeighbors(data,x,y)
				clusters.append(cluster)

	return clusters

# Find peaks in clusters with COM calculations
def findPeaks(clusters):
	peaks = []

	for cluster in clusters:
		dataSum = 0
		xCM = 0
		yCM = 0
		for pt in cluster:
			(data,x,y) = pt
			dataSum += data
			xCM += data*x
			yCM += data*y
		xCM /= dataSum
		yCM /= dataSum
		peaks.append((xCM,yCM))

		# Mark peaks
		pt = Circle(Point(xCM,yCM),3)
		pt.setFill(color_rgb(255,0,0))
		pt.setOutline(color_rgb(255,0,0))
		pt.draw(win)

	return peaks

def main():
	filename = "CDW_GonTaS2"
	data = np.loadtxt(filename+".txt")

	# Identify width and height of data array
	global size_x, size_y
	size_x = len(data[0])
	size_y = len(data)

	# Determine data value corresponding to white
	global max_point, threshold
	max_point = max(max(line) for line in data)
	threshold = 0.2*max_point

	global win
	win = GraphWin('CDW Data', size_x, size_y)
	win.setBackground('black')

	# Show CDW images
	bg = setBackground(filename, data)
	bg.draw(win)

	# Identify all clusters around lattice points
	clusters = findClusters(data)

	# Identify the peaks within lattice point clusters
	peaks = findPeaks(clusters)

	# Close window only after mouse click in window
	win.getMouse()
	win.close()

main()