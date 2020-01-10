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

# Check if two vectors represent the same vector if the data was perfect
def vectorsSimilar(v1,v2):
	(x1,y1) = v1
	(x2,y2) = v2
	v1Len = np.sqrt(sum(i**2 for i in v1))
	v2Len = np.sqrt(sum(i**2 for i in v2))
	# False if the vector lengths are too different
	if np.abs(v1Len-v2Len)/min(v1Len,v2Len) > 0.2:
		return False
	cosThet = (x1*x2+y1*y2)/(v1Len*v2Len)
	# False if the directions of the vectors are too different
	if cosThet < 0.8:
		return False
	return True

# Find primitive lattice vectors
def findLatticeVectors(peaks):
	vectors = []

	# Identify the 6 shortest vectors connecting to neighboring peaks from each peak (stored in vectors)
	for peak in peaks:
		peakVectors = []
		for peak2 in peaks:
			if peak2 != peak:
				test = (peak2[0]-peak[0],peak2[1]-peak[1])
				if len(peakVectors) == 6:
					testLen = np.sqrt(sum(x**2 for x in test))
					index,maxLen = 0,0
					for i, vector in enumerate(peakVectors):
						vectorLen = np.sqrt(sum(x**2 for x in vector))
						if vectorLen > maxLen:
							maxLen = vectorLen
							index = i
					if maxLen > testLen:
						peakVectors[index] = test
				else:
					peakVectors.append(test)
		for vec in peakVectors:
			(x0,y0) = peak
			(dx,dy) = vec
			(x,y) = (x0+dx,y0+dy)
			# ln = Line(Point(x0,y0), Point(x,y))
			# ln.setOutline(color_rgb(0,255,0))
			# ln.draw(win)
		vectors += peakVectors

	# Group together similar vectors (same if data is perfect, to identify primitive vectors as the most common)
	vectors2 = vectors
	vectorClusters = []

	while len(vectors2) > 0:
		vec = vectors2.pop(0)
		cluster = [vec]
		vListLen = len(vectors2)

		for j in range(vListLen):
			i = vListLen-1-j
			vec2 = vectors[i]
			if vectorsSimilar(vec,vec2):
				cluster.append(vec2)
				vectors2.pop(i)

		vectorClusters.append(cluster)

	# Identify the two groups with the largest number of vectors (associated with the primitive vectors)
	maxVecClusters = []

	vectorClusterLengths = []
	for vectorCluster in vectorClusters:
		vectorClusterLengths.append(len(vectorCluster))

	i1 = vectorClusterLengths.index(max(vectorClusterLengths))
	vectorClusterLengths.pop(i1)
	maxVecClusters.append(vectorClusters.pop(i1))
	i2 = vectorClusterLengths.index(max(vectorClusterLengths))
	vectorClusterLengths.pop(i2)
	maxVecClusters.append(vectorClusters.pop(i2))

	# Take averages within vector group to settle on primitive vectors
	primitives = []
	length1 = len(maxVecClusters[0])
	length2 = len(maxVecClusters[1])
	primitives.append((sum(x for (x,y) in maxVecClusters[0])/length1,sum(y for (x,y) in maxVecClusters[0])/length1))
	primitives.append((sum(x for (x,y) in maxVecClusters[1])/length2,sum(y for (x,y) in maxVecClusters[1])/length2))

	return primitives

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

	# Identify primitive vectors
	vectors = findLatticeVectors(peaks)

	# Show primitive vectors off of every peak
	for peak in peaks:
		(x,y) = peak
		(dx1,dy1) = vectors[0]
		(dx2,dy2) = vectors[1]
		ln1 = Line(Point(x,y),Point(x+dx1,y+dy1))
		ln1.setOutline(color_rgb(0,0,255))
		ln1.draw(win)
		ln2 = Line(Point(x,y),Point(x+dx2,y+dy2))
		ln2.setOutline(color_rgb(0,0,255))
		ln2.draw(win)

	# Close window only after mouse click in window
	win.getMouse()
	win.close()

main()