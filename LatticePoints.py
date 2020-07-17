import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import csv

from graphics import *
from PIL import Image as NewImage

# Declaring global variables_____________________
filename = "CDW_Data"
# Maximum value in data
max_point = 0
min_point = 0
# Threshold value for considering a data point part of a cluster around any lattice point (rather than the space between clusters)
threshold = 0
# Cluster radius, and maximum magnitude of all lattice vectors
averageRadius = 0
maxLatticeSpacing = 0
# Number of bins and cutoff radius for radial correlation function g(r)
bins, cutoff = 15,5
# Window and size of window
win = 0
size_x = 512
size_y = 512
scale = 1

def createDataImage(data):

	# Create window for drawing
	size_x = len(data[0])
	size_y = len(data)
	
	# Draw data
	img = NewImage.new("RGB", (size_x, size_y))
	putpixel = img.putpixel
	for y in range(size_y):
		for x in range(size_x):
			color = int((data[y][x]-min_point)/(max_point-min_point)*255)
			putpixel((x,y),(color,color,color))

	img.save(filename+".gif",'gif')

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
		# pt = Point(int(scale*x),int(scale*y))
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
		pt = Circle(Point(int(scale*xCM),int(scale*yCM)),1)
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
			# ln = Line(Point(int(scale*x0),int(scale*y0)), Point(int(scale*x),int(scale*y)))
			# ln.setOutline(color_rgb(0,255,0))
			# ln.draw(win)
		vectors += peakVectors

	# Group together similar vectors (same if data is perfect, to identify primitive vectors as the most common)
	vectors2 = vectors
	vectorClusters = []

	while len(vectors2) > 0:
		vec = vectors2.pop(0)
		if vec[1] < 0:
			(x,y) = vec
			vec = (-x,-y)
		cluster = [vec]
		vListLen = len(vectors2)

		for j in range(vListLen):
			i = vListLen-1-j
			vec2 = vectors[i]
			if vec2[1] < 0:
				(x,y) = vec2
				vec2 = (-x,-y)
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

def findRadius(cluster, peak):
	(x0,y0) = peak
	dataSum = 0
	radius = 0
	for pt in cluster:
		(data,x,y) = pt
		rSq = (x-x0)**2+(y-y0)**2
		dataSum += data
		radius += data*rSq
	radius /= dataSum
	radius = np.sqrt(radius)

	return radius,dataSum

def findAveRadius(clusters,peaks):
	dataSum = 0
	aveRadius = 0

	for i,cluster in enumerate(clusters):
		peak = peaks[i]
		(x,y) = peak
		if x+maxLatticeSpacing/2 > size_x or x-maxLatticeSpacing/2 < 0 or y+maxLatticeSpacing/2 > size_y or y-maxLatticeSpacing/2 < 0:
			continue
		radius,clusterData = findRadius(cluster, peak)
		aveRadius += clusterData*radius
		dataSum += clusterData

	aveRadius /= dataSum

	return aveRadius

def RDF(peaks,bins,cutoff):
	distances = []

	for peak in peaks:
		for peak2 in peaks:
			if peak != peak2:
				(x,y) = peak
				(x2,y2) = peak2
				distance = np.sqrt((x2-x)**2+(y2-y)**2)/maxLatticeSpacing
				if distance < cutoff:
					distances.append(distance)

	gR = [0 for i in range(bins)]
	N = len(peaks)
	binWidth = cutoff/bins
	for dist in distances:
		index = int(bins*dist/cutoff)
		a = np.floor(dist/binWidth)*binWidth
		b = a+binWidth
		weight = 1/(np.pi*(b**2-a**2))/N
		gR[index] += weight

	return gR

def storeRDF(gR):
	# Save RDF plot as an image
	plt.plot([cutoff/bins*i for i in range(bins)],gR)
	plt.title("Radial Distribution Function")
	plt.xlabel("r")
	plt.ylabel("g(r)")
	plt.savefig(filename+"_RDF.png")
	# Store RDF as a csv file
	data = [cutoff/bins*i for i in range(bins)],gR
	with open(filename+"_RDF.csv", 'w',newline='') as f:
		wr = csv.writer(f, quoting=csv.QUOTE_ALL)
		wr.writerow(zip(*data))

	img = NewImage.open(filename+"_RDF.png")
	width,height = img.size

	gRWin = GraphWin('g(r)', width, height)
	gRPlot = Image(Point(width/2,height/2), filename+"_RDF.png")
	gRPlot.draw(gRWin)
	gRWin.getMouse()
	gRWin.close()

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

	# Identify all clusters around lattice points
	clusters = findClusters(data)

	# Identify the peaks within lattice point clusters
	peaks = findPeaks(clusters)

	# Identify primitive vectors
	vectors = findLatticeVectors(peaks)
	global maxLatticeSpacing
	pLen = [np.sqrt(sum(x**2 for x in vectors[0]))]
	pLen.append(np.sqrt(sum(x**2 for x in vectors[1])))
	maxLatticeSpacing = max(pLen)
	print(maxLatticeSpacing)

	# Find approximate radius of each cluster
	global averageRadius
	averageRadius = findAveRadius(clusters, peaks)

	# Remove peaks too close to the edge
	peaks2 = []
	clusters2 = []

	for i,peak in enumerate(peaks):
		(x,y) = peak

		if x+2*averageRadius > size_x or x-2*averageRadius < 0 or y+2*averageRadius > size_y or y-2*averageRadius < 0:
			continue

		peaks2.append(peak)
		clusters2.append(clusters[i])

		# Mark peaks
		pt = Circle(Point(int(scale*x),int(scale*y)),5)
		pt.setFill(color_rgb(128,0,128))
		pt.setOutline(color_rgb(128,0,128))
		pt.draw(win)

		# Show primitive vectors off of every peak
		(dx1,dy1) = vectors[0]
		(dx2,dy2) = vectors[1]
		ln1 = Line(Point(int(scale*x),int(scale*y)),Point(int(scale*(x+dx1)),int(scale*(y+dy1))))
		ln1.setOutline(color_rgb(0,0,255))
		ln1.draw(win)
		ln2 = Line(Point(int(scale*x),int(scale*y)),Point(int(scale*(x+dx2)),int(scale*(y+dy2))))
		ln2.setOutline(color_rgb(0,0,255))
		ln2.draw(win)

	# Close window only after mouse click in window
	win.getMouse()
	win.close()

	with open(filename+"_Peaks.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for peak in peaks2:
			wr.writerow(peak)

	gR = RDF(peaks2,bins,cutoff)
	storeRDF(gR)

main()