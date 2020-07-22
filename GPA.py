import os
import sys

import numpy as np
import csv

from graphics import *
from PIL import Image as NewImage

from modules.bondFunctions import *
from modules.reciprocalLattice import *
from modules.hycht.py import *

# Declaring global variables_____________________
filename = "CDW_Data"
# Variables related to peaks
peaks = []
latticeSpacing = 0
# Window and size of window
win = 0
size_x = 512
size_y = 512
scale = 1
# Fourier transform width and height
#   FTWidth and FTHeight are used in the function retFT - more on what these are precisely in the comment for retFT
#   FTImgWidth and FTImgHeight are simply the size of the display window for the Fourier transform
FTWidth = 505
FTHeight = 505
FTImgWidth = 505
FTImgHeight = 505

braggFilteredData_return_type = "Phase"

pi = np.pi
sqrt = np.sqrt

# Generate an image file for the data
def createDataImage(data):

	# Create window for drawing
	size_x = len(data[0])
	size_y = len(data)

	# Determine data value corresponding to white in CDW image
	max_point = max(max(line) for line in data)
	min_point = min(min(line) for line in data)

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

# Stores a file for and returns the "j-matrix"
#   A j-matrix is a 2D matrix with the width and height of the CDW image
#   Each point in the j-matrix corresponds to a point in the image
#   The value of each point is the index of the closest peak as stored in filename+"_Peaks.csv"
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

# Returns the bond matrix, a 2D array that indicates how all the points are interconnected
#   The bond matrix has a width and height equal to the number of peaks
#   The indices (i,j) within the array identifies two peaks by index
#   The point (i,j) will be 1 if peaks i and j are connected and 0 if they are not
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

# Return an image that shows the Fourier transform in k-space
def FTImage(rec_data,num_peaks,return_type,extension):
	global FTImgWidth, FTImgHeight

	createRecDataImage(filename,rec_data,num_peaks,FTImgWidth,FTImgHeight,return_type,extension)

	if return_type == "Re":
		img = Image(Point(int(FTImgWidth/2.0),int(FTImgHeight/2.0)), filename+"_"+extension+"_Re.gif")
	elif return_type == "Im":
		img = Image(Point(int(FTImgWidth/2.0),int(FTImgHeight/2.0)), filename+"_"+extension+"_Im.gif")
	elif return_type == "Mag":
		img = Image(Point(int(FTImgWidth/2.0),int(FTImgHeight/2.0)), filename+"_"+extension+"_Mag.gif")

	return img

# Return an image that shows the Fourier transform in k-space
def IFTImage(braggFilteredData,return_type):
	global size_x, size_y, scale

	createReconstructedImage(filename,braggFilteredData,size_x,size_y,return_type)

	if return_type == "Re":
		img = NewImage.open(filename+"_BraggFiltered_Re.gif")
	elif return_type == "Im":
		img = NewImage.open(filename+"_BraggFiltered_Im.gif")
	elif return_type == "Phase":
		img = NewImage.open(filename+"_BraggFiltered_Phase.gif")
	else: # return_type == "Mag" by default
		img = NewImage.open(filename+"_BraggFiltered_Mag.gif")
	img = img.resize((int(size_x*scale),int(size_y*scale)),NewImage.ANTIALIAS)
	img.save(filename+"_BraggFilteredScaled.gif","gif")
	img = Image(Point(int(size_x*scale/2.0),int(size_y*scale/2.0)), filename+"_BraggFilteredScaled.gif")

	return img

def main():
	global filename, scale, braggFilteredData_return_type
	if len(sys.argv) > 1:
		filename = sys.argv[1]
	if len(sys.argv) > 2:
		scale = float(sys.argv[2])
	if len(sys.argv) > 3:
		braggFilteredData_return_type = sys.argv[3]
	if not os.path.isdir(filename):
		filename = "TestCDW_512px"
	if braggFilteredData_return_type not in ["Re","Im","Mag","Phase"]:
		braggFilteredData_return_type = "Phase"
	os.chdir(filename)
	data = np.loadtxt(filename+".txt")

	# Identify width and height of data array
	global size_x, size_y
	size_x = len(data[0])
	size_y = len(data)

	# Show CDW images
	global win
	win = GraphWin('CDW Data', int(size_x*scale), int(size_y*scale))
	win.setBackground('black')
	bg = setBackground(data)
	bg.draw(win)

	# Store discrete points for each peak in real space (if data exists)
	global peaks
	if os.path.exists(filename+"_Peaks.csv"):
		with open(filename+"_Peaks.csv",newline='') as file:
			reader = csv.reader(file,delimiter=',',quotechar='|')
			for row in reader:
				x,y = row
				x = float(x)
				y = float(y)
				peaks.append([x,y])
	for peak in peaks:
		x,y = peak
		pt = Point(scale*x,scale*y)
		mark = Circle(pt,2)
		mark.setFill(color_rgb(255,0,0))
		mark.setOutline(color_rgb(255,0,0))
		mark.draw(win)
	num_peaks = len(peaks)

	jMatrix = retJMatrix(peaks,size_x,size_y)
	bondMatrix = retTriangulation(jMatrix,num_peaks,size_x,size_y)

	# Calculate the lattice spacing (useful in determing the scale for many calculations)
	global latticeSpacing
	latticeSpacing = 0
	total_bonds = 0
	for i,peak in enumerate(peaks):
		x,y = peak
		for j,bonded in enumerate(bondMatrix[i]):
			if bonded == 1:
				x2,y2 = peaks[j]
				total_bonds += 1
				latticeSpacing += sqrt((x2-x)**2 + (y2-y)**2)
	latticeSpacing /= total_bonds

	# Generate and show Fourier transform
	global FTWidth, FTHeight
	global FTImgWidth, FTImgHeight
	rec_data = findFT(filename,peaks,size_x,size_y,FTWidth,FTHeight)
	imgFT = FTImage(rec_data,num_peaks,"Mag","FT")
	winFT = GraphWin('CDW Data - Fourier Transform', FTImgWidth, FTImgHeight)
	imgFT.draw(winFT)

	maskCenter = winFT.getMouse()
	G_x = 4*pi*(maskCenter.getX()-FTImgWidth/2.0)/FTImgWidth/sqrt(3)
	G_y = 4*pi*(maskCenter.getY()-FTImgHeight/2.0)/FTImgHeight/sqrt(3)
	print(latticeSpacing)
	print(4*pi/sqrt(3)/sqrt(G_x**2+G_y**2))
	graphCenter = Point(int(FTImgWidth/2.0),int(FTImgHeight/2.0))
	ln = Line(graphCenter,maskCenter)
	ln.setOutline(color_rgb(0,128,255))
	ln.draw(winFT)
	mark = Circle(maskCenter,3)
	mark.setFill(color_rgb(0,128,255))
	mark.setOutline(color_rgb(0,128,255))
	mark.draw(winFT)

	maskRadiusEndpoint = winFT.getMouse()
	maskRadius = sqrt((maskRadiusEndpoint.getX()-maskCenter.getX())**2+(maskRadiusEndpoint.getY()-maskCenter.getY())**2)
	radius = 4*pi*maskRadius/sqrt(FTImgWidth**2+FTImgHeight**2)/sqrt(3)
	mark = Circle(maskCenter,int(maskRadius))
	mark.setOutline(color_rgb(255,255,0))
	mark.draw(winFT)

	winFT.getMouse()
	winFT.close()

	# Apply mask on the Fourier Transform
	mask = stepFTMask
	rec_data_masked = mask(rec_data, G_x, G_y, radius)
	imgFT_masked = FTImage(rec_data_masked,num_peaks,"Mag","FT_masked")

	# Generate and show Fourier transform
	braggFilteredData = findInvFT(filename,rec_data_masked,size_x,size_y)
	imgIFT = IFTImage(braggFilteredData,braggFilteredData_return_type)
	winIFT = GraphWin('CDW Data - Inverse Fourier Transform', int(size_x*scale), int(size_y*scale))
	imgIFT.draw(winIFT)
	winIFT.getMouse()
	winIFT.close()

	win.getMouse()
	win.close()

main()