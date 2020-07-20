import os
import sys

import numpy as np
import csv

from graphics import *
from PIL import Image as NewImage

from modules.bondFunctions import *
from modules.reciprocalLattice import *

# Declaring global variables_____________________
filename = "CDW_Data"
min_point = 0
max_point = 0
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

pi = np.pi
sqrt = np.sqrt

# Generate an image file for the data
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

# Return a 2D array of the Fourier transformed image in k-space
#   The limits of the 2D array are always [-4pi/a, 4pi/a] in both the k_x and k_y directions
#   FTWidth and FTHeight correspond to the resolution you want in the Fourier transform
#   To confirm, this means the center of the image corresponds to k = (0,0)
def retFT(peaks,latticeSpacing):
	global FTImgWidth, FTImgHeight
	noFT = not os.path.exists(filename+"_ReciprocalLattice.csv") or not os.path.exists(filename+"_ReciprocalLattice_Im.csv")

	if noFT:
		return findFT(filename,peaks,size_x,size_y,latticeSpacing,FTWidth,FTHeight)

	rec_data_re = []
	rec_data_im = []
	with open(filename+"_ReciprocalLattice.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			rec_data_re.append(row)

	with open(filename+"_ReciprocalLattice_Im.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			rec_data_im.append(row)

	return rec_data_re, rec_data_im

# Return an image that shows the Fourier transform in k-space
def FTImage(rec_data,num_peaks,return_type):
	global FTImgWidth, FTImgHeight
	if return_type == "Re":
		noFTImage = not os.path.exists(filename+"_FT.gif")
	elif return_type == "Im":
		noFTImage = not os.path.exists(filename+"_FT_Im.gif")
	elif return_type == "Mag":
		noFTImage = not os.path.exists(filename+"_FT_Mag.gif")

	if noFTImage:
		createRecDataImage(filename,rec_data,num_peaks,FTImgWidth,FTImgHeight,return_type)

	if return_type == "Re":
		img = NewImage.open(filename+"_FT.gif")
	elif return_type == "Im":
		img = NewImage.open(filename+"_FT_Im.gif")
	elif return_type == "Mag":
		img = NewImage.open(filename+"_FT_Mag.gif")
	img = img.resize((FTImgWidth,FTImgHeight),NewImage.ANTIALIAS)
	img.save(filename+"_FTScaled.gif","gif")
	img = Image(Point(int(FTImgWidth/2.0),int(FTImgHeight/2.0)), filename+"_FTScaled.gif")

	return img

# Return a 2D array of the inverse Fourier transformed image in real space
#   The size of the array is equal to the size of the orignal image
def retInvFT(rec_data,latticeSpacing):
	global size_x, size_y
	noInvFT = not os.path.exists(filename+"_ReconstructedLattice_Re.csv") or not os.path.exists(filename+"_ReconstructedLattice_Im.csv")

	if noInvFT:
		return findInvFT(filename,rec_data,latticeSpacing,size_x,size_y)

	data_re = []
	data_im = []
	with open(filename+"_ReconstructedLattice_Re.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			data_re.append(row)

	with open(filename+"_ReconstructedLattice_Im.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			data_im.append(row)

	return data_re, data_im

# Return an image that shows the Fourier transform in k-space
def IFTImage(braggFilteredData,return_type):
	global size_x, size_y, scale

	createReconstructedImage(filename,braggFilteredData,size_x,size_y,return_type)

	if return_type == "Re":
		img = NewImage.open(filename+"_BraggFiltered_Re.gif")
	elif return_type == "Im":
		img = NewImage.open(filename+"_BraggFiltered_Im.gif")
	elif return_type == "Mag":
		img = NewImage.open(filename+"_BraggFiltered_Mag.gif")
	img = img.resize((int(size_x*scale),int(size_y*scale)),NewImage.ANTIALIAS)
	img.save(filename+"_BraggFilteredScaled.gif","gif")
	img = Image(Point(int(size_x*scale/2.0),int(size_y*scale/2.0)), filename+"_BraggFilteredScaled.gif")

	return img

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

	# Determine data value corresponding to white in CDW image
	global max_point, min_point
	max_point = max(max(line) for line in data)
	min_point = min(min(line) for line in data)

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
	# print(latticeSpacing)

	# Generate and show Fourier transform
	rec_data = retFT(peaks,latticeSpacing)
	global FTWidth, FTHeight
	global FTImgWidth, FTImgHeight
	FTWidth = len(rec_data[0][0])
	FTHeight = len(rec_data[0])
	imgFT = FTImage(rec_data,num_peaks,"Mag")
	winFT = GraphWin('CDW Data - Fourier Transform', FTImgWidth, FTImgHeight)
	imgFT.draw(winFT)

	maskCenter = winFT.getMouse()
	kx_center = 8.0*pi*(maskCenter.getX()-FTImgWidth/2.0)/FTImgWidth/latticeSpacing
	ky_center = 8.0*pi*(maskCenter.getY()-FTImgHeight/2.0)/FTImgHeight/latticeSpacing
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
	radius = 8.0*pi*maskRadius/FTImgWidth/latticeSpacing
	mark = Circle(maskCenter,int(maskRadius))
	mark.setOutline(color_rgb(255,255,0))
	mark.draw(winFT)

	winFT.getMouse()
	winFT.close()

	# Apply mask on the Fourier Transform
	mask = stepFTMask
	rec_data_masked = mask(rec_data, latticeSpacing, kx_center, ky_center, radius)

	# Generate and show Fourier transform
	braggFilteredData = retInvFT(rec_data_masked,latticeSpacing)
	imgIFT = IFTImage(braggFilteredData,"Mag")
	winIFT = GraphWin('CDW Data - Inverse Fourier Transform', int(size_x*scale), int(size_y*scale))
	imgIFT.draw(winIFT)
	winIFT.getMouse()
	winIFT.close()

	win.getMouse()
	win.close()

main()