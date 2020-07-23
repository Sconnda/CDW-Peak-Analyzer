import os
import sys

import numpy as np
import csv

from graphics import *
from PIL import Image as NewImage

from modules.bondFunctions import *
from modules.reciprocalLattice import *
from modules.hycht import *

# Declaring global variables_____________________
filename = "CDW_Data"
# Variables related to peaks
peaks = []
# latticeSpacing = 0
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

# Return two masked FT images
def selectG(peaks,num_peaks):
	global filename, size_x, size_y
	global FTWidth, FTHeight
	global FTImgWidth, FTImgHeight

	rec_data = findFT(filename,peaks,size_x,size_y,FTWidth,FTHeight)
	imgFT = FTImage(rec_data,num_peaks,"Mag","FT")
	winFT = GraphWin('Fourier Transform', FTImgWidth, FTImgHeight)
	imgFT.draw(winFT)

	maskCenter1 = winFT.getMouse()
	G1_x = 4*pi*(maskCenter1.getX()-FTImgWidth/2.0)/FTImgWidth/sqrt(3)
	G1_y = 4*pi*(maskCenter1.getY()-FTImgHeight/2.0)/FTImgHeight/sqrt(3)
	# print(latticeSpacing)
	# print(4*pi/sqrt(3)/sqrt(G1_x**2+G1_y**2))
	graphCenter1 = Point(int(FTImgWidth/2.0),int(FTImgHeight/2.0))
	ln1 = Line(graphCenter1,maskCenter1)
	ln1.setOutline(color_rgb(0,128,255))
	ln1.draw(winFT)
	mark1 = Circle(maskCenter1,3)
	mark1.setFill(color_rgb(0,128,255))
	mark1.setOutline(color_rgb(0,128,255))
	mark1.draw(winFT)

	maskRadiusEndpoint1 = winFT.getMouse()
	maskRadius1 = sqrt((maskRadiusEndpoint1.getX()-maskCenter1.getX())**2+(maskRadiusEndpoint1.getY()-maskCenter1.getY())**2)
	radius1 = 4*pi*maskRadius1/sqrt(FTImgWidth**2+FTImgHeight**2)/sqrt(3)
	mark1 = Circle(maskCenter1,int(maskRadius1))
	mark1.setOutline(color_rgb(255,255,0))
	mark1.draw(winFT)

	maskCenter2 = winFT.getMouse()
	G2_x = 4*pi*(maskCenter2.getX()-FTImgWidth/2.0)/FTImgWidth/sqrt(3)
	G2_y = 4*pi*(maskCenter2.getY()-FTImgHeight/2.0)/FTImgHeight/sqrt(3)
	# print(latticeSpacing)
	# print(4*pi/sqrt(3)/sqrt(G2_x**2+G2_y**2))
	graphCenter2 = Point(int(FTImgWidth/2.0),int(FTImgHeight/2.0))
	ln2 = Line(graphCenter2,maskCenter2)
	ln2.setOutline(color_rgb(0,128,255))
	ln2.draw(winFT)
	mark2 = Circle(maskCenter2,3)
	mark2.setFill(color_rgb(0,128,255))
	mark2.setOutline(color_rgb(0,128,255))
	mark2.draw(winFT)

	maskRadiusEndpoint2 = winFT.getMouse()
	maskRadius2 = sqrt((maskRadiusEndpoint2.getX()-maskCenter2.getX())**2+(maskRadiusEndpoint2.getY()-maskCenter2.getY())**2)
	radius2 = 4*pi*maskRadius2/sqrt(FTImgWidth**2+FTImgHeight**2)/sqrt(3)
	mark2 = Circle(maskCenter2,int(maskRadius2))
	mark2.setOutline(color_rgb(255,255,0))
	mark2.draw(winFT)

	winFT.getMouse()
	winFT.close()

	# Apply mask on the Fourier Transform
	mask = stepFTMask
	rec_data_masked1 = mask(rec_data, G1_x, G1_y, radius1)
	imgFT_masked1 = FTImage(rec_data_masked1,num_peaks,"Mag","FT_masked1")
	rec_data_masked2 = mask(rec_data, G2_x, G2_y, radius2)
	imgFT_masked2 = FTImage(rec_data_masked2,num_peaks,"Mag","FT_masked2")

	return (rec_data_masked1,rec_data_masked2), ((G1_x,G1_y),(G2_x,G2_y))

# Return an image that shows the inverse Fourier transform in real space
def IFTImage(braggFilteredData,return_type,extension):
	global size_x, size_y, scale

	createReconstructedImage(filename,braggFilteredData,size_x,size_y,return_type,extension)

	if return_type == "Re":
		img = NewImage.open(filename+"_"+extension+"_Re.gif")
	elif return_type == "Im":
		img = NewImage.open(filename+"_"+extension+"_Im.gif")
	elif return_type == "Phase":
		img = NewImage.open(filename+"_"+extension+"_RawPhase.gif")
	else: # return_type == "Mag" by default
		img = NewImage.open(filename+"_"+extension+"_Mag.gif")
	img = img.resize((int(size_x*scale),int(size_y*scale)),NewImage.ANTIALIAS)
	img.save(filename+"_"+extension+"_Scaled.gif","gif")
	img = Image(Point(int(size_x*scale/2.0),int(size_y*scale/2.0)), filename+"_"+extension+"_Scaled.gif")

	return img

def fieldImage(fieldData,folder,extension):
	global filename, size_x, size_y, scale

	if not os.path.isdir(folder):
		os.mkdir(folder)

	createFieldImage(filename,fieldData,folder,extension,size_x,size_y)

	img = NewImage.open(folder+"/"+filename+"_"+folder+"_"+extension+".gif")
	img = img.resize((int(size_x*scale),int(size_y*scale)),NewImage.ANTIALIAS)
	img.save(folder+"/"+filename+"_"+folder+"_"+extension+"_Scaled.gif","gif")
	img = Image(Point(int(size_x*scale/2.0),int(size_y*scale/2.0)), folder+"/"+filename+"_"+folder+"_"+extension+"_Scaled.gif")

	winField = GraphWin(folder+": "+extension, int(size_x*scale), int(size_y*scale))
	img.draw(winField)
	winField.getMouse()
	winField.close()

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

	# Show original image
	global win
	win = GraphWin('Original lattice', int(size_x*scale), int(size_y*scale))
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

	# # Calculate the lattice spacing (useful in determing the scale for many calculations)
	# jMatrix = retJMatrix(peaks,size_x,size_y)
	# bondMatrix = retTriangulation(jMatrix,num_peaks,size_x,size_y)

	# global latticeSpacing
	# latticeSpacing = 0
	# total_bonds = 0
	# for i,peak in enumerate(peaks):
	# 	x,y = peak
	# 	for j,bonded in enumerate(bondMatrix[i]):
	# 		if bonded == 1:
	# 			x2,y2 = peaks[j]
	# 			total_bonds += 1
	# 			latticeSpacing += sqrt((x2-x)**2 + (y2-y)**2)
	# latticeSpacing /= total_bonds

	rec_data_masked,G = selectG(peaks,num_peaks)
	rec_data_masked1, rec_data_masked2 = rec_data_masked
	G1, G2 = G

	# Generate and show Bragg filtered data for both values of g
	braggFilteredData1 = findInvFT(filename,rec_data_masked1,"BraggFiltered1",size_x,size_y)
	imgIFT1 = IFTImage(braggFilteredData1,braggFilteredData_return_type,"BraggFiltered1")
	winIFT1 = GraphWin('Bragg Filtered Data: g1', int(size_x*scale), int(size_y*scale))
	imgIFT1.draw(winIFT1)
	winIFT1.getMouse()
	winIFT1.close()

	braggFilteredData2 = findInvFT(filename,rec_data_masked2,"BraggFiltered2",size_x,size_y)
	imgIFT2 = IFTImage(braggFilteredData2,braggFilteredData_return_type,"BraggFiltered2")
	winIFT2 = GraphWin('Bragg Filtered Data: g2', int(size_x*scale), int(size_y*scale))
	imgIFT2.draw(winIFT2)
	winIFT2.getMouse()
	winIFT2.close()

	# Generate phase image
	phaseData1 = findPhaseData(filename,size_x,size_y,G1,"BraggFiltered1")
	phaseData2 = findPhaseData(filename,size_x,size_y,G2,"BraggFiltered2")
	imgPhaseField1 = fieldImage(phaseData1,"BraggFiltered1","Phase")
	imgPhaseField2 = fieldImage(phaseData2,"BraggFiltered2","Phase")

	# Generate g(r)
	gR_Data1 = find_gR_Data(filename,size_x,size_y,"BraggFiltered1")
	gR_Data2 = find_gR_Data(filename,size_x,size_y,"BraggFiltered2")
	img_gR_Field1_x = fieldImage(gR_Data1[0],"BraggFiltered1","g(r)_x")
	img_gR_Field1_y = fieldImage(gR_Data1[1],"BraggFiltered1","g(r)_y")
	img_gR_Field2_x = fieldImage(gR_Data2[0],"BraggFiltered2","g(r)_x")
	img_gR_Field2_y = fieldImage(gR_Data2[1],"BraggFiltered2","g(r)_y")

	# # Generate displacement image
	# displacementData = findDisplacementData(filename,phaseData,size_x,size_y)
	# imgDisplacementField = fieldImage(displacementData,"DisplacementField1")
	# winDisplacementField = GraphWin('Displacement Field', int(size_x*scale), int(size_y*scale))
	# imgDisplacementField.draw(winDisplacementField)
	# winDisplacementField.getMouse()
	# winDisplacementField.close()

	win.getMouse()
	win.close()

main()