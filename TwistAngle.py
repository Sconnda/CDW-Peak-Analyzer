import os
import sys

from math import log10, floor
import numpy as np
import csv

from graphics import *
from PIL import Image as NewImage

from modules.hycht import *

# Declaring global variables_____________________
filename = "CDW_Data"
# Variables related to peaks
peaks = []
# Window and size of window
win = 0
size_x = 512
size_y = 512
scale = 1
# Adjust value of a to reflect lattice spacing in pixels
# Lattice spacing of graphene: 2.464 Angstroms
a = 34.5

pi = np.pi
sqrt = np.sqrt
asin = np.arcsin
floor = np.floor

# Round to three significant figures
def round_to_3(x):
	return round(x,2-int(floor(log10(abs(x)))))

# Identify numpy array elements that are outliers
def zero_one_matrix(A):
	B = floor(A)-0.5
	B /= abs(B)
	B += 1
	B /= 2
	return B

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

def fieldImage(fieldData,folder,extension):
	global filename, size_x, size_y, scale

	if not os.path.isdir(folder):
		os.mkdir(folder)

	createFieldImage(filename,fieldData,folder,extension,size_x,size_y)

	img = NewImage.open(folder+"/"+filename+"_"+folder+"_"+extension+".gif")
	img = img.resize((int((size_x)*scale),int(size_y*scale)),NewImage.ANTIALIAS)
	img.save(folder+"/"+filename+"_"+folder+"_"+extension+"_Scaled.gif","gif")
	img = Image(Point(int((size_x)*scale/2.0),int(size_y*scale/2.0)), folder+"/"+filename+"_"+folder+"_"+extension+"_Scaled.gif")

	winField = GraphWin(folder+": "+extension, int((size_x)*scale), int(size_y*scale))
	img.draw(winField)
	winField.getMouse()
	winField.close()

	return img

def twistAngleData(folder,gR_Data):
	global a
	b = 2/a/sqrt(3)
	zero_one_gR = zero_one_matrix(gR_Data/b/2)
	thetaData = 2*asin(gR_Data/b/2*(1-zero_one_gR)+zero_one_gR)*180/pi
	print("Average twist angle: "+str(round_to_3(thetaData.mean()))+"\N{DEGREE SIGN}")
	print("Twist angle standard deviation: "+str(round_to_3(np.std(thetaData)))+"\N{DEGREE SIGN}")
	thetaData = thetaData.tolist()

	with open(filename+"_"+folder+"_TwistAngle.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in thetaData:
			wr.writerow(row)

	return thetaData

def main():
	global filename, scale, a
	if len(sys.argv) > 1:
		filename = sys.argv[1]
	if len(sys.argv) > 2:
		scale = float(sys.argv[2])
	if len(sys.argv) > 3:
		a = float(sys.argv[3])
	if not os.path.isdir(filename):
		filename = "TestCDW_512px"
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

	# Load g(r)
	gR_Data1 = []
	with open(filename+"_BraggFiltered1_g(r).csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			gR_Data1.append(row)
	gR_Data1 = np.array(gR_Data1)

	gR_Data2 = []
	with open(filename+"_BraggFiltered2_g(r).csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			gR_Data2.append(row)
	gR_Data2 = np.array(gR_Data2)

	# Generate twist angle
	print("FIRST RECIPROCAL LATTICE VECTOR")
	theta1 = twistAngleData("BraggFiltered1",gR_Data1)
	print("\n")
	print("SECOND RECIPROCAL LATTICE VECTOR")
	theta2 = twistAngleData("BraggFiltered2",gR_Data2)
	print("\n")
	delTheta = abs(np.array(theta2)-np.array(theta1)).tolist()
	with open(filename+"_FieldImages_TwistAngleDifference.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in delTheta:
			wr.writerow(row)

	img_theta1 = fieldImage(theta1,"BraggFiltered1","TwistAngle")
	img_theta2 = fieldImage(theta2,"BraggFiltered2","TwistAngle")
	img_delTheta = fieldImage(delTheta,"FieldImages","TwistAngleDifference")

	win.getMouse()
	win.close()

main()