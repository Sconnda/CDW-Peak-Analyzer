import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import math
import random
import csv

from graphics import *
from PIL import Image as NewImage

pi = np.pi
sqrt = np.sqrt
sin = np.sin
cos = np.cos
atan2 = np.arctan2
exp = np.exp

# SKANDA

def findPhaseData(filename,size_x,size_y,G,extension):
	G_x = G[0]
	G_y = G[1]
	rawPhaseData = []
	with open(filename+"_"+extension+"_RawPhase.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			rawPhaseData.append(row)
	rawPhaseData = np.array(rawPhaseData)

	X,Y = np.mgrid[0:size_x:1, 0:size_y:1]

	# Y = np.flip(Y,1)

	phaseData = (rawPhaseData-2*pi*G_x*X-2*pi*G_y*Y+pi)%(2*pi)-pi

	with open(filename+"_"+extension+"_Phase.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in phaseData:
			wr.writerow(row)

	return phaseData

def gradient(data,x,y,size_x,size_y):
	grad_x = 0
	grad_y = 0

	if x == 0:
		grad_x = data[y][1]-data[y][0]
	elif x == size_x-1:
		grad_x = data[y][x]-data[y][x-1]
	else:
		grad_x = (data[y][x+1]-data[y][x-1])/2

	if y == 0:
		grad_y = data[1][x]-data[0][x]
	elif y == size_y-1:
		grad_y = data[y][x]-data[y-1][x]
	else:
		grad_y = (data[y+1][x]-data[y-1][x])/2

	grad = (grad_x,grad_y)
	return grad

def find_gR_Data(filename,size_x,size_y,extension):
	rawPhaseData = []
	with open(filename+"_"+extension+"_RawPhase.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			rawPhaseData.append(row)
	rawPhaseData = np.array(rawPhaseData)

	gR_x = [[0 for x in range(size_x)] for y in range(size_y)]
	gR_y = [[0 for x in range(size_x)] for y in range(size_y)]
	for y in range(size_y):
		for x in range(size_x):
			grad = gradient(rawPhaseData,x,y,size_x,size_y)
			gR_x[y][x] = 2*pi*grad[0]
			gR_y[y][x] = 2*pi*grad[1]

	with open(filename+"_"+extension+"_g(r)_x.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in gR_x:
			wr.writerow(row)

	with open(filename+"_"+extension+"_g(r)_y.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in gR_y:
			wr.writerow(row)

	return gR_x,gR_y

# DAN

# Save displacement image data as a size_x by size_y csv file and return an array including the data. One of the inputs should be the phase image, I'll work on that
def findDisplacementData(filename,phaseData,size_x,size_y,g1,g2):
        
        #dummy g1 and g2 for now
        g1 = np.array([[15],[13]])
        g2 = np.array([[7],[4]])

        #make a g matrix
        gMatrix = np.concatenate((g1,g2),axis=1)
        
        #transpose and invert it according to the paper to get the a matrix
        gMatrixTranspose = np.transpose(gMatrix)
        aMatrix = np.linalg.inv(gMatrixTranspose)

        #seperate a1 and a2 vectors
        a1 = np.array([[aMatrix[0,0]],[aMatrix[1,0]]])
        a2 = np.array([[aMatrix[0,1]],[aMatrix[1,1]]])

        

	return 0

# Create, save, and return an image for any real field in a graphics window
def createFieldImage(filename,fieldData,folder,extension,size_x,size_y):
	field = []
	with open(filename+"_"+folder+"_"+extension+".csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			field.append(row)

	# Find extrema
	min_point = -pi
	max_point = pi

	# Draw data
	img = NewImage.new("RGB", (size_x, size_y))
	putpixel = img.putpixel
	for y in range(size_y):
		for x in range(size_x):
			color = int((field[y][x]-min_point)/(max_point-min_point)*255)
			putpixel((x,y),(color,color,color))

	img.save(folder+"/"+filename+"_"+folder+"_"+extension+".gif",'gif')

	return img
