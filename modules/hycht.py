import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import math
import random
import csv

from graphics import *
from PIL import Image as NewImage
from PIL import ImageDraw as NewImageDraw
from PIL import ImageFont as NewImageFont

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

	Y,X = np.mgrid[0:size_y:1, 0:size_x:1]

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
def findDisplacementData(filename,extension1,extension2,g1,g2):
        

        PhaseData1 = []
        with open(filename+"_"+extension1+"_Phase.csv",newline='') as file:
            reader = csv.reader(file,delimiter=',',quotechar='|')
            for row in reader:
                for i,x in enumerate(row):
                    row[i] = float(x)
                PhaseData1.append(row)
        
        PhaseData2 = []
        with open(filename+"_"+extension2+"_Phase.csv",newline='') as file:
            reader = csv.reader(file,delimiter=',',quotechar='|')
            for row in reader:
                for i,x in enumerate(row):
                    row[i] = float(x)
                PhaseData2.append(row)



        #dummy g1 and g2 for now
        #g1 = np.array([[15],[13]])
        #g2 = np.array([[7],[4]])

        #make a g matrix
        gMatrix = np.vstack((g1,g2))
        
        #transpose and invert it according to the paper to get the a matrix
        gMatrixTranspose = np.transpose(gMatrix)
        aMatrix = np.linalg.inv(gMatrixTranspose)

        #seperate a1 and a2 vectors
        a1 = np.array([[aMatrix[0,0]],[aMatrix[1,0]]])
        a2 = np.array([[aMatrix[0,1]],[aMatrix[1,1]]])

        

        #calculate the displacement field x and y components 
        displacementFieldX = a1[0]*PhaseData1+a2[0]*PhaseData2
        displacementFieldY = a1[1]*PhaseData1+a2[1]*PhaseData2

        #save them as a csv file, there is only 1 extension because the displacement field is a property of the lattice not the g vector, I just picked extension1 arbitrarily
        with open(filename+"_"+extension1+"_DisplacementFieldX.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in displacementFieldX:
                wr.writerow(row)

        with open(filename+"_"+extension1+"_DisplacementFieldY.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in displacementFieldY:
                wr.writerow(row)

        return displacementFieldX, displacementFieldY


def findDistortionField(filename,extension1,extension2,g1,g2,size_x,size_y)
        
        #begining of this is the same as findDisplacementField, this is done to avoid compounding errors by differentiating the phase instead of the displacement field
        PhaseData1 = []
        with open(filename+"_"+extension1+"_Phase.csv",newline='') as file:
            reader = csv.reader(file,delimiter=',',quotechar='|')
            for row in reader:
                for i,x in enumerate(row):
                    row[i] = float(x)
                PhaseData1.append(row)


        PhaseData2 = []
        with open(filename+"_"+extension2+"_Phase.csv",newline='') as file:
            reader = csv.reader(file,delimiter=',',quotechar='|')
            for row in reader:
                for i,x in enumerate(row):
                    row[i] = float(x)
                PhaseData2.append(row)

        #make a g matrix
        gMatrix = np.vstack((g1,g2))
        
        #transpose and invert it according to the paper to get the a matrix
        gMatrixTranspose = np.transpose(gMatrix)
        aMatrix = np.linalg.inv(gMatrixTranspose)

        #seperate a1 and a2 vectors, 1st term is x second is y
        a1 = np.array([[aMatrix[0,0]],[aMatrix[1,0]]])
        a2 = np.array([[aMatrix[0,1]],[aMatrix[1,1]]])


        #create 0 matricies to hold the derivatives
        
        dPg1dx = [[0 for x in range(size_x)] for y in range(size_y)]
        dPg1dy = [[0 for x in range(size_x)] for y in range(size_y)]
        dPg2dx = [[0 for x in range(size_x)] for y in range(size_y)]
        dPg2dy = [[0 for x in range(size_x)] for y in range(size_y)]


        #basically a modified gradient to match appendix D of hycht paper one for g1 phase data and one for g2 phasedata
        
        #g1 stuff
        for y in range(size_y):
            for x in range(size_x):
                 
                dpdx = 0
                dpdy = 0

                if x == 0:
                    dpdx = sin(PhaseData1[y][1]-PhaseData1[y][0])
                elif x == size_x-1:
                    dpdx = -sin(PhaseData1[y][x-1]-PhaseData1[y][x])
                else:
                    dpdx = (sin(PhaseData1[y][x+1]-PhaseData1[y][x])-sin(PhaseData1[y][x-1]-PhaseData1[y][x]))/2

                if y == 0:
                    dpdy = sin(PhaseData1[1][x]-PhaseData1[0][x])
                elif y == size_y-1:
                    dpdy = -sin(PhaseData1[y-1][x]-PhaseData1[y][x])
                else:
                    dpdy = (sin(PhaseData1[y+1][x]-PhaseData1[y][x])-sin(PhaseData1[y-1][x]-PhaseData1[y][x]))/2

                dPg1dx[y][x] = dpdx
                dPg1dy[y][x] = dpdy

                if x<10 and y<10:
                    print("G1: dpdx is " + dpdx + " dpdy is " + dpdy)




        #g2 stuff
        for y in range(size_y):
            for x in range(size_x):
                 
                dpdx = 0
                dpdy = 0

                if x == 0:
                    dpdx = sin(PhaseData2[y][1]-PhaseData2[y][0])
                elif x == size_x-1:
                    dpdx = -sin(PhaseData2[y][x-1]-PhaseData2[y][x])
                else:
                    dpdx = (sin(PhaseData2[y][x+1]-PhaseData2[y][x])-sin(PhaseData2[y][x-1]-PhaseData2[y][x]))/2

                if y == 0:
                    dpdy = sin(PhaseData2[1][x]-PhaseData2[0][x])
                elif y == size_y-1:
                    dpdy = -sin(PhaseData2[y-1][x]-PhaseData2[y][x])
                else:
                    dpdy = (sin(PhaseData2[y+1][x]-PhaseData2[y][x])-sin(PhaseData2[y-1][x]-PhaseData2[y][x]))/2

                dPg2dx[y][x] = dpdx
                dPg2dy[y][x] = dpdy

                if x<10 and y<10:
                    print("G2: dpdx is " + dpdx + " dpdy is " + dpdy)


        #make a zero matrix for the distortion matrix
        distortionField = np.zeros((2,2,size_y,size_x))
        
        #turn the phase derivatives into a matrix to be dotted with the amatrix
        phaseDerivativeMatrix = np.array([[dPg1dx,dPg1dy],[dPg2dx,dPg2dy]])

        
        #calculate the distortion field and return it as a numpy array of shape (2,2,size_y,size_x)
        for y in range(size_y):
            for x in range(size_x):
                distortionField[:,:,y,x] = -(1/(2*pi))*(aMatrix.dot(phaseDerivativeMatrix[:,:,y,x]))
        

        #this bit should save still working on it 
        with open(filename+"_"+extension1+"_DisplacementFieldX.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in displacementFieldX:
                wr.writerow(row)
        

        return distortionField



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
	if extension == "Phase":
		min_point = -2*pi
		max_point = 2*pi
	else:
		min_point = min(min(field))
		max_point = max(max(field))

	# Draw data
	img = NewImage.new("RGB", (size_x, size_y))
	putpixel = img.putpixel

	for y in range(size_y):
		for x in range(size_x):
			color = int((field[y][x]-min_point)/(max_point-min_point)*255)
			putpixel((x,y),(color,color,color))

	img.save(folder+"/"+filename+"_"+folder+"_"+extension+".gif",'gif')

	return img
