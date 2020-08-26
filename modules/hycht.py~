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

def mag(re,im):
    return sqrt(re**2+im**2)

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

def phaseGradient(data,x,y,size_x,size_y):
	grad_x = 0
	grad_y = 0

	if x == 0:
		grad_x = sin(data[y][1]-data[y][0])
	elif x == size_x-1:
		grad_x = sin(data[y][x]-data[y][x-1])
	else:
		grad_x = (sin(data[y][x+1]-data[y][x])+sin(data[y][x]-data[y][x-1]))/2

	if y == 0:
		grad_y = sin(data[1][x]-data[0][x])
	elif y == size_y-1:
		grad_y = sin(data[y][x]-data[y-1][x])
	else:
		grad_y = (sin(data[y+1][x]-data[y][x])+sin(data[y][x]-data[y-1][x]))/2

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

    mag_vect = np.vectorize(mag)
    gR_x = [[0 for x in range(size_x)] for y in range(size_y)]
    gR_y = [[0 for x in range(size_x)] for y in range(size_y)]

    for y in range(size_y):
        for x in range(size_x):
            grad = phaseGradient(rawPhaseData,x,y,size_x,size_y)
            gR_x[y][x] = grad[0]/(2*pi)
            gR_y[y][x] = grad[1]/(2*pi)

    gR = mag_vect(np.array(gR_x), np.array(gR_y))

    with open(filename+"_"+extension+"_g(r)_x.csv", 'w',newline='') as f:
        wr = csv.writer(f)
        for row in gR_x:
            wr.writerow(row)

    with open(filename+"_"+extension+"_g(r)_y.csv", 'w',newline='') as f:
        wr = csv.writer(f)
        for row in gR_y:
            wr.writerow(row)

    with open(filename+"_"+extension+"_g(r).csv", 'w',newline='') as f:
        wr = csv.writer(f)
        for row in gR:
            wr.writerow(row)

    return gR_x,gR_y

# DAN

# Save displacement image data as a size_x by size_y csv file and return an array including the data. One of the inputs should be the phase image, I'll work on that
def findDisplacementData(filename,g1,g2):

        PhaseData1 = []
        with open(filename+"_BraggFiltered1_Phase.csv",newline='') as file:
            reader = csv.reader(file,delimiter=',',quotechar='|')
            for row in reader:
                for i,x in enumerate(row):
                    row[i] = float(x)
                PhaseData1.append(row)
        
        PhaseData2 = []
        with open(filename+"_BraggFiltered2_Phase.csv",newline='') as file:
            reader = csv.reader(file,delimiter=',',quotechar='|')
            for row in reader:
                for i,x in enumerate(row):
                    row[i] = float(x)
                PhaseData2.append(row)



        #dummy g1 and g2 for now
        #g1 = np.array([[15],[13]])
        #g2 = np.array([[7],[4]])

        #make a g matrix
        gTMatrix = np.vstack((g1,g2))
        
        #transpose and invert it according to the paper to get the a matrix
        aMatrix = np.linalg.inv(gTMatrix)

        #seperate a1 and a2 vectors
        a1 = np.array([[aMatrix[0,0]],[aMatrix[1,0]]])
        a2 = np.array([[aMatrix[0,1]],[aMatrix[1,1]]])

        

        #calculate the displacement field x and y components 
        displacementFieldX = -(1/(2*pi))*(a1[0]*PhaseData1+a2[0]*PhaseData2)
        displacementFieldY = -(1/(2*pi))*(a1[1]*PhaseData1+a2[1]*PhaseData2)

        #save them as a csv file, there is only 1 extension because the displacement field is a property of the lattice not the g vector, I just picked extension1 arbitrarily
        with open(filename+"_FieldImages_DisplacementFieldX.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in displacementFieldX:
                wr.writerow(row)

        with open(filename+"_FieldImages_DisplacementFieldY.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in displacementFieldY:
                wr.writerow(row)

        return displacementFieldX, displacementFieldY


def findDistortionField(filename,extension1,extension2,g1,g2,size_x,size_y):
 
        # Begining of this is the same as findDisplacementField, this is done to avoid compounding errors by differentiating the phase instead of the displacement field
        PhaseData1 = []
        with open(filename+"_BraggFiltered1_Phase.csv",newline='') as file:
            reader = csv.reader(file,delimiter=',',quotechar='|')
            for row in reader:
                for i,x in enumerate(row):
                    row[i] = float(x)
                PhaseData1.append(row)


        PhaseData2 = []
        with open(filename+"_BraggFiltered2_Phase.csv",newline='') as file:
            reader = csv.reader(file,delimiter=',',quotechar='|')
            for row in reader:
                for i,x in enumerate(row):
                    row[i] = float(x)
                PhaseData2.append(row)

        #make a g matrix
        gTMatrix = np.vstack((g1,g2))
        
        # Invert it according to the paper to get the a matrix
        #   Note: vstack arranges the matrix so that it is already transposed as in the paper
        aMatrix = np.linalg.inv(gTMatrix)

        #create 0 matricies to hold the derivatives        
        dPg1dx = [[0 for x in range(size_x)] for y in range(size_y)]
        dPg1dy = [[0 for x in range(size_x)] for y in range(size_y)]
        dPg2dx = [[0 for x in range(size_x)] for y in range(size_y)]
        dPg2dy = [[0 for x in range(size_x)] for y in range(size_y)]


        #basically a modified gradient to match appendix D of hycht paper one for g1 phase data and one for g2 phasedata
        
        for y in range(size_y):
            for x in range(size_x):
                #g1 stuff
                grad1 = phaseGradient(PhaseData1,x,y,size_x,size_y)
                dPg1dx[y][x] = grad1[0]
                dPg1dy[y][x] = grad1[1]

                #g2 stuff
                grad2 = phaseGradient(PhaseData2,x,y,size_x,size_y)
                dPg2dx[y][x] = grad2[0]
                dPg2dy[y][x] = grad2[1]

        #make a zero matrix for the distortion matrix
        distortionField = np.zeros((2,2,size_y,size_x))
        
        #turn the phase derivatives into a matrix to be dotted with the amatrix
        phaseDerivativeMatrix = np.array([[dPg1dx,dPg1dy],[dPg2dx,dPg2dy]])
        
        #calculate the distortion field and return it as a numpy array of shape (2,2,size_y,size_x)
        for y in range(size_y):
            for x in range(size_x):
                distortionField[:,:,y,x] = -(1/(2*pi))*(aMatrix.dot(phaseDerivativeMatrix[:,:,y,x]))

        return distortionField


#returns the strain and rotation matricies for every point as numpy arrays of shape (2,2,size_y,size_x) 
def findStrainAndRotation(filename,distortionField,size_x,size_y):
        
        #make the an array to store the transpose 
        transposeDistortion = np.zeros((2,2,size_y,size_x))

        #make the transpose
        transposeDistortion[0,0,:,:] = distortionField[0,0,:,:]
        transposeDistortion[1,1,:,:] = distortionField[1,1,:,:]
        transposeDistortion[0,1,:,:] = distortionField[1,0,:,:]
        transposeDistortion[1,0,:,:] = distortionField[0,1,:,:]

        #the strain field at every r (for small distortions) is 1/2(distortion + transposeDistortion)
        strain = (transposeDistortion+distortionField)/2

        #the rotation matrix at every r (for small distortions) is 1/2(distortion - transposeDistortion)
        rotation = (distortionField-transposeDistortion)/2

        #create size_x by size_y arrays to store the components of the strain and rotation fields
        strainXX = strain[0,0,:,:]
        strainXY = strain[0,1,:,:]
        strainYX = strain[1,0,:,:]
        strainYY = strain[1,1,:,:]

        #rotation xx and yy should be all zero but including them to make sure nothings going wrong

        rotationXX = rotation[0,0,:,:]
        rotationXY = rotation[0,1,:,:]
        rotationYX = rotation[1,0,:,:]
        rotationYY = rotation[1,1,:,:]
        
        
        #now we push them all out as csv files

        #strain ones
        with open(filename+"_FieldImages_StrainXX.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in strainXX:
                wr.writerow(row)


        with open(filename+"_FieldImages_StrainXY.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in strainXY:
                wr.writerow(row)

        with open(filename+"_FieldImages_StrainYX.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in strainYX:
                wr.writerow(row)

        with open(filename+"_FieldImages_StrainYY.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in strainYY:
                wr.writerow(row)

        #rotation ones

        with open(filename+"_FieldImages_RotationXX.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in rotationXX:
                wr.writerow(row)

        with open(filename+"_FieldImages_RotationXY.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in rotationXY:
                wr.writerow(row)

        with open(filename+"_FieldImages_RotationYX.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in rotationYX:
                wr.writerow(row)

        with open(filename+"_FieldImages_RotationYY.csv", 'w',newline='') as f:
            wr = csv.writer(f)
            for row in rotationYY:
                wr.writerow(row)

        return strain, rotation


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
        min_point = -pi
        max_point = pi
    elif extension == "TwistAngle" or extension == "TwistAngleDifference":
        field_np = np.array(field)
        min_point = 0
        max_point = field_np.mean()+np.std(field_np)
    elif extension in ["StrainXX","StrainYX","StrainXY","StrainYY","RotationXX","RotationYX","RotationXY","RotationYY"]:
        min_point = -.1
        max_point = .1
    elif extension in ["g(r)_x","g(r)_y"]:
        min_point = min(min(field))
        max_point = max(max(field))
        min_point = min(min_point,-max_point)
        max_point = max(-min_point,max_point)
    else:
        min_point = min(min(field))
        max_point = max(max(field))
    print("Black: "+str(min_point))
    print("White: "+str(max_point))

    # Draw data
    img = NewImage.new("RGB", (size_x, size_y))
    putpixel = img.putpixel
    

    if max_point!=min_point:
        for y in range(size_y):
            for x in range(size_x):
                color = int((field[y][x]-min_point)/(max_point-min_point)*255)
                putpixel((x,y),(color,color,color))
    elif max_point==min_point:
        for y in range(size_y):
            for x in range(size_x):
                color = int(field[y][x]-min_point)
                putpixel((x,y),(color,color,color))


    img.save(folder+"/"+filename+"_"+folder+"_"+extension+".gif",'gif')

    return img
