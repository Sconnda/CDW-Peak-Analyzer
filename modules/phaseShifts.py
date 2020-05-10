import math
import numpy as np
import csv
import cv2 as cv

beta = 1

def prob_flip(C,rand):
	x = 2*rand-1
	return np.exp(beta*C*x)*beta*C/(2*np.sinh(beta*C))

def nudge(phaseShift,width,height):
	i = int(2*width*np.random.sample())
	j = int(2*height*np.random.sample())
	C = 0
	for di in (-1,0,1):
		for dj in (-1,0,1):
			if di==0 and dj==0:
				continue
			if i+di < 0 or i+di >= 2*width:
				continue
			if j+dj < 0 or j+dj >= 2*height:
				continue
			C += 2*phaseShift[i+di][j+dj][0]-1
	x = np.random.sample()
	p_x = prob_flip(C,x)
	flip = np.random.sample()
	for k in range(3):
		if flip < p_x:
			phaseShift[i][j][k] = x

	return phaseShift

def genPhaseShift(stress_factor,lattice_spacing,width,height):
	print("Generating phase shifts...")
	phaseShift = np.zeros([2*width,2*height,3],float)
	for i in range(2*width):
		for j in range(2*height):
			gray = np.random.sample()
			for k in range(3):
				phaseShift[i][j][k] = gray

	cv.imshow('phaseShift',phaseShift)
	key = cv.waitKey(0)
	cv.destroyAllWindows()

	phaseShift2 = np.zeros([2*width,2*height,3],float)
	for i in range(2*width):
		for j in range(2*height):
			for k in range(3):
				phaseShift2[i][j][k] = phaseShift[i][j][k]

	while True:
		phaseShift2 = nudge(phaseShift2,width,height)
		cv.imshow('phaseShift',phaseShift2)
		key = cv.waitKey(1) & 0xFF
		if key == ord('q'):
			break

	cv.destroyAllWindows()

def dummy():
	ave = [0,0,0]
	randomWalk = stress_factor*lattice_spacing
	pDone = 0
	for i in range(2*width):
		if 10*int(10*(i+1)/(2*width)) > pDone:
			pDone = 10*int(10*(i+1)/(2*width))
			print(str(pDone)+"%")
		for j in range(2*height):
			for k in [0,2]:
				if i == 0 and j == 0:
					phaseShift[i][j][k] = 0.5
					ave[k] += 0.5
				elif j == 0:
					shift = phaseShift[i-1][j][k]+randomWalk*(2*np.random.sample()-1)
					phaseShift[i][j][k] = shift
					ave[k] += shift
				elif i == 0:
					shift = phaseShift[i][j-1][k]+randomWalk*(2*np.random.sample()-1)
					phaseShift[i][j][k] = shift
					ave[k] += shift
				else:
					shift1 = phaseShift[i][j-1][k]
					shift2 = phaseShift[i-1][j][k]
					weight = np.random.sample()
					shift = weight*shift1+(1-weight)*shift2+randomWalk*2**0.5*(2*np.random.sample()-1)
					phaseShift[i][j][k] = shift
					ave[k] += shift

	for k in [0,2]:
		ave[k] /= 4*width**2+4*height**2
		for i in range(2*width):
			for j in range(2*height):
				phaseShift[i][j][k] -= ave[k]

	phaseShift2 = np.zeros([width,height,3],float)
	for i in range(width):
		for j in range(height):
			for k in [0,2]:
				phaseShift2[i][j][k] = phaseShift[i+int(width/2)][j+int(height/2)][k]/lattice_spacing+0.5

	cv.imshow('phaseShift',phaseShift2)
	key = cv.waitKey(0)
	cv.destroyAllWindows()

genPhaseShift(0.3,10,100,100)