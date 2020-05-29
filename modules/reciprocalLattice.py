import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import math
import random
import csv

from graphics import *
from PIL import Image as NewImage

# Declaring global variables_____________________
filename = "CDW_Data"
# Cutoff radius for orientational order correlation function G6(r)
cutoff = 100
# Window and size of window
win = 0
size_x = 512
size_y = 512

pi = np.pi
sin = np.sin
cos = np.cos

def findFT(filename,peaks, latticeSpacing, width, height):
	print("Calculating Fourier Transform...")
	rec_data = [[0 for x in range(width)] for y in range(height)]
	pDone = 0
	for y in range(height):
		for x in range(width):
			G_x = 4*pi*(2*x/width-1)/latticeSpacing
			G_y = 4*pi*(2*y/height-1)/latticeSpacing
			n = [0,0]
			for peak in peaks:
				R_x,R_y = peak
				n[0] += cos(G_x*R_x+G_y*R_y)
				n[1] -= sin(G_x*R_x+G_y*R_y)
			rec_data[y][x] = n[0]
		if int(100*float(y+1)/height) > pDone:
			pDone = int(100*float(y+1)/height)
			print(str(pDone)+"%")

	with open(filename+"_ReciprocalLattice.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in rec_data:
			wr.writerow(row)

	return rec_data

def createRecDataImage(filename,rec_data,num_peaks):

	# Create window for drawing
	width = len(rec_data[0])
	height = len(rec_data)

	min_point = num_peaks
	max_point = -num_peaks
	ave = 0
	# Find extrema
	for y in range(height):
		for x in range(width):
			# Ignoring central peak, necessary when the peaks are very smeared at non-zero G
			if abs(x-width/2) < 20 and abs(y-height/2) < 20:
				continue
			val = rec_data[y][x]
			ave += val
			if val < min_point:
				min_point = val
			if val > max_point:
				max_point = val
	ave /= (width*height)

	# Draw data
	img = NewImage.new("RGB", (width, height))
	putpixel = img.putpixel
	for y in range(height):
		for x in range(width):
			# Ignoring central peak, necessary when the peaks are very smeared at non-zero G
			if abs(x-width/2) < 20 and abs(y-height/2) < 20:
				continue
			color = max(int((rec_data[y][x]-ave)/(max_point-ave)*255),0)
			putpixel((x,y),(color,color,color))

	img.save(filename+"_FT.gif",'gif')
	img = Image(Point(int(width/2),int(height/2)), filename+"_FT.gif")

	return img