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
sqrt = np.sqrt
sin = np.sin
cos = np.cos
atan2 = np.arctan2
exp = np.exp

def findFT(filename,peaks,size_x,size_y, latticeSpacing, width, height):
	# print("Calculating Fourier Transform...")
	data = [[0 for x in range(size_x)] for y in range(size_y)]
	for peak in peaks:
		x,y = peak
		x = min(int(round(x)),size_x-1)
		y = min(int(round(y)),size_y-1)
		data[y][x] = 1

	# data = np.loadtxt(filename+".txt")

	data_np = np.array(data)
	rec_data = np.fft.fft2(data_np)
	rec_data = np.fft.fftshift(rec_data)
	rec_data_re = rec_data.real.tolist()
	rec_data_im = rec_data.imag.tolist()

	# rec_data_re = [[0 for x in range(width)] for y in range(height)]
	# rec_data_im = [[0 for x in range(width)] for y in range(height)]
	# pDone = 0
	# for y in range(height):
	# 	for x in range(width):
	# 		# The image should have a width and height each corresponding to a range of (-4*pi/a, 4*pi/a)
	# 		k_x = 4*pi*(2.0*x/width-1)/latticeSpacing
	# 		k_y = 4*pi*(2.0*y/height-1)/latticeSpacing
	# 		n = [0,0]
	# 		for peak in peaks:
	# 			R_x,R_y = peak
	# 			n[0] += cos(k_x*R_x+k_y*R_y)
	# 			n[1] -= sin(k_x*R_x+k_y*R_y)
	# 		rec_data_re[y][x] = n[0]
	# 		rec_data_im[y][x] = n[1]
	# 	if round(0.1*int(1000*float(y+1)/height),2) > pDone:
	# 		pDone = round(0.1*int(1000*float(y+1)/height),2)
	# 		print(str(pDone)+"%")

	with open(filename+"_ReciprocalLattice.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in rec_data_re:
			wr.writerow(row)

	with open(filename+"_ReciprocalLattice_Im.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in rec_data_im:
			wr.writerow(row)

	return rec_data_re, rec_data_im

def stepFTMask(rec_data, latticeSpacing, kx_center, ky_center, radius):
	rec_data_re, rec_data_im = rec_data
	width = len(rec_data_re[0])
	height = len(rec_data_re)
	rec_data_filtered_re = [[0 for x in range(width)] for y in range(height)]
	rec_data_filtered_im = [[0 for x in range(width)] for y in range(height)]
	for i in range(height):
		ky = 8.0*pi*(i-height/2.0)/height/latticeSpacing
		for j in range(width):
			kx = 8.0*pi*(j-width/2.0)/width/latticeSpacing
			if (kx-kx_center)**2+(ky-ky_center)**2 < radius**2:
				rec_data_filtered_re[i][j] = rec_data_re[i][j]
				rec_data_filtered_im[i][j] = rec_data_im[i][j]
	return rec_data_filtered_re,rec_data_filtered_im

def gaussianFTMask(rec_data, latticeSpacing, kx_center, ky_center, radius):
	sigma = radius/sqrt(2)
	rec_data_re, rec_data_im = rec_data
	width = len(rec_data_re[0])
	height = len(rec_data_re)
	rec_data_filtered_re = [[0 for x in range(width)] for y in range(height)]
	rec_data_filtered_im = [[0 for x in range(width)] for y in range(height)]
	for i in range(height):
		ky = 8.0*pi*(i-height/2.0)/height/latticeSpacing
		for j in range(width):
			kx = 8.0*pi*(j-width/2.0)/width/latticeSpacing
			mask = exp(-((kx-kx_center)**2+(ky-ky_center)**2)/(2*sigma**2))
			rec_data_filtered_re[i][j] = mask*rec_data_re[i][j]
			rec_data_filtered_im[i][j] = mask*rec_data_im[i][j]
	return rec_data_filtered_re,rec_data_filtered_im

def highPassFilterMask(kxSc,kySc):
	return (1-exp(-pi*(kxSc**2+kySc**2)/2))

def globalFTFilter(rec_data):
	width = len(rec_data[0])
	height = len(rec_data)
	rec_data_filtered = rec_data
	for i in range(height):
		for j in range(width):
			kxSc = 2.0*(j-width/2.0)/width
			kySc = 2.0*(i-height/2.0)/height
			mask = highPassFilterMask
			rec_data_filtered[i][j] = mask(kxSc,kySc)*rec_data[i][j]
	return rec_data_filtered

def createRecDataImage(filename,rec_data,num_peaks,width,height,return_type):
	dataWidth = len(rec_data[0][0])
	dataHeight = len(rec_data[0])

	if return_type == "Re":
		rec_data = rec_data[0]
	elif return_type == "Im":
		rec_data = rec_data[1]
	elif return_type == "Mag":
		rec_data_mag = [[0 for x in range(dataWidth)] for y in range(dataHeight)]
		for i in range(dataHeight):
			for j in range(dataWidth):
				rec_data_mag[i][j] = sqrt(rec_data[0][i][j]**2+rec_data[1][i][j]**2)
		rec_data = rec_data_mag

	min_point = num_peaks
	max_point = -num_peaks
	ave = 0

	# # Pass through a high-pass filter to increase the brightness of the edges
	# rec_data = globalFTFilter(rec_data)

	# Find extrema
	for y in range(dataHeight):
		for x in range(dataWidth):
			# # Ignoring central peak, necessary when the peaks are very smeared at non-zero G
			# if abs(x-dataWidth/2) < 20 and abs(y-dataHeight/2) < 20:
			# 	continue

			val = rec_data[y][x]
			ave += val
			if val < min_point:
				min_point = val
			if val > max_point:
				max_point = val
	ave /= (dataWidth*dataHeight)

	# Draw data
	img = NewImage.new("RGB", (width, height))
	putpixel = img.putpixel
	for Y in range(height):
		for X in range(width):
			x = int(X*float(dataWidth)/width)
			y = int(Y*float(dataHeight)/height)
			# # Ignoring central peak, necessary when the peaks are very smeared at non-zero G
			# if abs(x-dataWidth/2) < 20 and abs(y-dataHeight/2) < 20:
			# 	continue
			color = min(int((rec_data[y][x]-min_point)/(max_point-min_point)*255),255)
			putpixel((X,Y),(color,color,color))

	if return_type == "Re":
		img.save(filename+"_FT.gif",'gif')
		img = Image(Point(int(width/2),int(height/2)), filename+"_FT.gif")
	elif return_type == "Im":
		img.save(filename+"_FT_Im.gif",'gif')
		img = Image(Point(int(width/2),int(height/2)), filename+"_FT_Im.gif")
	elif return_type == "Mag":
		img.save(filename+"_FT_Mag.gif",'gif')
		img = Image(Point(int(width/2),int(height/2)), filename+"_FT_Mag.gif")

	return img

def singlePointIFT(rec_data_mag,rec_data_phase,width,height,latticeSpacing,x,y):
	# The image should have a width and height each corresponding to a range of (-4*pi/a, 4*pi/a)
	n = [0,0]
	for i in range(height):
		for j in range(width):
			k_y = 4*pi*(2*i/height-1)/latticeSpacing
			k_x = 4*pi*(2*j/width-1)/latticeSpacing
			rec_data_phase[i][j] += k_x*x+k_y*y
			n[0] += rec_data_mag[i][j]*cos(rec_data_phase[i][j])
			n[1] += rec_data_mag[i][j]*sin(rec_data_phase[i][j])
	if y%10 == 0 and x %10 == 0:
		print((x,y))
	return n

def findInvFT(filename,rec_data,latticeSpacing,size_x,size_y):
	# print("Calculating Inverse Fourier Transform...")
	# data_re = [[0 for x in range(size_x)] for y in range(size_y)]
	# data_im = [[0 for x in range(size_x)] for y in range(size_y)]
	# pDone = 0

	# # Store rec_data as magnitudes and phases, which are easier to deal with here
	# width = len(rec_data[0][0])
	# height = len(rec_data[0])
	# rec_data_mag = [[0 for x in range(width)] for y in range(height)]
	# for i in range(height):
	# 	for j in range(width):
	# 		rec_data_mag[i][j] = sqrt(rec_data[0][i][j]**2+rec_data[1][i][j]**2)
	# rec_data_phase = [[0 for x in range(width)] for y in range(height)]
	# for i in range(height):
	# 	for j in range(width):
	# 		rec_data_phase[i][j] = atan2(rec_data[1][i][j],rec_data[0][i][j])

	rec_data_re = np.array(rec_data[0])
	rec_data_im = np.array(rec_data[1])
	rec_data = rec_data_re + 1j*rec_data_im

	data = np.fft.ifft2(rec_data)
	# data = np.fft.fftshift(data)
	data_re = data.real.tolist()
	data_im = data.imag.tolist()

	# singlePointIFT_vect = np.vectorize(singlePointIFT)
	# singlePointIFT_vect.excluded.add(0)
	# singlePointIFT_vect.excluded.add(1)
	# X = [[x for x in range(size_x)] for y in range(size_y)]
	# Y = [[y for x in range(size_x)] for y in range(size_y)]
	# data_re,data_im = singlePointIFT_vect(rec_data_mag,rec_data_phase,width,height,latticeSpacing,X,Y)
	# data_re = data_re.tolist()
	# data_im = data_im.tolist()

	# for y in range(size_y):
	# 	for x in range(size_x):
	# 		n = [0,0]
	# 		for i in range(height):
	# 			for j in range(width):
	# 				k_y = 4*pi*(2*i/height-1)/latticeSpacing
	# 				k_x = 4*pi*(2*j/width-1)/latticeSpacing
	# 				rec_data_phase[i][j] += k_x*x+k_y*y
	# 				n[0] += rec_data_mag[i][j]*cos(rec_data_phase[i][j])
	# 				n[1] += rec_data_mag[i][j]*sin(rec_data_phase[i][j])
	# 		data_re[y][x] = n[0]
	# 		data_im[y][x] = n[1]
	# 	# Print percent done while running
	# 	if round(0.1*int(1000*float(y+1)/size_y),2) > pDone:
	# 		pDone = round(0.1*int(1000*float(y+1)/size_y),2)
	# 		print(str(pDone)+"%")

	with open(filename+"_ReconstructedLattice_Re.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in data_re:
			wr.writerow(row)

	with open(filename+"_ReconstructedLattice_Im.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in data_im:
			wr.writerow(row)

	return data_re,data_im

def mag(re,im):
	return sqrt(re**2+im**2)

def createReconstructedImage(filename,braggFilteredData,size_x,size_y,return_type):
	if return_type == "Im":
		braggFilteredData = braggFilteredData[1]
	elif return_type == "Mag":
		braggFilteredData_mag = [[0 for x in range(size_x)] for y in range(size_y)]
		mag_vect = np.vectorize(mag)
		braggFilteredData_mag = mag_vect(braggFilteredData[0],braggFilteredData[1])
		# for i in range(size_y):
		# 	for j in range(size_x):
		# 		braggFilteredData_mag[i][j] = sqrt(braggFilteredData[0][i][j]**2+braggFilteredData[1][i][j]**2)
		braggFilteredData = braggFilteredData_mag.tolist()
	else: # return_type of "Re" is taken default
		return_type = "Re"
		braggFilteredData = braggFilteredData[0]

	# Find extrema
	min_point = min(min(braggFilteredData))
	max_point = max(max(braggFilteredData))

	# Draw data
	img = NewImage.new("RGB", (size_x, size_y))
	putpixel = img.putpixel
	for y in range(size_y):
		for x in range(size_x):
			color = int((braggFilteredData[y][x]-min_point)/(max_point-min_point)*255)
			putpixel((x,y),(color,color,color))

	if return_type == "Re":
		img.save(filename+"_BraggFiltered_Re.gif",'gif')
		img = Image(Point(int(size_x/2),int(size_y/2)), filename+"_BraggFiltered_Re.gif")
	if return_type == "Im":
		img.save(filename+"_BraggFiltered_Im.gif",'gif')
		img = Image(Point(int(size_x/2),int(size_y/2)), filename+"_BraggFiltered_Im.gif")
	elif return_type == "Mag":
		img.save(filename+"_BraggFiltered_Mag.gif",'gif')
		img = Image(Point(int(size_x/2),int(size_y/2)), filename+"_BraggFiltered_Mag.gif")

	return img