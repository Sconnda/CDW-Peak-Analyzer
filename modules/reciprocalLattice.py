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

def findFT(filename,peaks,size_x,size_y, width, height):
	# print("Calculating Fourier Transform...")
	data = [[0 for x in range(size_x)] for y in range(size_y)]
	for peak in peaks:
		x,y = peak
		x = min(int(round(x)),size_x-1)
		y = min(int(round(y)),size_y-1)
		data[y][x] = 1

	data = np.loadtxt(filename+".txt")
	data_list = data.tolist()
	data = np.array(data)
	min_data = min(min(data_list))
	data -= min_data

	rec_data = np.fft.fft2(data)
	rec_data = np.fft.fftshift(rec_data)
	# rec_data = np.flip(rec_data,0)
	rec_data_re = rec_data.real.tolist()
	rec_data_im = rec_data.imag.tolist()

	with open(filename+"_ReciprocalLattice.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in rec_data_re:
			wr.writerow(row)

	with open(filename+"_ReciprocalLattice_Im.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in rec_data_im:
			wr.writerow(row)

	return rec_data_re, rec_data_im

def stepFTMask(rec_data, kx_center, ky_center, radius):
	rec_data_re, rec_data_im = rec_data
	width = len(rec_data_re[0])
	height = len(rec_data_re)
	rec_data_filtered_re = [[0 for x in range(width)] for y in range(height)]
	rec_data_filtered_im = [[0 for x in range(width)] for y in range(height)]
	for i in range(height):
		ky = (i-height/2.0)/height
		for j in range(width):
			kx = (j-width/2.0)/width
			if (kx-kx_center)**2+(ky-ky_center)**2 < radius**2:
				rec_data_filtered_re[i][j] = rec_data_re[i][j]
				rec_data_filtered_im[i][j] = rec_data_im[i][j]
	return rec_data_filtered_re,rec_data_filtered_im

def gaussianFTMask(rec_data, kx_center, ky_center, radius):
	rec_data_re, rec_data_im = rec_data
	width = len(rec_data_re[0])
	height = len(rec_data_re)
	rec_data_filtered_re = [[0 for x in range(width)] for y in range(height)]
	rec_data_filtered_im = [[0 for x in range(width)] for y in range(height)]
	for i in range(height):
		ky = (i-height/2.0)/height
		for j in range(width):
			kx = (j-width/2.0)/width
			mask = exp(-4*pi*((kx-kx_center)**2+(ky-ky_center)**2)/(kx_center**2+ky_center**2))
			rec_data_filtered_re[i][j] = mask*rec_data_re[i][j]
			rec_data_filtered_im[i][j] = mask*rec_data_im[i][j]
	return rec_data_filtered_re,rec_data_filtered_im

def highPassFilterMask(kxSc,kySc):
	return (1-exp(-pi*(kxSc**2+kySc**2)/0.0001))

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

def createRecDataImage(filename,rec_data,num_peaks,width,height,return_type,extension):
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

	# Pass through a high-pass filter to increase the brightness of the edges
	rec_data = globalFTFilter(rec_data)

	# Find extrema
	min_point = min(min(line) for line in rec_data)
	max_point = max(max(line) for line in rec_data)

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
		img.save(filename+"_"+extension+"_Re.gif",'gif')
		img = Image(Point(int(width/2),int(height/2)), filename+"_"+extension+"_Re.gif")
	elif return_type == "Im":
		img.save(filename+"_"+extension+"_Im.gif",'gif')
		img = Image(Point(int(width/2),int(height/2)), filename+"_"+extension+"_Im.gif")
	elif return_type == "Mag":
		img.save(filename+"_"+extension+"_Mag.gif",'gif')
		img = Image(Point(int(width/2),int(height/2)), filename+"_"+extension+"_Mag.gif")

	return img

def mag(re,im):
	return sqrt(re**2+im**2)

def phase(re,im):
	return atan2(im,re)

def findInvFT(filename,rec_data,extension,size_x,size_y):
	rec_data_re = np.array(rec_data[0])
	rec_data_im = np.array(rec_data[1])
	rec_data = rec_data_re + 1j*rec_data_im

	rec_data = np.fft.fftshift(rec_data)

	data = np.fft.ifft2(rec_data)
	data_re = data.real.tolist()
	data_im = data.imag.tolist()

	data_mag = [[0 for x in range(size_x)] for y in range(size_y)]
	mag_vect = np.vectorize(mag)
	data_mag = mag_vect(data.real,data.imag).tolist()

	data_phase = [[0 for x in range(size_x)] for y in range(size_y)]
	phase_vect = np.vectorize(phase)
	data_phase = phase_vect(data.real,data.imag).tolist()

	with open(filename+"_"+extension+"_Re.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in data_re:
			wr.writerow(row)

	with open(filename+"_"+extension+"_Im.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in data_im:
			wr.writerow(row)

	with open(filename+"_"+extension+"_Mag.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in data_mag:
			wr.writerow(row)

	with open(filename+"_"+extension+"_RawPhase.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in data_phase:
			wr.writerow(row)

	return data_re,data_im

def createReconstructedImage(filename,braggFilteredData,size_x,size_y,return_type,extension):
	braggFilteredData_Re = braggFilteredData[0]

	braggFilteredData_Im = braggFilteredData[1]

	braggFilteredData_Mag = [[0 for x in range(size_x)] for y in range(size_y)]
	mag_vect = np.vectorize(mag)
	braggFilteredData_Mag = mag_vect(braggFilteredData[0],braggFilteredData[1]).tolist()

	braggFilteredData_Phase = [[0 for x in range(size_x)] for y in range(size_y)]
	phase_vect = np.vectorize(phase)
	braggFilteredData_Phase = phase_vect(braggFilteredData[0],braggFilteredData[1]).tolist()

	# REAL IMAGE__________________
	# Find extrema
	min_point_Re = min(min(braggFilteredData_Re))
	max_point_Re = max(max(braggFilteredData_Re))
	# Draw data
	img_Re = NewImage.new("RGB", (size_x, size_y))
	putpixel = img_Re.putpixel
	for y in range(size_y):
		for x in range(size_x):
			color = int((braggFilteredData_Re[y][x]-min_point_Re)/(max_point_Re-min_point_Re)*255)
			putpixel((x,y),(color,color,color))

	# IMAGINARY IMAGE__________________
	# Find extrema
	min_point_Im = min(min(braggFilteredData_Im))
	max_point_Im = max(max(braggFilteredData_Im))
	# Draw data
	img_Im = NewImage.new("RGB", (size_x, size_y))
	putpixel = img_Im.putpixel
	for y in range(size_y):
		for x in range(size_x):
			color = int((braggFilteredData_Im[y][x]-min_point_Im)/(max_point_Im-min_point_Im)*255)
			putpixel((x,y),(color,color,color))

	# MAGNITUDE IMAGE__________________
	# Find extrema
	min_point_Mag = min(min(braggFilteredData_Mag))
	max_point_Mag = max(max(braggFilteredData_Mag))
	# Draw data
	img_Mag = NewImage.new("RGB", (size_x, size_y))
	putpixel = img_Mag.putpixel
	for y in range(size_y):
		for x in range(size_x):
			color = int((braggFilteredData_Mag[y][x]-min_point_Mag)/(max_point_Mag-min_point_Mag)*255)
			putpixel((x,y),(color,color,color))

	# PHASE IMAGE__________________
	# Find extrema
	min_point_Phase = -pi
	max_point_Phase = pi
	# Draw data
	img_Phase = NewImage.new("RGB", (size_x, size_y))
	putpixel = img_Phase.putpixel
	for y in range(size_y):
		for x in range(size_x):
			color = int((braggFilteredData_Phase[y][x]-min_point_Phase)/(max_point_Phase-min_point_Phase)*255)
			putpixel((x,y),(color,color,color))

	img_Re.save(filename+"_"+extension+"_Re.gif",'gif')
	img_Im.save(filename+"_"+extension+"_Im.gif",'gif')
	img_Mag.save(filename+"_"+extension+"_Mag.gif",'gif')
	img_Phase.save(filename+"_"+extension+"_RawPhase.gif",'gif')

	if return_type == "Re":
		img = Image(Point(int(size_x/2),int(size_y/2)), filename+"_"+extension+"_Re.gif")
	elif return_type == "Im":
		img = Image(Point(int(size_x/2),int(size_y/2)), filename+"_"+extension+"_Im.gif")
	elif return_type == "Phase":
		img = Image(Point(int(size_x/2),int(size_y/2)), filename+"_"+extension+"_RawPhase.gif")
	else: # return_type == "Mag" by default
		img = Image(Point(int(size_x/2),int(size_y/2)), filename+"_"+extension+"_Mag.gif")

	return img
