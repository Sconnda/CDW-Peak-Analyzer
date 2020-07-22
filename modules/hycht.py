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

def findPhaseData(filename,size_x,size_y,G_x,G_y):
	data = [[0 for x in range(size_x)] for y in range(size_y)]
	for peak in peaks:
		x,y = peak
		x = min(int(round(x)),size_x-1)
		y = min(int(round(y)),size_y-1)
		data[y][x] = 1


	with open(filename+"_ReciprocalLattice.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in rec_data_re:
			wr.writerow(row)

	with open(filename+"_ReciprocalLattice_Im.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in rec_data_im:
			wr.writerow(row)

	return rec_data_re, rec_data_im

# DAN

# Save displacement image data as a size_x by size_y csv file and return an array including the data. One of the inputs should be the phase image, I'll work on that
def findDisplacementData(filename,phaseData,size_x,size_y):
	return 0

# Create, save, and return an image for any real field in a graphics window
def createFieldImage(filename,fieldData,extension,size_x,size_y):
	return 0