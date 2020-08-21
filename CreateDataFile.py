import os
import sys

import numpy as np

import csv

from PIL import Image as NewImage

# Declaring global variables_____________________
filename = "CDW_Data"
# Window and size of window
win = 0
size_x = 512
size_y = 512
scale = 1

def main():
	global filename, scale, braggFilteredData_return_type
	if len(sys.argv) > 1:
		filename = sys.argv[1]
	if len(sys.argv) > 2:
		scale = float(sys.argv[2])
	if not os.path.isdir(filename):
		filename = "TestCDW_512px"
	os.chdir(filename)

	img = NewImage.open(filename+".gif")

	# Identify width and height of data array
	global size_x, size_y
	size_x,size_y = img.size
	print("________________________________")
	print("Width: " + str(size_x))
	print("Height: " + str(size_y))
	print("________________________________")

	imgL = img.convert("L")

	# From luminosity data
	data = np.array(imgL).tolist()

	# # From red to blue color scale
	# dataL = np.array(imgL).tolist()
	# dataL = np.array(dataL)-255
	# imgHSV = img.convert("HSV")
	# dataH = np.array(imgHSV)[:,:,0].tolist()
	# dataH = np.array(dataH)-128
	# print(dataL)
	# print(np.amax(dataL))
	# print(np.amin(dataL))
	# dataH = dataH/abs(dataH)
	# data = dataH*dataL

	row1 = ["1","2"]
	with open(filename+".txt", 'w',newline='') as f:
		for row in data:
			rowP = [str(x) for x in row]
			f.write("\t".join(rowP))
			f.write("\n")

main()