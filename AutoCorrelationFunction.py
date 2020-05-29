import os
import sys

import numpy as np
import matplotlib.pyplot as plt
# import math
# import random

import csv

# from graphics import *
# from PIL import Image as NewImage

filename = "TestCDW_512px"

def store(data):
	# Plot data
	X = []
	Y = []
	Y_max = data[0][1]
	for line in data:
		X.append(line[0])
		Y.append(line[1]/Y_max)
	plt.plot(X,Y,'bo')
	# plt.ylim((0, 1))
	plt.title("_Translational Order Correlation Function")
	plt.xlabel("r")
	plt.ylabel("G(r)")
	plt.show()
	# Store G(r) as a csv file
	with open(filename+"_Translational_130.csv", 'w',newline='') as f:
		wr = csv.writer(f, quoting=csv.QUOTE_ALL)
		for i in range(len(X)):
			x = X[i]
			y = Y[i]
			wr.writerow((x,y))

def main():
	global filename
	if len(sys.argv) > 1:
		filename = sys.argv[1]
	if not os.path.isdir(filename):
		filename = "TestCDW_512px"
	os.chdir(filename)

	data = np.loadtxt(filename+"_Translational_130.txt")
	store(data)

main()