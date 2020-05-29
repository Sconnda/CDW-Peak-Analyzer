import math
import csv

def findJMatrix(filename,peaks,size_x,size_y):
	print("Finding j-Matrix...")
	num_peaks = len(peaks)
	jMatrix = [[0 for i in range(size_x)] for j in range(size_y)]
	pDone = 0
	for y in range(size_y):
		for x in range(size_x):
			dmin = math.hypot(size_x-1, size_y-1)
			j = -1
			for i in range(num_peaks):
				xPeak,yPeak = peaks[i]
				d = math.hypot(xPeak-x, yPeak-y)
				if d < dmin:
					dmin = d
					j = i
			jMatrix[y][x] = j
		if int(100*(y+1)/size_y) > pDone:
			pDone = int(100*(y+1)/size_y)
			print(str(pDone)+"%")

	with open(filename+"_jMatrix.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in jMatrix:
			wr.writerow(row)

	return jMatrix

def triangulation(jMatrix,num_peaks,size_x,size_y):
	print("Finding triangulation...")
	bondMatrix = [[0 for i in range(num_peaks)] for j in range(num_peaks)]
	pDone = 0
	for y in range(size_y):
		for x in range(size_x):
			j = int(jMatrix[y][x])
			for dy in [-1,0,1]:
				for dx in [-1,0,1]:
					if x+dx < 0 or x+dx >= size_x or y+dy < 0 or y+dy >= size_y:
						continue
					i = int(jMatrix[y+dy][x+dx])
					if i != j:
						bondMatrix[j][i] = 1
						bondMatrix[i][j] = 1
		if 10*int(10*(y+1)/size_y) > pDone:
			pDone = 10*int(10*(y+1)/size_y)
			print(str(pDone)+"%")

	with open("VortexLattice_Triangulation.csv", 'w',newline='') as f:
		wr = csv.writer(f)
		for row in bondMatrix:
			wr.writerow(row)

	return bondMatrix