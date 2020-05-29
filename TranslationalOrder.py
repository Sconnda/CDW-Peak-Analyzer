import os
import sys
import keyboard

import numpy as np
import matplotlib.pyplot as plt
import math
import random
import csv

from graphics import *
from PIL import Image as NewImage

from modules.bondFunctions import findJMatrix
from modules.bondFunctions import triangulation
from modules.reciprocalLattice import findFT
from modules.reciprocalLattice import createRecDataImage

# Declaring global variables_____________________
filename = "CDW_Data"
# Cutoff radius for orientational order correlation function G6(r)
cutoff = 100
# Window and size of window
win = 0
size_x = 512
size_y = 512
scale = 1

pi = np.pi
sqrt = np.sqrt
sin = np.sin
cos = np.cos
ceil = np.ceil

findingVectorsG = True
inverseVectorsG = False

def retJMatrix(peaks,size_x,size_y):
	noJMatrix = not os.path.exists(filename+"_jMatrix.csv")

	if noJMatrix:
		return findJMatrix(filename,peaks,size_x,size_y)

	jMatrix = []
	with open(filename+"_jMatrix.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			jMatrix.append(row)

	return jMatrix

def retTriangulation(jMatrix,n_peaks,size_x,size_y):
	noTriangulation = not os.path.exists(filename+"_Triangulation.csv")

	if noTriangulation:
		return triangulation(jMatrix,n_peaks,size_x,size_y)

	bondMatrix = []
	with open(filename+"_Triangulation.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			bondMatrix.append(row)

	return bondMatrix

def retFT(peaks,latticeSpacing):
	noFT = not os.path.exists(filename+"_ReciprocalLattice.csv")

	if noFT:
		return findFT(filename,peaks,latticeSpacing,505,505)

	rec_data = []
	with open(filename+"_ReciprocalLattice.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			rec_data.append(row)

	return rec_data

def background():
	noImage = not os.path.exists(filename+".gif")

	if noImage:
		img = NewImage.new('RGB',(int(size_x*scale),int(size_y*scale)),(255,255,255))
		img.save(filename+"_Scaled.gif","gif")
		img = Image(Point(int(scale*size_x/2),int(scale*size_y/2)), filename+"_Scaled.gif")
		return img

	img = NewImage.open(filename+".gif")
	img = img.resize((int(scale*size_x),int(scale*size_y)),NewImage.ANTIALIAS)
	img.save(filename+"_Scaled.gif","gif")
	img = Image(Point(int(scale*size_x/2),int(scale*size_y/2)), filename+"_Scaled.gif")
	return img

def FTImage(rec_data,num_peaks):
	noFTImage = not os.path.exists(filename+"_FT.gif")

	if noFTImage:
		createRecDataImage(filename,rec_data,num_peaks)

	img = NewImage.open(filename+"_FT.gif")
	img = img.resize((505,505),NewImage.ANTIALIAS)
	img.save(filename+"_FTScaled.gif","gif")
	img = Image(Point(253,253), filename+"_FTScaled.gif")
	return img

def storeReciprocalLatticeVectors():
	global findingVectorsG
	findingVectorsG = False

def findReciprocalLatticeVectors(peaks,num_peaks,latticeSpacing):
	global inverseVectorsG
	rec_data = retFT(peaks,latticeSpacing)
	width = len(rec_data[0])
	height = len(rec_data)
	imgFT = FTImage(rec_data,num_peaks)
	winFT = GraphWin('CDW Data - Fourier Transform (Select Lattice Vectors)', 505, 505)
	imgFT.draw(winFT)

	keyboard.add_hotkey('s',lambda: storeReciprocalLatticeVectors())
	vectorsG = []
	center = Point(253,253)
	searchSizeX = int(ceil(30*width/505.0))
	searchSizeY = int(ceil(30*height/505.0))
	mark = Circle(center,3)
	mark.setFill(color_rgb(255,255,255))
	mark.setOutline(color_rgb(255,255,255))
	mark.draw(winFT)
	while findingVectorsG:
		searchCenter = winFT.getMouse()

		r0 = (int(searchCenter.getX()*width/505.0),int(searchCenter.getY()*height/505.0))
		maxVal = rec_data[r0[1]][r0[0]]
		maxPoint = r0
		for dy in range(searchSizeY):
			for dx in range(searchSizeX):
				if r0[1]+dy+1 < height and r0[0]+dx+1 < width:
					val = rec_data[r0[1]+dy+1][r0[0]+dx+1]
					if val > maxVal:
						maxVal = val
						maxPoint = (r0[0]+dx+1,r0[1]+dy+1)
				if r0[1]-dy-1 >= 0 and r0[0]-dx-1 >= 0:
					val = rec_data[r0[1]-dy-1][r0[0]-dx-1]
					if val > maxVal:
						maxVal = val
						maxPoint = (r0[0]-dx-1,r0[1]-dy-1)

		pt = Point(int(maxPoint[0]*505.0/width),int(maxPoint[1]*505.0/height))

		if findingVectorsG == False:
			winFT.close()
			break

		mark = Circle(pt,3)
		ln = Line(center,pt)
		if not inverseVectorsG:
			vectorsG.append([(pi*(pt.getX()-center.getX())/(63*latticeSpacing),pi*(pt.getY()-center.getY())/(63*latticeSpacing))])
			mark.setFill(color_rgb(255,0,0))
			mark.setOutline(color_rgb(255,0,0))
			ln.setOutline(color_rgb(255,0,0))
		else:
			vectorsG[-1].append((pi*(pt.getX()-center.getX())/(63*latticeSpacing),pi*(pt.getY()-center.getY())/(63*latticeSpacing)))
			mark.setFill(color_rgb(0,255,0))
			mark.setOutline(color_rgb(0,255,0))
			ln.setOutline(color_rgb(0,255,0))
		mark.draw(winFT)
		ln.draw(winFT)
		inverseVectorsG = not inverseVectorsG

	return vectorsG

def psi(origin, peak,vecG):
	x0,y0 = origin
	x,y = peak
	re = 0
	im = 0
	for G in vecG:
		Gx,Gy = G
		re += cos(Gx*(x-x0)+Gy*(y-y0))
		im += sin(Gx*(x-x0)+Gy*(y-y0))
	re /= 2
	im /= 2
	return (re,im)

def translationalOrder(peaks,edges,vectorsG,latticeSpacing):
	N_pairs = len(vectorsG)*(len(vectorsG)-1)*2

	G_r_re = [[0 for i in range(cutoff)] for j in range(N_pairs)]
	G_r_im = [[0 for i in range(cutoff)] for j in range(N_pairs)]
	N_r = [[0 for i in range(cutoff)] for j in range(N_pairs)]

	print("Calculating translational order...")
	percentDone = 0
	pairIndex = 0
	for vSet1 in range(len(vectorsG)):
		for vSet2 in range(len(vectorsG)):
			if vSet2 <= vSet1:
				continue
			for v1 in vectorsG[vSet1]:
				for v2 in vectorsG[vSet2]:
					vecG = [v1,v2]
					for i,peak in enumerate(peaks):
						pDone = 5*int(20*(float(pairIndex)+float(i+1)/len(peaks))/N_pairs)
						if pDone > percentDone:
							percentDone = pDone
							print(str(pDone)+"% Done")

						isEdge = i in edges
						if isEdge:
							continue

						re,im = (1,0)

						(x,y) = peak
						for j,peak2 in enumerate(peaks):
							if peak == peak2:
								continue

							isEdge = j in edges
							if isEdge:
								continue

							(x2,y2) = peak2
							distance = sqrt((x2-x)**2+(y2-y)**2)/latticeSpacing
							if distance >= cutoff or distance < 0.5:
								continue

							re2,im2 = psi(peak,peak2,vecG)
							g_r = (re*re2+im*im2, re2*im-re*im2)

							rBin = int(distance-0.5)
							G_r_re[pairIndex][rBin] += g_r[0]
							G_r_im[pairIndex][rBin] += g_r[1]
							N_r[pairIndex][rBin] += 1

					for rBin in range(cutoff):
						if N_r[pairIndex][rBin] > 0:
							G_r_re[pairIndex][rBin] = G_r_re[pairIndex][rBin]/N_r[pairIndex][rBin]
							G_r_im[pairIndex][rBin] = G_r_im[pairIndex][rBin]/N_r[pairIndex][rBin]

					pairIndex += 1

	return (G_r_re,G_r_im)

def storeG(G_r):
	maxLen = len(G_r[0])
	for index in range(len(G_r)):
		G_r0 = 2*G_r[index][0]-G_r[index][1]
		for i in range(maxLen):
			G_r[index][i] /= G_r0
		# Save G6(r) plot as an image
		plt.plot([i+1 for i in range(maxLen)],G_r[index],'bo')
		plt.ylim((0, 1))
		plt.title("Translational Order Correlation Function")
		plt.xlabel("r")
		plt.ylabel("G(r)")
		plt.savefig(filename+"_TranslationalOrder.png")

	img = NewImage.open(filename+"_TranslationalOrder.png")
	width,height = img.size

	GWin = GraphWin('G(r)', width, height)
	GPlot = Image(Point(width/2,height/2), filename+"_TranslationalOrder.png")
	GPlot.draw(GWin)
	GWin.getMouse()
	GWin.close()

	# Store G6(r) as a csv file
	with open(filename+"_TranslationalOrder.csv", 'w',newline='') as f:
		wr = csv.writer(f, quoting=csv.QUOTE_ALL)
		for i in range(maxLen):
			row = [i+1]
			for index in range(len(G_r)):
				row.append(G_r[index][i])
			wr.writerow(row)

	img = NewImage.open(filename+"_TranslationalOrder.png")
	width,height = img.size

def main():
	global filename, scale, drawVoronoi
	if len(sys.argv) > 1:
		filename = sys.argv[1]
	if len(sys.argv) > 2:
		scale = float(sys.argv[2])
	if not os.path.isdir(filename):
		filename = "TestCDW_512px"
	os.chdir(filename)

	data = np.loadtxt(filename+".txt")
	global size_x, size_y
	size_x = len(data[0])
	size_y = len(data)

	peaks = []
	with open(filename+"_Peaks.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			x,y = row
			x = float(x)
			y = float(y)
			peaks.append([x,y])
	num_peaks = len(peaks)

	global win
	jMatrix = retJMatrix(peaks,size_x,size_y)
	win = GraphWin('CDW Data', int(size_x*scale), int(size_y*scale))
	img = background()
	img.draw(win)

	for peak in peaks:
		x,y = peak
		pt = Circle(Point(scale*x,scale*y),2)
		pt.setFill(color_rgb(255,0,0))
		pt.setOutline(color_rgb(255,0,0))
		pt.draw(win)

	defects = []
	bondMatrix = retTriangulation(jMatrix,len(peaks),size_x,size_y)
	for j,peak in enumerate(peaks):
		x,y = peak
		num_bonds = sum(bondMatrix[j])
		if num_bonds != 6:
			defects.append(peak)

	latticeSpacing = 0
	total_bonds = 0
	edges = []
	for i,peak in enumerate(peaks):
		x,y = peak
		for j,bonded in enumerate(bondMatrix[i]):
			if bonded == 1:
				x2,y2 = peaks[j]
				total_bonds += 1
				latticeSpacing += sqrt((x2-x)**2 + (y2-y)**2)
		if y < 0.05*size_y or y > 0.95*size_y:
			edges.append(i)
		elif x < 0.05*size_x or x > 0.95*size_x:
			edges.append(i)

	latticeSpacing /= total_bonds

	for defect in defects:
		x,y = defect
		pt = Circle(Point(scale*x,scale*y),5)
		pt.setFill(color_rgb(255,255,0))
		pt.setOutline(color_rgb(255,255,0))
		pt.draw(win)		

	for j in range(len(edges)):
		i = edges[j]
		x,y = peaks[i]
		pt = Circle(Point(scale*x,scale*y),5)
		pt.setFill(color_rgb(0,255,0))
		pt.setOutline(color_rgb(0,255,0))
		pt.draw(win)

	vectorsG = findReciprocalLatticeVectors(peaks,num_peaks,latticeSpacing)
	print("Lattice constant from FFT reciprocal lattice vectors:")
	for vPair in vectorsG:
		v = vPair[0]
		print(4*pi/(sqrt(3)*(v[0]**2+v[1]**2)**0.5))
	print("Lattice constant from FFT reciprocal lattice vectors (inverse):")
	for vPair in vectorsG:
		v = vPair[1]
		print(4*pi/(sqrt(3)*(v[0]**2+v[1]**2)**0.5))
	print("Lattice constant from bond length average:")
	print(latticeSpacing)
	G_r = translationalOrder(peaks,edges,vectorsG,latticeSpacing)

	win.getMouse()
	win.close()

	G_r_Clipped = G_r[0]
	for i in range(len(G_r_Clipped)):
		G_r_Clipped[i] = G_r_Clipped[i][0:20]

	storeG(G_r_Clipped)

main()