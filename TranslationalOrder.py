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
# Borders of image to select edge peaks
border_x = 0.1
border_y = 0.07
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

FTWidth = 1001
FTHeight = 1001
FTImgWidth = 505
FTImgHeight = 505

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
	global FTImgWidth, FTImgHeight
	noFT = not os.path.exists(filename+"_ReciprocalLattice.csv") or not os.path.exists(filename+"_ReciprocalLattice_Im.csv")

	if noFT:
		return findFT(filename,peaks,latticeSpacing,FTWidth,FTHeight)

	rec_data_re = []
	rec_data_im = []
	with open(filename+"_ReciprocalLattice.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			rec_data_re.append(row)

	with open(filename+"_ReciprocalLattice_Im.csv",newline='') as file:
		reader = csv.reader(file,delimiter=',',quotechar='|')
		for row in reader:
			for i,x in enumerate(row):
				row[i] = float(x)
			rec_data_im.append(row)

	return rec_data_re, rec_data_im

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

# Return an image that shows the Fourier transform in k-space
def FTImage(rec_data,num_peaks,return_type):
	global FTImgWidth, FTImgHeight
	if return_type == "Re":
		noFTImage = not os.path.exists(filename+"_FT.gif")
	elif return_type == "Im":
		noFTImage = not os.path.exists(filename+"_FT_Im.gif")
	elif return_type == "Mag":
		noFTImage = not os.path.exists(filename+"_FT_Mag.gif")

	if noFTImage:
		createRecDataImage(filename,rec_data,num_peaks,FTImgWidth,FTImgHeight,return_type)

	if return_type == "Re":
		img = NewImage.open(filename+"_FT.gif")
		img = img.resize((FTImgWidth,FTImgHeight),NewImage.ANTIALIAS)
		img.save(filename+"_FTScaled.gif","gif")
		img = Image(Point(int(FTImgWidth/2.0),int(FTImgHeight/2.0)), filename+"_FTScaled.gif")
	elif return_type == "Im":
		img = NewImage.open(filename+"_FT_Im.gif")
		img = img.resize((FTImgWidth,FTImgHeight),NewImage.ANTIALIAS)
		img.save(filename+"_FTScaled.gif","gif")
		img = Image(Point(int(FTImgWidth/2.0),int(FTImgHeight/2.0)), filename+"_FTScaled.gif")
	elif return_type == "Mag":
		img = NewImage.open(filename+"_FT_Mag.gif")
		img = img.resize((FTImgWidth,FTImgHeight),NewImage.ANTIALIAS)
		img.save(filename+"_FTScaled.gif","gif")
		img = Image(Point(int(FTImgWidth/2.0),int(FTImgHeight/2.0)), filename+"_FTScaled.gif")

	return img

def storeReciprocalLatticeVectors():
	global findingVectorsG
	findingVectorsG = False

def findReciprocalLatticeVectors(peaks,num_peaks,latticeSpacing):
	global inverseVectorsG, FTImgWidth, FTImgHeight
	# Store the Fourier transform of the peaks. The function "retFT" returns a tuple with the real and imaginary parts
	rec_data = retFT(peaks,latticeSpacing)
	rec_data_mag = [[0 for x in range(dataWidth)] for y in range(dataHeight)]
	for i in range(dataHeight):
		for j in range(dataWidth):
			rec_data_mag[i][j] = sqrt(rec_data[0][i][j]**2+rec_data[1][i][j]**2)
	rec_data = rec_data_mag

	width = len(rec_data[0][0])
	height = len(rec_data[0])
	imgFT = FTImage(rec_data,num_peaks,"Mag")
	winFT = GraphWin('CDW Data - Fourier Transform (Select Lattice Vectors)', FTImgWidth, FTImgHeight)
	imgFT.draw(winFT)

	keyboard.add_hotkey('s',lambda: storeReciprocalLatticeVectors())
	vectorsG = []
	center = Point(int(FTImgWidth/2.0),int(FTImgHeight/2.0))
	searchSizeX = int(ceil(60.0*width/FTImgWidth))
	searchSizeY = int(ceil(60.0*height/FTImgHeight))
	mark = Circle(center,3)
	mark.setFill(color_rgb(255,255,255))
	mark.setOutline(color_rgb(255,255,255))
	mark.draw(winFT)
	while findingVectorsG:
		searchCenter = winFT.getMouse()

		if findingVectorsG == False:
			winFT.close()
			break

		bounds = [Point(searchCenter.getX()-int(searchSizeX*float(FTImgWidth)/(2*width)),searchCenter.getY()-int(searchSizeY*float(FTImgHeight)/(2*height))),Point(searchCenter.getX()-int(searchSizeX*float(FTImgWidth)/(2*width)),searchCenter.getY()+int(searchSizeY*float(FTImgHeight)/(2*height))),Point(searchCenter.getX()+int(searchSizeX*float(FTImgWidth)/(2*width)),searchCenter.getY()-int(searchSizeY*float(FTImgHeight)/(2*height))),Point(searchCenter.getX()+int(searchSizeX*float(FTImgWidth)/(2*width)),searchCenter.getY()+int(searchSizeY*float(FTImgHeight)/(2*height)))]
		for bound in bounds:
			boundPoint = Circle(bound,3)
			boundPoint.setFill(color_rgb(0,0,255))
			boundPoint.draw(winFT)

		r0 = (int(searchCenter.getX()*width/float(FTImgWidth)),int(searchCenter.getY()*height/float(FTImgHeight)))

		# # Reciprocal lattice vector from brightest peak in FT
		# maxPoint = [0,0]
		# maxVal = 0
		# for dy in range(searchSizeY):
		# 	for dx in range(searchSizeX):
		# 		x = r0[0]-int(searchSizeX/2)+dx
		# 		y = r0[1]-int(searchSizeY/2)+dy
		# 		if x <= 0 or y <= 0 or x > width or y > height:
		# 			continue
		# 		val = rec_data[y][x]
		# 		if val > maxVal:
		# 			maxVal = val
		# 			maxPoint[0] = x
		# 			maxPoint[1] = y
		# endPoint = (maxPoint[0]*float(FTImgWidth)/width,maxPoint[1]*float(FTImgHeight)/height)

		# # Reciprocal lattice vector from average location of peak weighted by brightness in FT
		# cmPoint = [0,0]
		# mass = 0
		# for dy in range(searchSizeY):
		# 	for dx in range(searchSizeX):
		# 		x = maxPoint[0]-int(searchSizeX/2)+dx
		# 		y = maxPoint[1]-int(searchSizeY/2)+dy
		# 		if x <= 0 or y <= 0 or x > width or y > height:
		# 			continue
		# 		val = rec_data[y][x]
		# 		mass += val
		# 		cmPoint[0] += val*x
		# 		cmPoint[1] += val*y
		# cmPoint[0] /= mass
		# cmPoint[1] /= mass
		# endPoint = (cmPoint[0]*float(FTImgWidth)/width,cmPoint[1]*float(FTImgHeight)/height)

		# Reciprocal lattice vector from minimizing variation in psi
		optPoint = [0,0]
		minVar = 100
		print(searchSizeY)
		for dy in range(searchSizeY):
			for dx in range(searchSizeX):
				x = r0[0]-int(searchSizeX/2)+dx
				y = r0[1]-int(searchSizeY/2)+dy
				if x <= 0 or y <= 0 or x > width or y > height:
					continue
				vecG = [(8*pi*(float(x)/width-0.5)/latticeSpacing,8*pi*(float(y)/height-0.5)/latticeSpacing),(0,0)]
				var = psi_var(peaks,vecG)[0]
				if var < minVar:
					minVar = var
					optPoint[0] = x
					optPoint[1] = y
			print(dy)
		endPoint = (optPoint[0]*float(FTImgWidth)/width,optPoint[1]*float(FTImgHeight)/height)

		pt = Point(int(endPoint[0]),int(endPoint[1]))

		if findingVectorsG == False:
			winFT.close()
			break

		mark = Circle(pt,3)
		ln = Line(center,pt)
		if not inverseVectorsG:
			vectorsG.append([(8*pi*(endPoint[0]-float(FTImgWidth)/2)/(float(FTImgWidth)*latticeSpacing),8*pi*(endPoint[1]-float(FTImgHeight)/2)/(float(FTImgHeight)*latticeSpacing))])
			mark.setFill(color_rgb(255,0,0))
			mark.setOutline(color_rgb(255,0,0))
			ln.setOutline(color_rgb(255,0,0))
		else:
			vectorsG[-1].append((8*pi*(endPoint[0]-float(FTImgWidth)/2)/(float(FTImgWidth)*latticeSpacing),8*pi*(endPoint[1]-float(FTImgHeight)/2)/(float(FTImgHeight)*latticeSpacing)))
			mark.setFill(color_rgb(0,255,0))
			mark.setOutline(color_rgb(0,255,0))
			ln.setOutline(color_rgb(0,255,0))
		mark.draw(winFT)
		ln.draw(winFT)
		inverseVectorsG = not inverseVectorsG

	# Only for perfect CDW, lattice spacing of 30 pixels
	# G_mag = 2*pi/(15*sqrt(3))
	# dG_scale = 0
	# theta = pi/6
	# dG = [dG_scale*G_mag*cos(theta),dG_scale*G_mag*sin(theta)]
	# vectorsG = []
	# vectorsG.append([(dG[0],G_mag+dG[1]),(0,-G_mag+dG[1])])
	# vectorsG.append([(G_mag*sin(pi/3)+dG[0],G_mag*cos(pi/3)+dG[1]),(-G_mag*sin(pi/3)+dG[0],-G_mag*cos(pi/3)+dG[1])])
	# vectorsG.append([(G_mag*sin(pi/3)+dG[0],-G_mag*cos(pi/3)+dG[1]),(-G_mag*sin(pi/3)+dG[0],G_mag*cos(pi/3)+dG[1])])

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

def psi_var(peaks, vecG):
	psi_var = 0
	psi_peaks = [(0,0) for peak in peaks]

	psi_ave = [0,0]
	for i,peak in enumerate(peaks):
		psi_peak = psi(peaks[0], peak, vecG)
		psi_ave[0] += psi_peak[0]
		psi_ave[1] += psi_peak[1]
		psi_peaks[i] = psi_peak
	psi_ave[0] /= len(peaks)
	psi_ave[1] /= len(peaks)

	for i,peak in enumerate(psi_peaks):
		psi_var += (psi_peaks[i][0]-psi_ave[0])**2 + (psi_peaks[i][1]-psi_ave[1])**2
	psi_var /= len(psi_peaks)

	return psi_var,psi_ave


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
	print(size_x)
	print(size_y)

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
		if y < border_y*size_y or y > (1-border_y)*size_y:
			edges.append(i)
		elif x < border_x*size_x or x > (1-border_x)*size_x:
			edges.append(i)
	latticeSpacing /= total_bonds

	latticeSpacingVar = 0
	for i,peak in enumerate(peaks):
		x,y = peak
		for j,bonded in enumerate(bondMatrix[i]):
			if bonded == 1:
				x2,y2 = peaks[j]
				total_bonds += 1
				latticeSpacingVar += (sqrt((x2-x)**2 + (y2-y)**2)-latticeSpacing)**2
	latticeSpacingVar /= total_bonds
	print("________________")
	print("Lattice spacing (px): "+str(latticeSpacing))
	print("Variance of lattice spacing (px^2): "+str(latticeSpacingVar))
	print("Number of bonds used: "+str(total_bonds))
	print("________________")

	# latticeSpacing = 22.425

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

	# Store optimal vectors as a csv file
	with open(filename+"_OptimalQVectors.csv", 'w',newline='') as f:
		wr = csv.writer(f, quoting=csv.QUOTE_ALL)
		for vPair in vectorsG:
			row = [vPair[0][0],vPair[0][1],vPair[1][0],vPair[1][1]]
			wr.writerow(row)

	latticeSpacing = 22.425/2

	G_r = translationalOrder(peaks,edges,vectorsG,latticeSpacing)

	win.getMouse()
	win.close()

	G_r_Clipped = G_r[0]
	# for i in range(len(G_r_Clipped)):
	# 	G_r_Clipped[i] = G_r_Clipped[i][0:20]

	N_pairs = len(vectorsG)*(len(vectorsG)-1)*2
	psi_var_v = [0 for j in range(N_pairs)]
	pairIndex = 0
	for vSet1 in range(len(vectorsG)):
		for vSet2 in range(len(vectorsG)):
			if vSet2 <= vSet1:
				continue
			for v1 in vectorsG[vSet1]:
				for v2 in vectorsG[vSet2]:
					vecG = [v1,v2]
					psi_var_v[pairIndex] = psi_var(peaks,vecG)
					pairIndex += 1

	print(psi_var_v)

	storeG(G_r_Clipped)

main()