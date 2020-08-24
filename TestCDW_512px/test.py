import os,sys
import numpy as np

data = np.load("TestCDW_512px.txt")
b = np.fft.fft2(data)
c = np.fft.ifft2(b)

d = (c.real-a).tolist()

print(min(min(d)))
print(max(max(d)))