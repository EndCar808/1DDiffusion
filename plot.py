import matplotlib
matplotlib.use('Tkagg')
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import itertools

data = np.loadtxt("./Data/C-Data/variables.txt")
print(data.shape)
print(data)

exact = np.loadtxt("./Data/C-Data/exact.txt")
print(exact.shape)

soln = np.loadtxt("./Data/C-Data/solution.txt")
print(soln.shape)



X_MAX = int(data[0]);
dx    = float(data[1]);
dt    = float(data[2]);
alpha = float(data[3]);
T     = float(data[4]);
numtsteps = int(T / dt)

x = np.asarray([ -np.pi + i*dx for i in range(0, X_MAX)])


for i in range(0, numtsteps):	 
	print("Iter: {}/{} -- Error: {}".format(i, numtsteps, np.linalg.norm(exact[i, :] - soln[i, :])))