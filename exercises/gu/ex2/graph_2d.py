#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 15:40:19 2024

@author: urmagt2
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

path = '/home/urmagt2/ull_tp_2425/exercises/gu/ex2/'

data = np.loadtxt(path + 'output.dat')

tt = data[:,0]
numsteps = len(tt)
numparts = int(np.shape(data[:,1:])[1]/3)
parts = data[:,1:]
parts = parts.reshape((numsteps, numparts, 3))

xs = parts[:,:,0]
ys = parts[:,:,1]
zs = parts[:,:,2]

ax = plt.subplot(projection = '3d')

for ii in range(numparts):
    ax.plot(xs[:,ii],ys[:,ii],zs[:,ii])

snap = data[0][1:].reshape(numparts,3)
xx, yy, zz = snap.T