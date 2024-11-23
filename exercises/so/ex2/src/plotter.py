#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 22:34:56 2024

@author: oscarna
"""

"""
-Initial condition visualization
-...
"""

import numpy as np
import matplotlib.pyplot as plt

import os

directory = '/Users/oscarna/Documents/Compu/ull_tp_2425/exercises/so/ex2'

os.chdir(directory)

tree = True

N_ic = 5
create_ic = False
update_pars = True
make = False
execute = True
ploteo = True
d3 = True
fnum = 1


dt, dt_out, t_end = 0.01, 10, 100
epsilon = 0.1
theta = 1

# input_file = "data/stars.txt" 
input_file = "data/random_bodies.txt"
create_bodies = create_ic
N_bodies = N_ic
radius = 10
output_file = "output/output.txt"

lines_custom = ["# Simulation parameters", 
                "dt = %.2f       # Time step"%dt,
                "dt_out = %.2f    # Printing time step"%dt_out,
                "t_end = %.2f     # Final execution time"%t_end, 
                "epsilon = %.5f   # Softening scale"%epsilon,
                "theta = %.1f    # l/D ratio compared to theta"%theta,
                "",
                "# File names",
                "input_file = \"%s\""%input_file,
                "create_bodies = .%s. "%create_bodies,
                "N_bodies = %i "%N_bodies,
                "radius = %.f"%radius,
                "output_file = \"%s\""%output_file]

if update_pars:
    with open('custom.par', 'w') as file:
        for l in lines_custom:
            file.write(l+ '\n')

if execute:
    if make:
        os.system("make clean")
        os.system("make")
    os.system("./ex2 custom.par")

# Data reading and plotting



if ploteo:
    data  = np.loadtxt(output_file)
    plt.close(fnum)
    
    if d3 == False:
        fig, ax = plt.subplots(2,2, num = fnum, figsize = (15,15))
    if d3:
        fig, ax = plt.subplots(num = fnum, figsize = (10,12), subplot_kw= {'projection': '3d'})
        
    for i in range(N_bodies):
        if N_bodies <=10:
            color = plt.colormaps['tab10'](i)
        else:
            color = 'black'
        
        if d3 == False:
            ax[0,0].plot(data[:,1+3*i],  data[:,2+3*i], color = color, marker = '.', ls = '', label = 'm%i'%(i+1), alpha = 0.5)
            ax[0,1].plot(data[:,2+3*i],  data[:,3+3*i], color = color, marker = '.', ls = '', label = 'm%i'%(i+1), alpha = 0.5)
            ax[1,0].plot(data[:,1+3*i],  data[:,3+3*i], color = color, marker = '.', ls = '', label = 'm%i'%(i+1), alpha = 0.5)
        if d3:
            ax.plot(data[:,1+3*i],  data[:,2+3*i], data[:,3+3*i], color = color, marker = '.', ls = '', label = 'm%i'%(i+1), alpha = 0.5)

    for i in range(2):
        for j in range(2):
            if d3 == False:
                ax[i,j].minorticks_on()
                ax[i,j].grid()
                ax[i,j].set_xlim(-radius,radius)
                ax[i,j].set_ylim(-radius,radius)

    fig.suptitle('N-body sim with $\\epsilon$=%.1e\nTotal time: %.2f\nTotal number of iterations: %i\nBarnes-Hut: %s     $\\theta$=%.2f'%(epsilon,t_end, t_end/dt, tree, theta))
    if d3 == False:
        ax[0,0].set_xlabel('x')
        ax[0,1].set_xlabel('y')
        ax[1,0].set_xlabel('x')
        ax[0,0].set_ylabel('y')
        ax[0,1].set_ylabel('z')
        ax[1,0].set_ylabel('z')
    if d3:
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_xlim(-radius, radius)
        ax.set_ylim(-radius, radius)
        ax.set_zlim(-radius, radius)
        ax.minorticks_on()
        ax.grid()
    
    if N_bodies <=10:
        if d3 == False:
            ax[1,0].legend()
        
    plt.tight_layout()
