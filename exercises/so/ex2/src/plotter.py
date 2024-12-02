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

def ploteo_func(output_file, N_bodies, t_end, dt, epsilon, tree, theta,
                trail = True, movie = False,
                d3 = False, fnum=1):
    data  = np.loadtxt(output_file)
    n_it = data.shape[0]
    plt.close(fnum)
    
    if d3 == False:
        fig, ax = plt.subplots(2,2, num = fnum, figsize = (15,15))
    if d3:
        fig, ax = plt.subplots(num = fnum, figsize = (10,12), subplot_kw= {'projection': '3d'})
    
    if movie:
        n_frames = data.shape[0]
    else:
        n_frames = 1
        
    
    
    # if trail == True:
    #     trail_n = None
    # if trail == False:
    #     trail_n = -1
    # else:
    #     trail_n = trail
    
    if trail == True:
        init_frame = None
    if trail == False:
        init_frame = -1
    
    for i_f in range(n_frames):
        if movie:
            end_frame = i_f+1
            if i_f == n_frames-1:
                end_frame = None
            if type(trail) != bool:
                if i_f >= trail:
                    init_frame =  end_frame - trail
                else:
                    init_frame = None
        else:
            end_frame = None
            init_frame = n_it-trail
            i_f = -1
        
            
        if d3:
            ax.clear()
        else:
            for i, j in zip([0,0,1,1],[0,1,0,1]):
                ax[i,j].clear()
              
        for i in range(N_bodies):
            if N_bodies <=10:
                color = plt.colormaps['tab10'](i)
            else:
                color = 'black'
            
            if d3 == False:
                ax[0,0].plot(data[i_f,1+3*i],  data[i_f,2+3*i], color = color, marker = '.', ls = '', label = 'm%i'%(i+1), alpha = 0.7)
                ax[0,0].plot(data[i_f,2+3*i],  data[i_f,3+3*i], color = color, marker = '.', ls = '', label = 'm%i'%(i+1), alpha = 0.7)
                ax[0,0].plot(data[i_f,1+3*i],  data[i_f,3+3*i], color = color, marker = '.', ls = '', label = 'm%i'%(i+1), alpha = 0.7)
                ax[0,0].plot(data[init_frame:end_frame,1+3*i],  data[init_frame:end_frame,2+3*i], color = color, marker = '', lw = 0.7, alpha = 0.5)
                ax[0,1].plot(data[init_frame:end_frame,2+3*i],  data[init_frame:end_frame,3+3*i], color = color, marker = '', lw = 0.7, alpha = 0.5)
                ax[1,0].plot(data[init_frame:end_frame,1+3*i],  data[init_frame:end_frame,3+3*i], color = color, marker = '', lw = 0.7, alpha = 0.5)
            if d3:
                ax.plot(data[i_f,1+3*i],  data[i_f,2+3*i], data[i_f,3+3*i], color = color, marker = '.', ls = '', label = 'm%i'%(i+1), alpha = 0.7)
                ax.plot(data[init_frame:end_frame,1+3*i],  data[init_frame:end_frame,2+3*i], data[init_frame:end_frame,3+3*i], color = color, marker = '', lw = 0.7, alpha = 0.5)

        if movie:
            fig.suptitle('N-body sim with $\\epsilon$=%.1e\nTime: %.2f  Iteration: %i   dt=%.2e\nBarnes-Hut: %s     $\\theta$=%.2f'%(epsilon,dt*i, i, dt, tree, theta))
    
    if movie == False:
        fig.suptitle('N-body sim with $\\epsilon$=%.1e\nTime: %.2f  Iteration: %i   dt=%.2e\nBarnes-Hut: %s     $\\theta$=%.2f'%(epsilon, n_it*dt, n_it, dt, tree, theta))
        
    for i in range(2):
        for j in range(2):
            if d3 == False:
                ax[i,j].minorticks_on()
                ax[i,j].grid()
                ax[i,j].set_xlim(-radius,radius)
                ax[i,j].set_ylim(-radius,radius)
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
        
    plt.tight_layout()

ploteo_func(output_file, N_bodies, t_end, dt, epsilon, tree, theta,
            trail = 10000, movie = False,
            d3 = d3, fnum=fnum)