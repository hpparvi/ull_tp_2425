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

N_ic = 10
create_ic = False
update_pars = True
make = False
execute = True
ploteo = True
fnum = 1



dt, dt_out, t_end = 0.01, 10, 10
epsilon = 0.1

# input_file = "data/stars.txt" 
input_file = "data/random_bodies.txt"
create_bodies = create_ic
N_bodies = N_ic
output_file = "output/output.txt"

lines_custom = ["# Simulation parameters", 
                "dt = %.2f       # Time step"%dt,
                "dt_out = %.2f    # Printing time step"%dt_out,
                "t_end = %.2f     # Final execution time"%t_end, 
                "epsilon = %.2f   # Softening scale"%epsilon,
                "",
                "# File names",
                "input_file = \"%s\""%input_file,
                "create_bodies = .%s. "%create_bodies,
                "N_bodies = %i "%N_bodies,
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
    
    fig, ax = plt.subplots(num = fnum, figsize = (6,6))
    for i in range(N_bodies):
        ax.plot(data[:,1+3*i],  data[:,2+3*i], color = plt.colormaps['tab10'](i), marker = '+', ls = '', label = 'm%i'%(i+1), alpha = 0.5)
    
    ax.minorticks_on()
    ax.grid()

ax.set_title('N-body sim with $\\epsilon$=%.1e\nTotal time: %.2f\nTotal number of iterations: %i'%(epsilon,t_end, t_end/dt))
ax.set_xlabel('x')
ax.set_ylabel('y')

ax.set_xlim(-1,1)
ax.set_ylim(-1,1)

ax.legend()
plt.tight_layout()
