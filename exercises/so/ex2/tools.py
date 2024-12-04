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

import os, glob, imageio

directory = '/Users/oscarna/Documents/Compu/ull_tp_2425/exercises/so/ex2'

os.chdir(directory)

tree = True

N_ic = 5
create_ic = False
update_pars = False
make = False
execute = False
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

N_threads = 12

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
                "output_file = \"%s\""%output_file,
                "",
                "# OpenMP settings",
                "N_threads = %i   # Number of threads to be used in parallelization"%N_threads]

if update_pars:
    with open('custom.par', 'w') as file:
        for l in lines_custom:
            file.write(l+ '\n')

if execute:
    if make:
        os.system("make clean")
        os.system("make")
    os.system("./ex2 custom.par")

#%% Data reading and plotting

def ploteo_func(output_file, N_bodies, t_end, dt, epsilon, tree, theta, radius,
                trail = True,
                movie = False, dt_movie = False,
                make_frames = True,
                make_mp4 = True, sim_name = 'animation', duration = 5,
                snapshot = False,
                d3 = False, fnum=1,
                saving_dir = '/Users/oscarna/Documents/Compu/ull_tp_2425/exercises/so/ex2/frames'):
    
    data  = np.loadtxt(output_file)
    n_it = data.shape[0]
    n_bodies_total = (data.shape[1]-1)/3
    
    if 'frames' not in os.listdir():
        os.system('mkdir frames')
    
    if N_bodies == 'all':
        N_bodies = n_bodies_total
        
    if make_frames:
        plt.close(fnum)
        
        if d3 == False:
            fig, ax = plt.subplots(2,2, num = fnum, figsize = (12,12))
    
        if d3:
            fig, ax = plt.subplots(num = fnum, figsize = (10,12), subplot_kw= {'projection': '3d'})
       
        if d3 == False:
            for i in range(2):
                for j in range(2):
                        ax[i,j].minorticks_on()
                        ax[i,j].grid()
                        ax[i,j].set_axisbelow(True)
                        ax[i,j].set_xlim(-radius,radius)
                        ax[i,j].set_ylim(-radius,radius)
                        
            ax[0,0].set_xlabel('x')
            ax[0,1].set_xlabel('y')
            ax[1,0].set_xlabel('x')
            ax[0,0].set_ylabel('y')
            ax[0,1].set_ylabel('z')
            ax[1,0].set_ylabel('z')
            
            ax[1,1].set_axis_off()
    
        if d3:
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            ax.set_xlim(-radius, radius)
            ax.set_ylim(-radius, radius)
            ax.set_zlim(-radius, radius)
            ax.minorticks_on()
            ax.grid()    
       

        if movie:
            plt.ioff()
            if dt_movie != False:
                n_frames = int(t_end/dt_movie)
                n_res = round(dt_movie / dt)
            else:
                n_frames = data.shape[0]
                n_res = 1
        else:
            n_frames = 1
            n_res = 1
            
        if trail == True:
            init_frame = None
        elif trail == False:
            init_frame = -1
        
        for i_f in range(n_frames):
            if movie:
                if dt_movie != False:
                    if i_f != range(n_frames)[-1]:
                        i_f = int(i_f * n_res)
                    else:
                        i_f = n_it-1
                    if type(trail) != bool:
                        if i_f ==0:
                            trail = trail*n_res
                end_frame = i_f+1
                if i_f/n_res == n_frames-1:
                    end_frame = None
                if type(trail) != bool:
                    if i_f >= trail:
                        if end_frame!=None:
                            init_frame =  end_frame - trail
                    else:
                        init_frame = None
            
            else:
                if snapshot != False:
                    i_f = snapshot
                    n_it = i_f
                    end_frame = i_f+1
                else:
                    end_frame = None
                
                if type(trail) != bool:
                    init_frame = n_it-trail
                i_f = -1
            
            if snapshot != False:
                i_f = snapshot
            
                
            if d3:
                ax.lines[:] = []  #.clear()
            else:
                for i, j in zip([0,0,1,1],[0,1,0,1]):
                    for line in ax[i,j].lines:
                        line.remove()
                    for collection in ax[i,j].collections:
                        collection.remove()
                    # ax[i,j].lines[:] = []  #.clear()
            
            if N_bodies <=10:
                color = plt.colormaps['tab10'](np.arange(N_bodies))
            else:
                color = 'black'        
            
            if d3 == False:
                ax[0,0].scatter(data[i_f,1::3],  data[i_f,2::3], color = color, marker = '.', ls = '', alpha = 0.7)
                ax[0,1].scatter(data[i_f,2::3],  data[i_f,3::3], color = color, marker = '.', ls = '', alpha = 0.7)
                ax[1,0].scatter(data[i_f,1::3],  data[i_f,3::3], color = color, marker = '.', ls = '', alpha = 0.7)
                if trail != False:
                    if N_bodies<=10:
                        for i in range(N_bodies):
                            ax[0,0].plot(data[init_frame:end_frame:n_res,1+3*i],  data[init_frame:end_frame:n_res,2+3*i],color = color[i], marker = '', lw = 0.7, alpha = 0.5)
                            ax[0,1].plot(data[init_frame:end_frame:n_res,2+3*i],  data[init_frame:end_frame:n_res,3+3*i],color = color[i], marker = '', lw = 0.7, alpha = 0.5)
                            ax[1,0].plot(data[init_frame:end_frame:n_res,1+3*i],  data[init_frame:end_frame:n_res,3+3*i],color = color[i], marker = '', lw = 0.7, alpha = 0.5)
                    else:
                        ax[0,0].plot(data[init_frame:end_frame:n_res,1::3],  data[init_frame:end_frame:n_res,2::3], color = color, marker = '', lw = 0.7, alpha = 0.5)
                        ax[0,1].plot(data[init_frame:end_frame:n_res,2::3],  data[init_frame:end_frame:n_res,3::3], color = color, marker = '', lw = 0.7, alpha = 0.5)
                        ax[1,0].plot(data[init_frame:end_frame:n_res,1::3],  data[init_frame:end_frame:n_res,3::3], color = color, marker = '', lw = 0.7, alpha = 0.5)
          
            if d3:
                ax.scatter(data[i_f,1::3],  data[i_f,2::3], data[i_f,3::3], color = color, marker = '.', ls = '', alpha = 0.7)
                if trail != False:
                    ax.plot(data[init_frame:end_frame:n_res,1::3],  data[init_frame:end_frame:n_res,2::3], data[init_frame:end_frame:n_res,3::3], color = color, marker = '', lw = 0.7, alpha = 0.5)
                    
            if movie:
                fig.suptitle('N-body sim with $N=$%i;   $\\epsilon$=%.1e\nTime: %.2f;   Iteration: %i;   dt=%.2e\nBarnes-Hut: %s;   $\\theta$=%.2f'%(n_bodies_total, epsilon,dt*i_f, i_f, dt, tree, theta))
                plt.savefig(os.path.join(saving_dir, 'frames/frame_'+f'{int(i_f/n_res):04}'))
                if i_f/n_res%50 == 0:
                    print('Frame '+f'{int(i_f/n_res):04} saved.')
                
        if movie == False:
            fig.suptitle('N-body sim with $N=$%i;   $\\epsilon$=%.1e\nTime: %.2f;   Iteration: %i;   dt=%.2e\nBarnes-Hut: %s;   $\\theta$=%.2f'%(n_bodies_total, epsilon, t_end, n_it, dt, tree, theta))
        
        if movie:
            plt.ion()
    if movie:        
        if make_mp4:
            video_name = os.path.join(saving_dir, sim_name+'.mp4')
            frames = sorted(glob.glob(os.path.join(saving_dir,'frames/frame_')+'*'))
            gif_duration = duration/len(frames)
            figures = []
            for f in frames:
                fram = imageio.imread(f)
                figures.append(fram)
            output_params = {'fps': 1/gif_duration, 'codec': 'libx264'}
            imageio.mimwrite(video_name, figures, **output_params)
        
    plt.tight_layout()

# N_bodies = 1000
# ploteo_func(output_file, N_bodies, t_end, dt, epsilon, tree, theta,
#             snapshot = 100,
#             trail = 10, movie = False,
#             d3 = d3, fnum=fnum)

N_bodies = 3
output_file = "output/output.txt"
dt, dt_out, t_end = 0.001, 10, 500
epsilon = 0.0
theta = 0
radius = 1.1



ploteo_func(output_file, N_bodies, t_end, dt, epsilon, tree, theta, radius,
            # snapshot = 1001, 
            dt_movie = 0.03,
            trail = True, 
            movie = False,
            saving_dir = '/Users/oscarna/Documents/Compu/ull_tp_2425/exercises/so/ex2/',
            make_frames = True,
            make_mp4 = True,
            sim_name = '3body',
            duration = 10,
            d3 = False, fnum=1)
