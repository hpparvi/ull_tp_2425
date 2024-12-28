# N-Body Barnes-Hut with OpenMP

Fortran implementation of an N-body code using Barnes-Hut algorithm and OpenMP.   

## Introduction

This code is based on Ángel de Vicente's notes "Técnicas Avanzadas de Programación 2018-2019", pgs. 74-81. The code has been modified in some important ways:
- The simulation parameters are specified in a configuration file (`.txt`) file to be provided when the program is executed. 
- Initial conditions are generated with a python script and read from a file (`.dat`).
- Instead of using a single script, we made the program modular.
- We use custom types for particle positions and velocities intead of relying on arrays. 
- We parallelized it with OpenMP.  


## Description

### Initial Conditions 
We generated two types of initial conditions (ICs) to test our code:
- Uniform Sphere Collapse. Generated with `make_ics_uniform_sphere.py`. Particles are initialized uniformly distributed in a sphere with zero velocities.
- Broucke A2 three body orbit:  Generated with `plot_broucke_a2.py`. Same as exercise 1. 

In the Uniform Sphere Collapse case we can set as many particles as we want. There the difference in performance between parallel and not parallel code is clearly seen.  

### Simulation Specifications

The parameters of the simulation are specified with a configuration file to be passed in terminal when the program is executed. Two examples for each IC can be found in `config_files/`. The parameters to be specified are:
- `G`: Gravitational constant
- `T`: Total integration time
- `DT`: Time step
- `EPSILON`: Smoothing length
- `THETA`: Barnes-Hut 'precision' parameter
- `N_SNAPSHOTS`: Number of snapshots to save
- `ICS_FILE`: Path to initial conditions to be read
- `SAVE_FOLDER`: Where to save teh data

### Code structure

We have two working fortran codes at:
- `fortran/parallel`
- `fortran/not_parallel`
  
Names are self explanatory. Parallelization is done with OpenMP. The two codebases are basically identical except for the `omp` flags in loops during integration of trajectories, and differences in the makefile to allow parallelization.

In both cases the source code is found in `/src` and it has the following modules:
- `geometry.f90`: (Module) Same as exercise 1.
- `particle.f90`: (Module) Same as exercise 1.
- `i_o_utils.f90`: (Module) Contains functions for:
  - reading the configuration file
  - reading tyhe ICs
  - create folder where we are saving snapshots
  - routine for saving snapshots
- `octree.f90` (Module) Main module with principal functions and subroutines for the implementation of Barnes-Hut algorithm. We have followed Ángel de Vicente's implementation closely, but using the custom types in geometry and particle modules instead of arrays. This part of the code is heavily commented. 
- `ex2.f90` (Program) Main program to compile. 


### Example usage

Make sure you have an ICS file (check `exercises/rp/ex2/ics/`) and a configuration file (check `exercises/rp/ex2/config_files/`)

Place yourself at `exercises/rp/ex2/fortran/parallel/`. 

Compile the code:

`>> make clean`

`>> make`

Execute 

`>> ./bin/ex2`

You will be prompted the following: 'Provide the path of config file:'
Paste the FULL path of the config file (use any of the two at ) and press enter. Then integration should start. 

You can read and plot the orbits with python scripts:
- `plot_broucke_a2.py` makes a nice figure of a three body orbit
- `plot_uniform_collapse.ipynb` makes a video of the spherical collapse 

Figures and videos are saved in `exercises/rp/ex2/figures/`

## Results

The specification of the machnie we used to run this code is:

12 CPUs:  AMD Ryzen 5 5600H with Radeon Graphics

We know the code is working because the orbit at `exercises/rp/ex2/figures/broucke_a2.png` has the shape it should have. 

We can also assertain the validity of the code by looking at the video of the spherical collapse `exercises/rp/ex2/figures/uniform_collapse.mp4`. I suspect there is something off but I think it has something to do with the parameter `THETA`. 

Regarding the difference between parallel vs. not parallell. It obviously depends on the number of bodies. For the three body orbit there is not much difference. For the spherical collapse with `N=1000` the difference is a more than twice speed up using parallel code.