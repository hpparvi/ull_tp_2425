========= About the variables in ex3 =========

Brief manual to understand what the variables 
do in the ex3.f90.

rank: the number that IDs the node (process)

comsize: how many nodes there are. It is set by
the mpirun -n flag. 

ierr: integer that stores types of errors or lack
thereof.

MPI_vector3d [and others]: the MPI necessary
types of the user-defined variables so that MPI
knows how to deal with them.

npart_per_node: the number of particles to send to
each node, excluding the last node in case they are 
not exactly divisible. 

npart_last_node: the remainder of the division, 
so that the last node receives however many 
particles are leftover. Might be unequal but
for smallish numbers of particles, it should
not be too bad. 

counts_nodes: the array that tells each node how
many particles to tackle

displacements: how many elements to move between 
node-assigned chunks

first, last: for each rank, there determine the indices
in between to grab the particles for each node.

i, j, k: loop indices

n: particle number

dt: current time instant in the simulation

t_end: time where the simulation ends

dt_out: time at which outputs are saved

t_out: time tracker that goes back to zero when the
output is saved

time_counter: integer to keep track of loops

total_timesteps: total number of iterations for the
loops

particles: particle array

particles_node: array that contains the subset of
particles for each node to carry out the leapfrog
method of integration (is different for each node).

a: accelerations

a_node: same as particles_node but for the 
accelerations.

openstatus_input: to check whether the inputs have 
been opened correctly

openstatus_output: to check whether the outputs have
been opened correctly

readstatus: to check if the things are being read

datafile: name of input file

resultfile: name of output file

count_begin: timer start in counts

count_end: timer end in counts

count_rate: conversion between counts and seconds

head: main cell of the tree

temp_cell: temporary cell to keep track of the current cell


================ About the parallelization ================

Some comments about the way I parallelized this exercise:

It is not possible to parallelize by sending the function
Calculate_Forces a subarray of the particles, because in order
to get the forces it is necessary to know the position (via the
tree) of ALL the particles. Therefore it is necessary to modify
the original Calculate_Forces function to take two extra
arguments: a start and end for the loop indices, and then, 
passing the entire particles array, each process calculates the 
total forces on that subset of particles.

============ Discussion of elapsed time ============

Below is an updated table to show the elapsed time for 
different several cases (excluding 3 particles because
I set the number of nodes to 8 for all runs for a fair
comparison to OpenMP). Again I excluded the printing
statements and writing to file. 

|     Case 1: no strategy     |       |       |       |       |       |                      |
|:---------------------------:|:-----:|:-----:|:-----:|:-----:|:-----:|:--------------------:|
| N = 3                       | 0.049 | 0.056 | 0.057 | 0.053 | 0.057 | 0.054 &plusmn; 0.003 |
| N = 20                      | 0.147 | 0.153 | 0.153 | 0.153 | 0.153 | 0.152 &plusmn; 0.002 |
| N = 50                      | 0.762 | 0.766 | 0.767 | 0.768 | 0.762 | 0.765 &plusmn; 0.003 |
| N = 80                      | 1.899 | 1.893 | 1.908 | 1.903 | 1.913 | 1.903 &plusmn; 0.007 |
| N = 150                     | 6.587 | 6.584 | 6.587 | 6.581 | 6.568 | 6.581 &plusmn; 0.007 |
| N = 400                     | 46.37 | 46.39 | 46.34 | 46.38 | 46.51 |  46.40 &plusmn; 0.06 |
|      Case 2: Barnes-Hut     |       |       |       |       |       |                      |
| N = 3                       | 0.172 | 0.172 | 0.171 | 0.170 | 0.167 | 0.170 &plusmn; 0.002 |
| N = 20                      | 0.263 | 0.260 | 0.261 | 0.261 | 0.260 | 0.261 &plusmn; 0.001 |
| N = 50                      | 0.970 | 0.972 | 0.966 | 0.979 | 0.969 | 0.971 &plusmn; 0.004 |
| N = 80                      | 1.825 | 1.831 | 1.837 | 1.838 | 1.838 | 1.834 &plusmn; 0.005 |
| N = 150                     | 4.790 | 4.838 | 4.788 | 4.819 | 4.825 |  4.81 &plusmn; 0.02  |
| N = 400                     | 23.32 | 23.45 | 23.40 | 23.37 | 23.41 |  23.39 &plusmn; 0.04 |
| Case 3: Barnes-Hut + OpenMP |       |       |       |       |       |                      |
| N = 3                       | 5.119 | 5.080 | 5.071 | 5.070 | 5.075 |  5.08 &plusmn; 0.02  |
| N = 20                      | 0.705 | 0.711 | 0.709 | 0.705 | 0.708 | 0.708 &plusmn; 0.002 |
| N = 50                      | 1.164 | 1.153 | 1.146 | 1.143 | 1.147 | 1.151 &plusmn; 0.007 |
| N = 80                      | 1.590 | 1.581 | 1.575 | 1.603 | 1.596 |  1.59 &plusmn; 0.01  |
| N = 150                     | 2.973 | 2.978 | 2.968 | 2.960 | 2.977 | 2.971 &plusmn; 0.007 |
| N = 400                     | 9.421 | 9.386 | 9.411 | 9.420 | 9.432 |  9.41 &plusmn; 0.02  |
|   Case 4: Barnes-Hut + MPI  |       |       |       |       |       |                      |
| N = 3                       |   -   |   -   |   -   |   -   |   -   |          -           |
| N = 20                      | 0.637 | 0.609 | 0.542 | 0.815 | 0.574 |  0.63 &plusmn; 0.09  |
| N = 50                      | 1.405 | 1.424 | 1.420 | 1.436 | 1.427 |  1.42 &plusmn; 0.01  |
| N = 80                      | 2.709 | 2.329 | 2.304 | 2.332 | 2.309 |   2.4 &plusmn; 0.2   |
| N = 150                     | 4.916 | 4.911 | 5.166 | 4.917 | 4.885 |   5.0 &plusmn; 0.1   |
| N = 400                     | 18.45 | 18.73 | 17.79 | 18.57 | 17.88 |  18.3 &plusmn; 0.4   |
