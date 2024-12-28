========= About the variables in ex3 =========

Brief manual to understand what the variables 
do in the ex3.f90.

rank: the number that IDs the node (process)

comsize: how many nodes there are. It is set by
the mpirun -n flag. 

ierr: integer that stores types of errors or lack
thereof.

offsets_vector: this stores the offsets in memory
of the vector3d type. It tells MPI how many bytes
to grab for each of the vector3d components.

example_vector: just a vector to get the sizes
of the components in bytes for the offsets.

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

vel_msg_per: the velocities to send to the first nodes

vel_msg_last: the velocities to send to the last node

a_msg_per: the accelerations to send to the first nodes

a_msg_last: the accelerations to send to the last node

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

a: accelerations

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

1. It is not possible to parallelize by sending the function
Calculate_Forces a subarray of the particles, because in order
to get the forces it is necessary to know the position (via the
tree) of ALL the particles. Therefore it is necessary to modify
the original Calculate_Forces function to take two extra
arguments: a start and end for the loop indices, and then, 
passing the entire particles array, each process calculates the 
total forces on that subset of particles.
