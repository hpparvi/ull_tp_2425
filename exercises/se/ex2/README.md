This folder contains the solution for the second assignment of the 
Técnicas de Programación course (2024-2025). The exercise is contained 
in the ex2.f90 file, while the barnes_hut.f90, geometry.f90 and 
particle.f90 files contain modules that are necessary for the proper
execution of ex2.f90. The exercise solves an N-body problem with the 
aid of a Barnes-Hut algorithm and OpenMP parallelization. 

The initial conditions are contained in the input files, which are
named: input_2b, input_3b and input_20b, all as .txt files, containing 
data for 2, 3 and 20 particles, respectively. The one with 3 bodies is
a stable solution for the three-body problem. Another set of random
conditions may be created using the Python notebook write_particle_file.

The output is stored in output.txt, which is used in the plot_nbody
notebook. This notebook shows the trajectories in 2D, 3D and an
animation.

=========== TO COMPILE AND RUN ===========

For non-parallel version, use the regular makefile, by writing in
your terminal either:
>> make 
or: 
>> make -f makefile

For parallelized version, use the corresponding makefile by writing:
>> make -f makefile_parallel

Both should run just fine, as the parallel sections are written
with magic commands that will be interpreted as comments if the flag
-fopenmp is not used when compiling. This is what the different 
makefiles account for. 

To run the exercise, use:
>> ./ex2 
or 
>> make test
as usual. 


=========== Discussion of parallelization ===========

