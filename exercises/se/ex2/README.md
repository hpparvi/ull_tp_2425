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

First, ALWAYS do:

> make clean

to remove previous compilations. 


For non-parallel version, use the regular makefile, by writing in
your terminal either:

> make

or: 

> make -f makefile

For parallelized version, use the corresponding makefile by writing:
> make -f makefile_parallel

Both should run just fine, as the parallel sections are written
with magic commands that will be interpreted as comments if the flag
-fopenmp is not used when compiling. This is what the different 
makefiles account for. 

To run the exercise, use:
> ./ex2

or

> make test

as usual. 

To compile the copy of exercise 1 for the elapsed time tests, do make 
clean and then:
> make -f makefile_ex1


=========== Discussion of parallelization ===========

Parallelization is not optimal for very low numbers of particles, 
as it is more expensive to open parallel blocks and assign the 
threads than it is to perform the very few calculations. 

This is also the reason why some blocks of code are not parallelized
(such as nullifying pointers or deleting the tree (see below).

Additionally, the Barnes-Hut algorithm is not necessarily best in some 
cases. For low numbers of particles, using a certain theta might cause 
errors that are not worth the time it takes to calculate the accelerations
due to every single particle (this can cause some stable 3-body problems
to deviate from stability, as some of my classmates found).

Reasons why I paralellized what I did, and not other sections:
- The tree cannot be parallelized, because different threads will reach 
the same point at the same time (not always, but there is an increased 
chance for a higher number of particles and for longer integrations) and
then the error "SHOULD NOT BE HERE" will be triggered because the head 
node is not updated in time. 
- Nullifying, deleting, other quick operations are not worth 
parallelizing given that the exchange of information and distribution of 
tasks between the threads takes longer than the actual calculations. 
- Other places are too complicated to parallelize with OpenMP

What I did parallelize was the update of accelerations and the force 
calculation, which is done particle by particle (the auxiliary function 
is not parallelized for the reasons above). 

=========== Discussion of elapsed time ===========

A Python notebook is included to plot the elapsed time for several
different cases, called elapsed_time.ipynb. 
For these tests, I removed the printing statements from
the parallelized code (commenting out) because printing takes a 
significant amount of time when running code. I also commented out the 
writing to file block. I then wrote down the results manually.




These are the results on a table:

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
| Case 3: parallel Barnes-Hut |       |       |       |       |       |                      |
| N = 3                       | 5.119 | 5.080 | 5.071 | 5.070 | 5.075 |  5.08 &plusmn; 0.02  |
| N = 20                      | 0.705 | 0.711 | 0.709 | 0.705 | 0.708 | 0.708 &plusmn; 0.002 |
| N = 50                      | 1.164 | 1.153 | 1.146 | 1.143 | 1.147 | 1.151 &plusmn; 0.007 |
| N = 80                      | 1.590 | 1.581 | 1.575 | 1.603 | 1.596 |  1.59 &plusmn; 0.01  |
| N = 150                     | 2.973 | 2.978 | 2.968 | 2.960 | 2.977 | 2.971 &plusmn; 0.007 |
| N = 400                     | 9.421 | 9.386 | 9.411 | 9.420 | 9.432 |  9.41 &plusmn; 0.02  |
