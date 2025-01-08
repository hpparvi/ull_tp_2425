## Exercise 3: N-body simulation using Barnes-Hut algorithm, parallel programming with MPI

The Barnes–Hut simulation is an approximation algorithm for performing an N-body simulation. The computational time is reduced since the interactions of all the particles with each other are not calculated. Instead, the space where the particles are placed (we can think about it as a box), known as the main cell, is divided in octant (8 subcells). This division it's repeated until all the particles are placed in a cell (subcell of the multiple divisions done). At this moment you have created the tree. After that, you have to calculate the center of mass of each division in subcells done. In this way, for a cell placed in certain octant, it's only necessary to calculate direct interaction with the nearest particles. The rest of them are considered as a only one particle, where its mass is the addition of the mass of the conglomerate, and the center of mass corresponds the position of this particle. So the number of calculation needed is reducted significantly, from $O(n²)$ to $O(n  logn)$.

In this program, Barnes-Hut algorithm have been used in a three-dimensional space with n particle. The implementation is done in the module  named **tree** that contains the subroutines required to use the algorithm. 

In order to check if the results of the N-body simulation are correct, it has been used a simple python scripts to make plots and animations (2D and 3D), named as plot_file.py, animation_file.py and animation3d_file.py; it is not necessary to change the number of particles in the script (nor other parameters), it works with whatever number of particles introduced (if they are saved as output.dat and it has the same format: t px1 py1 pz1 ... pxn pyn pzn). 

In the folder init_files, there are some examples of initial conditions files (included example test of the notes). Meanwhile, part_gen folder include a little program to generate the input.dat with a random distribution of particles (with null velocities). 


This time, main code has been parallelised using MPI, which reduces computational time for a large number of particles. Main difference is that MPI allows to distribute work across processes and supports communication between nodes. Number of cores used can be modified in makefile (in case of less than four particles, must be necessary to change this number -but for simulations of few particles, like those, parallel or sequencial doesn't make a significant variation). Also, as expected, if more cores are used, if it is faster the execution. But, sometimes, when we have a small number of particles, communication between nodes slows down the calculations, and using only one core (sequencial programming) it's more efficient.

Tests have been done under the same time conditions ($dt = 0.01, t_{end} = 0.1$), as were done for OpenMP. Few time steps are used, because time cannot be parallelised in the simulation, so it is not so important to demonstrate the difference in computation for each time step. First testing record (10 particles) shows that the parallelised program works slower than without parallelise, because in this case takes more time to call OpenMP and start to parallelise than calculate directly the interactions. Also, in the last test record (N=10000), we observe that time using MPI is almost double than with OpenMP, because OpenMP simulations were done in CCA computer (has 8 cores). Meanwhile, simulations with MPI were done in a Virtual Machine with only 4 cores.
     
| N | Sequencial time (s) |Time with OpenMp (s) | Time using MPI | 
|--|--|--|--|
| 10 | 3e-3 | 8e-3 | 1.1e-2 |
| 100 | 2.7e-2 | 9.9e-3 | 8.0e-3 |
| 1000 | 1.54 | 0.416 | 0.40 |
| 10000 | 158.47 | 39.79 | 72.06 |
