# Exercise 2: N-body simulation using Barnes-Hut algorithm

The Barnes–Hut simulation is an approximation algorithm for performing an N-body simulation. The computational time is reduced since the interactions of all the particles with each other are not calculated. Instead, the space where the particles are placed (we can think about it as a box), known as the main cell, is divided in octant (8 subcells). This division it's repeated until all the particles are placed in a cell (subcell of the multiple divisions done). At this moment you have created the tree. After that, you have to calculate the center of mass of each division in subcells done. In this way, for a cell placed in certain octant, it's only necessary to calculate direct interaction with the nearest particles. The rest of them are considered as a only one particle, where its mass is the addition of the mass of the conglomerate, and the center of mass corresponds the position of this particle. So the number of calculation needed is reducted significantly, from $O(n²)$ to $O(n  logn)$.

In this program, Barnes-Hut algorithm have been used in a three-dimensional space with n particle. The implementation is done in the module  named **tree** that contains the subroutines required to use the algorithm.


The code has been parallelized using OpenMP, which reduces the computational time for a large number of particles. Tests have been done under the same time conditions ($dt = 0.01, t_{end} = 0.5$). Few time steps are used since time cannot be parallelized in the simulation so it is not so important to demonstrate the difference in computation for each time step. 
     
| N | Time (without OpenMp) | Time (with OpenMP) |
|--|--|--|
| 100 |  |  |
| 1000 |  |  |
| 10000 | sc |  |
