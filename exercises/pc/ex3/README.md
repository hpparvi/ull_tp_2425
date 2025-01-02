# Barnes-Hut algorithm
This code implements the Barnes-Hut algorithm, an approximation method for $N$-body simulations that optimizes the calculation of forces between particles leading to reduce the computational complexity from $\textit{O}(N^2)$ to $\textit{O}(N \ logN)$. The algorithm achieves this by recursively grouping the $N$ bodies into regions and organizing them in an octree (3D) or a quadtree (2D). Each node in the tree represents a specific region of space: the root node corresponds to the entire space, while its child nodes represent subdivisions into 8 octants. This process continues recursively, dividing each region into octants until every subdivision contains one or no particles. To calculate the forces between particles, the tree is traversed hierarchically considering that particles in distant cells are treated as a single large particle centered at the center of mass of the cell, while particles in nearby cells are treated individually, and the process continues recursively. The decision to treat a cell as a single particles depends on $\theta$, a threshold parameter. Smaller values of $\theta$ treat more cells individually and provide more precise results, but the computational cost is incremented. 

# Code features
To perform an $N$-body simulation, the code follows these steps: 

1. Build the tree
2. Calculate the forces
3. Update positions and velocities of the particles using the calculated forces 
4. 

where the first two steps are approched by the Barnes-Hut algorithm. 

The code is structured into 6 modules and a main program, ```ex3.f90```. These modules are: 

* ```geometry.f90``` - stores basic mathematical objects and calculations
* ```particle.f90``` - stores properties about the particles such as position, velocity, acceleration and mass
* ```barnes.f90```   - stores the Barnes-Hut algorithm (tree building and forces calculation)
* ```calcs.f90```    - stores calculations of particle position, velocity and acceleration
* ```data.f90```     - allows to read and save data to files
* ```mpi_types.f90```- allows to create MPI dataypes for 3D vectors, 3D points and more complex structures like a 3D particle

The initial conditions for the simulations are stored in their corresponding input files in the directory named ```inputs```. To properly compile the code you need to move them to the directory where the code is. To compile and run the code you can use 
~~~
make -f Makefile_ex3
~~~
~~~
mpirun -n m ex3
~~~
To remove compilations use
~~~
make -f Makefile_ex3 clean
~~~

In case you don't want to use the default input files, you can generate one with random body data by using the ```random_input.f90``` program located in the ```generate input``` directory. To compile and run it write
~~~
make -f Makefile_input
~~~
To remove compilations use 
~~~
make -f Makefile_input clean
~~~

# Parallelization with MPI


Several time tests were made on a 4-cores computer using a time step of $dt = 0.01$, an output time step of $dt_{out} = 0.1$ and a final time of $t_{final} = 50$. The results are plotted in the figure ```times.png``` located in the ```plots``` directory, and are shown in the following table:

| Particles | Direct-sum time (ex1) [s]| Non-parallelized Barnes-Hut time [s]| Parallelized Barnes-Hut time with OpenMP [s]| Parallelized Barnes-Hut time with MPI [s]|
|:---------:|:-----:               |:--------:                        |:--------------:|:--------------:|
|     5     |   2   |     3    |        8       |                                   |           
|     10    |   2   |     4    |        7       |                                   |    
|     25    |   3   |     5    |        9       |                                   |    
|     50    |   3   |     6    |        9       |                                   |    
|    100    |   5   |     7    |       24       |                                   |    
|    250    |   27  |    29    |       34       |                                   |    
|    500    |  108  |    71    |       51       |                                   |    
|    1000   |  441  |    172   |       90       |                                   |    

For a small number of particles, the direct-sum algorithm is computationally faster as the Barnes-Hut has to build the tree and recursively computing the center of mass and total mass for each cell. However, for a larger number of particles the Barnes-Hut algorithm becomes more efficient due its optimized method for calculating forces between particles. As expected, the OpenMP-parallelized version of the Barnes-Hut algorithm takes more time to perform simulations when dealing with a small number of particles. This is because parallelization introduces overhead that outweights its benefits. On the contrary, for simulations involving a large number of particles, the parallelized code becomes computationally faster as the benefits of parallel processing -such as distributing workloads across multiple cores- outweigh the initial overhead. 