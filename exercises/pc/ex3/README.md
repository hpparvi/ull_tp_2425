
# Barnes-Hut algorithm
This code implements the Barnes-Hut algorithm, an approximation method for $N$-body simulations that optimizes the calculation of forces between particles leading to reduce the computational complexity from $\textit{O}(N^2)$ to $\textit{O}(N \ logN)$. The algorithm achieves this by recursively grouping the $N$ bodies into regions and organizing them in an octree (3D) or a quadtree (2D). Each node in the tree represents a specific region of space: the root node corresponds to the entire space, while its child nodes represent subdivisions into 8 octants. This process continues recursively, dividing each region into octants until every subdivision contains one or no particles. To calculate the forces between particles, the tree is traversed hierarchically considering that particles in distant cells are treated as a single large particle centered at the center of mass of the cell, while particles in nearby cells are treated individually, and the process continues recursively. The decision to treat a cell as a single particles depends on $\theta$, a threshold parameter. Smaller values of $\theta$ treat more cells individually and provide more precise results, but the computational cost is incremented. 

# Code features
The code is a distributed version of the Barnes-Hut algorithm for performing an $N$-body simulation using MPI, where each process handles a subset of the particles and performes parallel computations. The code follows these steps:

1. ```Module import```: the program begins by importing necessary modules, including the MPI-related ones.

2. ```MPI initialization```: the master process (rank = 0) initializes the MPI environment to set up communication between processes.

3. ```Creation of MPI datatype```: an MPI datatype for the particle3d type is created.

4. ```System clock initialization```: the master process initializes the system clock to measure the execution time of the simulation.

5. ```Reading input data```: the master process reads the simulation parameters and particle data from the input file. After reading the data, it distributes the simulation parameters to all processes and broadcasts particle information using MPI.

6. ```Build the initial tree and particle distribution```: all processes build the initial tree and calculate how many particles they should handle. The total number of particles is divided as equally possible among the processes. Any leftover particles are assigned to the first processes.

7. ```Initialization of particle accelerations```: each process sets the acceleration of its assigned particles to zero, and performs an all-to-all gather of particle data using MPI_Allgatherv to ensure that each process has access to the full dataset.

8. ```Force calculation and acceleration update```: each process calculates the forces between its particles using the Barnes-Hut tree algorithm. The forces are computed based on the tree of particles built in the previous step. Then, the processes gather particle data again using MPI_Allgatherv to update all particles' information. 

9. ```Opening output file```: the master process opens the output file to store the results of the simulation. 

10. ```Main simulation loop```: the main loop to update the properties of the particles begins. In each iteration: 

- Each process updates the velocity and position of its particles using the calculated forces. The particles' data is gather across all processes via MPI_Allgatherv. 
- After each update, all processes rebuild the tree to reflect the new particles' positions. 
- Each process resets the acceleration of its particles to zero and performs an all-to-all gather of particle data across all processes using MPI_Allgatherv again to ensure synchronization. 
- After resetting accelerations, each process calculates the new acceleration of its particles based on the updated tree structure. Then, each process updates the velocity of its particles. Once again, particle data is gathered across all processes to keep dataset synchronized. 
- Each process updates the simulation time. The master process checks if the current time has reached the output time. If it has, it saves the current time and the position of each particle to the ouput file. 

These steps inside the loop are repeated iteratively until the final time is reached. 

11. ```Finalization```: once the simulation reaches the final time, the master process closes the output file and calculates the total elapsed time for the simulation. Then, outputs the time taken for the simulation and an end of program message. 

12. ```Deallocation and MPI finalization```: finally, the custon MPI datatype for particle3d is deallocated and the MPI environment is finalized. 

The steps where the tree is build and the forces are calculated, are approched by the Barnes-Hut algorithm. 

The code is structured into 7 modules and a main program, ```ex3.f90```. These modules are: 

* ```geometry.f90``` - stores basic mathematical objects and calculations.
* ```particle.f90``` - stores properties about the particles such as position, velocity, acceleration and mass.
* ```barnes.f90```   - stores the Barnes-Hut algorithm (tree building and forces calculation).
* ```calcs.f90```    - stores calculations of particle position, velocity and acceleration.
* ```data.f90```     - allows to read and save data to files.
* ```mpi_types.f90```- allows to create MPI dataypes for 3D vectors, 3D points and more complex structures like a 3D particle.
* ```mpi_env.f90```  - stores subroutines that allow working in the MPI environment.

The initial conditions for the simulations are stored in their corresponding input files located in the directory named ```inputs```. To properly compile the code, you need to move these files to the directory where the code is located. Once the files are in place, you can compile the code by running
~~~
make -f Makefile_ex3
~~~
After compiling, you can execute it with the following command
~~~
mpirun -n m ex3
~~~
where *m* represents the number of processes you want to use. 

To remove compilations, run
~~~
make -f Makefile_ex3 clean
~~~

In case you don't want to use the default input files, you can generate one with random body data by using the ```random_input.f90``` program located in the ```generate input``` directory. To compile and execute it, run
~~~
make -f Makefile_input
~~~
To remove compilations, use 
~~~
make -f Makefile_input clean
~~~

# Runtime of the code

Several time tests were performed using different versions of the code on a 4-cores computer. The tests were carried out with a time step of $dt = 0.01$, an output time step of $dt_{out} = 0.1$, and a final time of $t_{final} = 50$. The results are plotted in the figure ```times.png``` located in the ```plots``` directory, and are shown in the following table:

| Particles | Direct-sum time (ex1) [s]| Non-parallelized Barnes-Hut time [s]| Parallelized Barnes-Hut time with OpenMP [s]| Parallelized Barnes-Hut time with MPI [s]|
|:---------:|:-----:               |:--------:                        |:--------------:|:--------------:|
|     5     |   2   |     3    |        8       |                  2                 |           
|     10    |   2   |     4    |        7       |                  2                 |    
|     25    |   3   |     5    |        9       |                  3                 |    
|     50    |   3   |     6    |        9       |                  4                 |    
|    100    |   5   |     7    |       24       |                  5                 |    
|    250    |   27  |    29    |       34       |                  9                 |    
|    500    |  108  |    71    |       51       |                  20                |    
|    1000   |  441  |    172   |       90       |                  51                |    

For a small number of particles, the direct-sum algorithm is computationally faster as the Barnes-Hut has to build the tree and recursively compute the center of mass and total mass for each cell. However, for a larger number of particles the Barnes-Hut algorithm becomes more efficient due its optimized method for calculating forces between particles. As expected, the OpenMP-parallelized version of the Barnes-Hut algorithm takes more time to perform simulations when dealing with a small number of particles. This is because parallelization introduces overhead that outweights its benefits. On the contrary, for simulations involving a large number of particles, the parallelized code becomes computationally faster as the benefits of parallel processing -such as distributing workloads across multiple cores- outweigh the initial overhead. 

The MPI-parallelized version of the Barnes-Hut algorithm significantly outperforms the OpenMP version in terms of execution speed. This is because OpenMP is designed to parallelize tasks within a single node, using a set of cores for these simulations. In contrast, MPI distributes the workload across multiple processes, making use of more resources and thereby reducing execution time. As a result, MPI handles large numbers of particles more effectively by distributing them across processes, preventing overload on a single node’s resources. OpenMP, on the other hand, struggles to manage larger particle sets within the limited resources of a single node.