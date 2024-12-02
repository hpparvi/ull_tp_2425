# Barnes-Hut algorithm
This code implements the Barnes-Hut algorithm, an approximation method for $N$-body simulations that optimizes the calculation of forces between particles leading to reduce the computational complexity from $\textit{O}(N^2)$ to $\textit{O}(N \ logN)$. The algorithm achieves this by recursively grouping the $N$ bodies into regions and organizing them in an octree (3D) or a quadtree (2D). Each node in the tree represents a specific region of space: the root node corresponds to the entire space, while its child nodes represent subdivisions into eight octants. This process continues recursively, dividing each region into octants until every subdivision contains one or no particles. To calculate the forces between particles, the tree is traversed hierarchically considering that particles in distant cells are treated as a single large particle centered at the center of mass of the cell, while particles in nearby cells are treated individually, and the process continues recursively. The decision to treat a cell as a single particles depends on $\theta$, a threshold parameter. Small er values of $\theta$ treat more cells individually and provide more precise results, but the computational cost is incremented. 

# Code features
To perform an $N$-body simulation, the code follows these steps: 

1. Build the tree
2. Calculate the forces
3. Update positions and velocities of the particles using the calculated forces 

where the first two steps are approched by the Barnes-Hut algorithm. 

The code is structured into 5 modules and a main program, ```ex2.f90```. These modules are: 

* ```geometry.f90``` - stores basic mathematical objects and calculations
* ```particle.f90``` - stores properties about the particles such as position, velocity, acceleration and mass
* ```barnes.f90``` - stores the Barnes-Hut algorithm (tree building and forces calculation)
* ```calcs.f90``` - stores calculations of particle position, velocity and acceleration
* ```data.f90```  - allows to read and save data to files

The initial conditions for the simulations are stored in their corresponding input files in the directory named ```inputs```. To properly compile the code you need to move them to its directory. To compile and run the code you can use 
~~~
make -f Makefile_ex2
~~~
or 
~~~
make -f Makefile_ex2_omp
~~~
depending if you want to use the code parallelized with OpenMP or not. To remove compilations use
~~~
make -f Makefile_ex2 clean
~~~
or 
~~~
make -f Makefile_ex2_omp clean
~~~

In case you don't want to use the default input files, you can generate one with random body data by using the ```random_input.f90``` program located in the ```generate input``` directory. To compile and run it write
~~~
make -f Makefile_random_input
~~~
To remove compilations use 
~~~
make -f Makefile_random_input clean
~~~

# Parallelization with OpenMP