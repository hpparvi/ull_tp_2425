This code uses the leapfrog integration method to simulate a N-body problem and the Barnes-Hut algorithm to optimize the calculations.

The main program, main.f90, utilizes the following modules:

	- geometry.f90: This module defines all the mathematical operations for 3D points and vectors.
	- definitions.f90: This module defines various types, including particles, pointers, ranges, and cells.
	- barnes_hut.f90: This module contains all the tree-related subroutines.
	- datatypes_mpi.f90: This module defines MPI datatypes for point3d, vector3d, and particle structures, enabling their use in parallel communication.
	- env_mpi.f90: This module provides essential subroutines for initializing, finalizing, and managing the MPI environment. 

The main program is based on the Calculate_forces subroutine, which computes the gravitational forces acting on each particle. It uses an input file called data_input.dat, which provides the initial data (time step, output data time step, total time, number of bodies, and the mass, position, and acceleration of each particle). The program also uses an output file to store the positions of each particle, printed and saved according to the requested format. Finally, the program repeats the calculations based on the initial parameters defined in the particulas.f90 program. For example:

	- dt = 0.01: The system advances in increments of 0.01 seconds for each iteration of the time loop. The smaller the value, the greater the precision.
	- dt_out = 1.0: The system prints results at this time interval, meaning the program outputs results every 1 second of simulation time.
	- t_end = 10: Total simulation time. In this case, the simulation runs for 10 seconds.

The particulas.f90 program generates the initial conditions for the simulation. It asks the user for the number of particles (n) and generates random values for mass (between 0 and 100), positions (between -100 and 100), and velocities (between -1 and 1) for each particle. These values are stored in the data_input.dat file. The ranges for mass, position, and velocity can be modified manually. Since the random_number function generates a number between 0 and 1, you can simply adjust the factors to change these ranges.

Please note that the time it takes to run the parallelized program depends not only on the parallelization but also on the number of cores and threads used. In my case, my computer has 2 cores, so the program may run faster on computers with higher computational capacity. Additionally, the performance depends on the value of theta, which is a proportionality factor between the distance of a particle to the center of mass of a cell and the size of the cell. Essentially, a smaller theta results in the tree being subdivided more, which increases the precision but also makes the computation slower and more resource-intensive. Conversely, a larger theta will result in fewer subdivisions, making the computation faster but less precise. This value of theta is relative to the positions of the particles; since my particles were within the range of -100 to 100, I chose a theta value of 50.
The only difference with the second exercise is that this code is parallelized with MPI, leading to the following changes in the main program:

	- Initialization: MPI is initialized, and the necessary datatypes are generated.
	- Master Thread: The master process starts timing the execution and reads the input data. It also allocates memory for the particles.
	- Data Broadcasting: The data and particle information are broadcast to all processes using MPI_Bcast. This ensures that all processes have the same initial data.
	- Particle Distribution: The particles_nodes subroutine distributes particles among processes, determining how many particles each process will handle based on the total number of processes.
	- Data Gathering: The gather_particles subroutine collects particle data from all processes and consolidates it into a global particle array.
	- Main Loop: The main loop uses the same subroutines as explained earlier in the sequential version of the code.

In this folder, you can find an image (MPI_timers.png) showing the time taken by the code to perform calculations for a set of particles. The results clearly demonstrate that the MPI parallelization is the fastest method among all four approaches.


INSTRUCTIONS:

To use the code properly, first create the initial conditions by running:

	$ gfortran particulas.f90 -o particulas
	$ ./particulas

Then, select the number of particles you want to simulate. After that, compile the main program and execute it by running:

	$ make
	$ ./main

The program will then run. Please note that these instructions are intended for Linux and may not apply to other operating systems. If you want to clear the .o, use the following instruction:

	$ make clean
