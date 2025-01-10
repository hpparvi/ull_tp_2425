# 1. CREATE INITIAL CONDITIONS

First of all, you have to create a .txt file with the initial conditions of the problem. The file's content must be organised in 7 columns: the first one contains the value of the masses of the particles, the next three columns contain the particles' initial coordinates, and on the other hand, the last three columns contain the particles' initial velocities. You can write this kind of .txt file manually or you can use initialconditions.f90, a program that asks you to specify the number of particles (N) and their initial velocities' absolute value. Then, this program creates N particles distributed uniformly on the surface of a sphere of radius 1, having each of them an initial velocity pointing towards the center of the sphere. All the particles' velocities have the same absolute value and all the masses are equal (mass = 1/N). 

# 2. RUN THE MAIN PROGRAM

Run the main program using Makefile. Then, write ./sim to execute the program.

The main program reads the information from initial_conditions.txt and asks you to specify the timestep value, output time value, value of the final time and the number of particles. Obviously, this number of particles must be equal to the number of particles mentioned above. After executing the program, we will get a file called output.dat that contains the next information: <br/>
<br/>
time  |  particle_id  |   x_coordinate   |  y_coordinate    |  z_coordinate
<br/>
Moreover, the execution time is printed in the console.

# 3. PLOT THE DATA

With for_plotting.ipynb you can plot the initial particle distribution using initial_conditions.txt and make a GIF of the evolution of the particles using output.dat.

# 4. EXPLANATION OF THE MAIN PROGRAM

This program, named ex3.f90, implements a parallel simulation using the Barnes-Hut model, which is used to simulate gravitational interactions between particles in a three-dimensional space. It leverages MPI (Message Passing Interface) for parallel execution across multiple processors. <br/>
<br/>
The program begins by initializing MPI to set up a distributed computing environment, allowing work to be distributed across multiple processes. In such a simulation, processes are allocated different particles, and each process computes the interactions between the particles it is assigned. Results are communicated between processes using MPI collective communication operations, such as MPI_Scatter and MPI_Gather.
# Step-by-Step Explanation: 
1. MPI Initialization: the functions MPI_INIT, MPI_COMM_RANK, and MPI_COMM_SIZE are called to initialize the parallel environment, retrieve the process rank, and get the total number of processes.
2. Reading Input Parameters: the process with rank 0 reads the user input values, such as the timestep size (dt), final simulation time (t_end), output time interval (dt_out), and the total number of particles (n).
3. Distributing Particles Among Processes: using MPI_BCAST, the value of n (total number of particles) is broadcast to all processes. Then, the number of particles assigned to each process is calculated and stored in particles_per_proc.
4. Preparing Buffers: the process with rank 0 allocates the buffers buffer_send and particles, which will contain the particle data to be sent and stored by each process. The initial particle information (mass, position, and velocity) is read from a text file and loaded into the particles array.
5. Distributing Data with MPI_Scatter: MPI_Scatter is used to distribute the particle data among different processes. Each process receives its portion of particles in the buffer_recv array.
6. Rebuilding Local Particles: aach process reconstructs the local particles using the data received in buffer_recv, assigning values for mass, position, velocity, and acceleration.
7. Building and Processing the Barnes-Hut Tree: a tree structure is used to compute gravitational interactions. This tree is initialized and particles are distributed across its cells. The tree efficiently represents gravitational interactions by dividing the space into cells containing particles. The processes of "finding cells" and "placing particles in cells" are carried out by the subroutines Find_Cell and Place_Cell.
8. Force Calculations and Updates: gravitational forces between particles are computed using the Barnes-Hut tree. The accelerations are calculated, and the particle velocities and positions are updated using the timestep dt through numerical integration.
9. Main Loop: a time-stepping loop simulates the gravitational interactions and updates the particles until the final time (t_end) is reached. At each timestep: 
The particle velocities and positions are updated. 
The tree is rebuilt, and gravitational interactions are recalculated. 
Results from each process are gathered using MPI_Gather and stored in an output file.
10. Execution Time Calculation: after the simulation ends, the total execution time is calculated using system_clock and printed by the process with rank 0.
11. Finalizing MPI: the program ends by calling MPI_FINALIZE to clean up and close the MPI environment.



