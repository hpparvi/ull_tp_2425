# 1. CREATE INITIAL CONDITIONS

First of all, you have to create a .txt file with the initial conditions of the problem. The file's content must be organised in 7 columns: the first one contains the value of the masses of the particles, the next three columns contain the particles' initial coordinates, and on the other hand, the last three columns contain the particles' initial velocities. You can write this kind of .txt file manually or you can use initialconditions.f90, a program that asks you to specify the number of particles (N) and their initial velocities' absolute value. Then, this program creates N particles distributed uniformly on the surface of a sphere of radius 1, having each of them an initial velocity pointing towards the center of the sphere. All the particles' velocities have the same absolute value and all the masses are equal (mass = 1/N). 

# 2. RUN THE MAIN PROGRAM

Run the main program using Makefile or doing this:  <br/>
<br/>
gfortran -c geometry.f90 particle.f90 barneshut.f90 <br/>
gfortran ex2.f90 geometry.o particle.o barneshut.o -o ex2 <br/>
./ex2 <br/>
<br/>

The main program reads the information from initial_conditions.txt and asks you to specify the timestep value, output time value, value of the final time and the number of particles. Obviously, this number of particles must be equal to the number of particles mentioned above. After executing the program, we will get a file called output.dat that contains the next information: <br/>
<br/>
time  |  particle_id  |   x_coordinate   |  y_coordinate    |  z_coordinate
<br/>
Moreover, the execution time is printed in the console.

# 3. PLOT THE DATA

With for_plotting.ipynb you can plot the initial particle distribution using initial_conditions.txt and make a GIF of the evolution of the particles using output.dat.
