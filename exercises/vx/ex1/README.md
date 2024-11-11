This code uses the leapfrog integration method to simulate a N-body problem. The main program ex1.f90 uses the modules geometry.f90 and particle.f90.
	- geometry.f90 is the module where all the mathematical operations are defined for 3d points and vectors.
	- particle.f90 is the module where the particles are defined as a type, which includes its mass, position and velocity.

The main program is based on this loop:
  DO i = 1, n
     DO j = i + 1, n
        rji = distance(p(i)%p, p(j)%p)
        r2 = rji%x**2 + rji%y**2 + rji%z**2
        r3 = r2 * sqrt(r2)
        a(i) = a(i) + ((p(j)%m * rji) / r3)
        a(j) = a(j) - ((p(i)%m * rji) / r3)
     END DO
  END DO
which computes the acceleration that every particle suffers from the gravity force because of the other particles. Also, it uses an input file called "data_input.dat" /
where the initial data (time step, output data time step, total time, number of bodies, and the mass, position, and acceleration of each particle) is given, as well as/
an output file where the positions of each particle are printed and saved following the pattern (1st particle, 2nd particle, ..., n particle, 1st particle, 2nd particle,/
...). Lastly, the program repeats the do while loop (calculates the simulation) based on the initial times you defined, for example:

 - dt = 0.01. The system will advance in increments of 0.01 seconds for each iteration of the time loop. The smaller the value, the greater the precision.
 - dt_out = 1.0. The system will print results in this time interval. In this case, the program will output results every 1 second of simulation time.
 - t_end = 10. Total time of the simulation. In this case, the simulation will run for 10 seconds.





