# Basic Dynamical Simulation, by Oscar Soler Pérez
Here is included the fortran code written to perform a direct numerical integration of a system consisting of gravitationally interacting particles using the leapfrog integration method.
The code has been inspired by the notes of the course Técnicas Avanzadas de Programación 2018-2019, page 34.

It contains:
- a geometry module, which defines point and vector data type and operations among them
- a particles module, which defines particle data type with mass, position, velocity and acceleration
- a leapfrog program with the numerical integration method

Executing the code with the command:
~~~
ex1 stars.txt
~~~
generates the *output.txt* file used to create figures 1 and 2.

Executing the *plotting.py* file with python (in an environment with gfortran up and running) also creates the output and then creates the figures.

The original numerical scheme by Ángel has not been removed from the code. It is used to compare the and check the results obtained by the new required method.

The *fig_2* files available contain comparisons among both methods. Initially single precision values were used, but in *fig_.._improved_..* calculations were performed with double precision in all variables.  
It can be seen that there is a difference among both methods, but negligible (e-11). Anyhow, it increases with iterations, so one has to keep it in mind.  
The cause of these differences was not found. Even though operations seem to be written in the same order (having the same rounding errors) the compiler might be changing them to an optimal order.

