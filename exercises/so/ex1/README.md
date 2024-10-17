# Basic Dynamical Simulation, by Oscar Soler Pérez
Here is included the fortran code written to perform a direct numerical integration of a system consisting of gravitationally interacting particles using the leapfrog integration method.
The code has been inspired by the notes of the course Técnicas Avanzadas de Programación 2018-2019, page 34.

It contains:
- a geometry module, which defines point and vector data type and operations among them
- a particles module, which defines particle data type with mass, position, velocity and acceleration
- a leapfrog program with the numerical integration method
