There are four distinct examples. To execute one of them, you need to change its extension to .dat.

1) initial_conditions.500
This file describes 500 particles of the same mass, placed following a uniform random distribution within a sphere of radius 1, and with 45% of the velocity they would have in a stable circular orbit if all the mass were at the center of mass. The result is a simulation in which a cloud of particles with diverse orbits forms, similar to what might be observed in globular clusters.

2) initial_conditions.inf
This file describes 3 particles that create a stable orbit, where the three particles orbit each other, forming the shape of an infinity symbol. The simulation runs until t=100t=100, maintaining stable orbits.

3) initial_conditions.3p
This file contains the positions and velocities of three particles of the same mass in a pseudo-stable orbit. Numerical errors and rounding in the method disrupt the stability, resulting in a "dance" where the gravitational interaction between the particles can be clearly observed. The method's precision can be improved by reducing the value of theta in the barnes_hut module.

4) initial_condition.2gal
This file contains the posicions and velocities of 5.000 particles,2500 placed following a uniform random distribution within a sphere of radius 1 centered in (-3,0,0) ant the other 2500 placed following a uniform random distribution within a sphere of radius 1 centered in (3,0,2) and with 45% of the velocity they would have in a stable circular orbit if all the mass were at the center of mass of each sfere. The result is a rought model of "colision" of "two galaxies"
