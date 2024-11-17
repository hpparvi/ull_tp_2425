## Input files
### initial_conditions.dat
The file initial_conditions.dat contains the data for the simulation. The format is as follows:\
line	contains\
1  dt\
2	t_end\
3	dt_out\
4	m1 px1 py1 pz1 vx1 vy1 vz1\
.\
.\
.\
n+3	mn pxn pyn pzn vxn vyn vzn\
with n the number of particles we are considering.

## Output files
### result.dat
The file result.dat contains the output data. Each line has the following format:\
t  px1  py1  pz1  px2  py2  pz2  ...  pxn  pyn  pzn
