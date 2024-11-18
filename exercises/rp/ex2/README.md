# N-Body Barnes-Hut with MPI

## Description

This code is based on Ángel de Vicente's notes "Técnicas Avanzadas de Programación 2018-2019", pgs. 74-81. The code has been modified in some important ways:
- The simulation parameters are specified in a `.txt` file to be provided when the program is executed. 
- Initial conditinos are read from a file.
- Instead of using a single script, we made the program modular.
- We use custom types for particle positions and velocities intead of relying on arrays. 
- We parallelized it with OpenMP.  


## Usage
Use a config.txt.