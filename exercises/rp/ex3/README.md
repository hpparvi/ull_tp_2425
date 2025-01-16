# Exercise 3: N-Body Barnes-Hut with MPI


This Fortran code implements a Barnes–Hut N-body simulation with MPI-based parallelization. Most of the code is the same or similar as in exercise 2. 

The following modules are unchanged: `geometry.f90`,`particle.f90`,`octree.f90`, `i_o_utils.f90`.

We added a `mpi_utils.f90` module to define custom MPI types. 

Also, instead of passing a configuration file via input, is it passed though `make run CONFIG='path/to/config_file.txt'` instead, and defaults to a config file with spherical collapse conditions with `N=1000`.

Initial conditions and plotting are the same as in exercise 2 (via python). 

## Description of code

1. Configuration & Initialization
   - Reads a config file specifying gravitational constants (G), total time (T), time-step (dt), path to initial conditions, etc
   - Reads initial conditions (ICS) for all particles and broadcasts them to every process.

2. Barnes–Hut Tree (Built on Every Rank)
   - *Each rank* constructs the same octree from the entire particle set.
   - This is not parallelized; only the force computation uses MPI.

3. Force Computation with MPI
   - The particles are partitioned among ranks based on a simple block distribution: each rank “owns” a subrange of the global array.
   - Each rank computes gravitational forces only for its assigned particles (calling compute_forces on that subrange).
   - The accelerations (partial results) are gathered from all ranks so that everyone sees the updated acceleration array.

4. Time Integration (Leapfrog)
   - After each partial leapfrog update, there is an MPI_Allgatherv step to synchronize the global bodies.

5. Output
   - Snapshots of the system are saved at specified intervals by rank 0.

Because the tree is built redundantly on each process, and there are frequent global communications, the code can be slower than serial N-body for moderate N. It is mainly a demonstration of how to parallelize force computations with MPI while leaving the tree construction serial. For large N, distributing the force work begins to offset communication costs.

## Compilation & Execution
- in the `rp/ex3/fortran/` subfolder, run with mpirun -np <num_procs> ./bin/ex3 <config_file>
- or run `make run NP=<num_procs> CONFIG=<config_file>`

