PROGRAM tree           !This program simulates the N-body problem using the Barnes-Hut algorithm.

  use mpi_f08          !This module provides access to the MPI library with modern Fortran
  use geometry         !This module defines 3D vector and point operations for vector3d and point3d types.
  use definitions      !This module defines key data structures for N-body simulations.
  use barneshut        !This module implements the Barnes-Hut algorithm.
  use mpi_types        !This module defines MPI datatypes for point3d, vector3d, and particle structures, enabling their use in parallel communication. 
  use mpi_env          !This module provides essential subroutines for initializing, finalizing, and managing the MPI environment.
  IMPLICIT NONE

  ! Variables for time measurement.
  integer(int64) :: start, finish, rate
  real :: elapsed_time

  ! Simulation parameters and arrays.
  integer(int64) :: i
  type(vector3d), allocatable :: a(:)           !Acceleration of each particle.
  type(particle3d), allocatable :: p(:)         !Array of particles.
  TYPE (CELL), POINTER :: head, temp_cell       !Head of the octree and a temporary cell pointer.
  type(vector3d) :: rji                         !Vector representing the distance between particles.

  INTEGER :: comsize, rank, ierr !Total number of processes, processes IDs and variable for error handling in MPI operations
  INTEGER :: i_start, i_end, n_local, n_extra !First and last index of the range of particles to process, number of particles per process and extra particles
  INTEGER, ALLOCATABLE :: n_nodes(:), displacements(:) !Number of particles per process and displacement array
  TYPE(MPI_Datatype) :: MPI_PARTICLE3D !MPI 3D particle

  character(len=30) :: data    !Input file containing initial particle data.
  character(len=30) :: orbits  !Output file for storing particle positions over time.

  data = 'data_input.dat'      !Assign input file name.
  orbits = 'data_output.dat'   !Assign output file name.



  !Initialize the MPI environment
  CALL Initialize_MPI(comsize, rank, ierr)

  !Create an MPI datatype for particle3D type
  CALL Create_MPI_PARTICLE3D(MPI_PARTICLE3D, ierr)

  IF (rank .EQ. 0) THEN 
     
     ! Read the simulation parameters and initial conditions from the input file.
     call system_clock(count_rate=rate)
     call system_clock(count=start)

     OPEN(10, file=data, status='old', action='read')

     ! Read the simulation time step, output time step, total time, and number of particles.
     READ(10, *) dt
     PRINT *, "The read value of the time step (dt) is ", dt

     READ(10, *) dt_out
     PRINT *, "The read value of the output time step (dt_out) is", dt_out

     READ(10, *) t_end
     PRINT *, "The read value of the total simulation time (t_end) is ", t_end

     READ(10, *) n
     PRINT *, "The read value of the number of bodies (n) is ", n

     ALLOCATE(a(n))  !Allocate memory for acceleration vectors.
     ALLOCATE(p(n))  !Allocate memory for particle structures.

     ! Read positions, velocities, and masses of each particle.
     DO i = 1, n
        READ(10, *) p(i)%m, p(i)%p%x, p(i)%p%y, p(i)%p%z, &            
             p(i)%v%x, p(i)%v%y, p(i)%v%z  
        PRINT *, "Particle", i, ":"
        PRINT *, "  Position:", p(i)%p
        PRINT *, "  Velocity:", p(i)%v
        PRINT *, "  Mass:    ", p(i)%m
     END DO

     CLOSE(10)


  END IF

  !Send simulation parameters from the input file to all processes
  CALL send_data(n, dt, dt_out, t_end, theta, ierr)

  IF (rank .NE. 0) THEN
     ALLOCATE(p(n)) !Allocate memory for particles array for all processes
  END IF

  !Send particle data to all processes
  CALL send_particles(p, MPI_PARTICLE3D, ierr)

  
  ! Initialize the Head Node of the Octree
  ALLOCATE(head)                 
  CALL Calculate_ranges(head, p) !Calculate the spatial range of the octree.
  head%type = 0                  !Initialize the head cell as empty.
  CALL Nullify_Pointers(head)    !Set all subcell pointers to null.

  ! Construct the Initial Octree
  DO i = 1, n
     CALL Find_Cell(head, temp_cell, p(i))  !Find the cell where the particle belongs.
     CALL Place_Cell(temp_cell, p(i), i)    !Place the particle in the appropriate cell.
  END DO
  CALL Borrar_empty_leaves(head)            !Remove empty leaf cells.
  CALL Calculate_masses(head, p)            !Calculate mass and center of mass for each cell.


  
  CALL particles_nodes(n, comsize, rank, ierr, i_start, i_end, n_local, n_extra, n_nodes, displacements) !Determine how particles are distributed among processes
  
  do i = i_start, i_end
     a(i) = vector3d(0.0, 0.0, 0.0)
  end do
  
  CALL gather_particles(p, i_start, i_end, n_nodes, displacements, MPI_PARTICLE3D, rank, ierr) !Gather particle data from all processes

  CALL Calculate_forces(head, n, p, rji, a, i_start, i_end, rank) !Each process calculates forces between its particles based on the tree
  CALL gather_particles(p, i_start, i_end, n_nodes, displacements, MPI_PARTICLE3D, rank, ierr) !Gather particle data from all processes

  
  IF (rank .EQ. 0) THEN
     OPEN(11, file=orbits, status='old', action='write')
  END IF

  ! Main Simulation Loop
  t_out = 0.0
  t = 0.0
  DO WHILE (t <= t_end)
     DO i = i_start, i_end
        p(i)%v = p(i)%v + a(i) * (dt / 2)
     END DO

     DO i = i_start, i_end
        p(i)%p = p(i)%p + p(i)%v * dt
     END DO

     ! Rebuild the octree after positions have changed.
     CALL Borrar_tree(head)
     CALL Calculate_ranges(head, p)
     head%type = 0
     CALL Nullify_Pointers(head)
     DO i = 1, n
        CALL Find_Cell(head, temp_cell, p(i))
        CALL Place_Cell(temp_cell, p(i), i)
     END DO
     CALL Borrar_empty_leaves(head)
     CALL Calculate_masses(head, p)
     
     do i = i_start, i_end
        a(i) = vector3d(0.0, 0.0, 0.0)
     end do

     CALL gather_particles(p, i_start, i_end, n_nodes, displacements, MPI_PARTICLE3D, rank, ierr) !Gather particle data from all processes

     CALL Calculate_forces(head, n, p, rji, a, i_start, i_end, rank) !Each process calculates forces between its particles based on the tree

     DO i = i_start, i_end
        p(i)%v = p(i)%v + a(i) * (dt / 2)
     END DO

     
     CALL gather_particles(p, i_start, i_end, n_nodes, displacements, MPI_PARTICLE3D, rank, ierr) !Gather particle data from all processes


     ! Write output particle positions.
     t_out = t_out + dt
     IF (t_out >= dt_out) THEN
        WRITE(11, '(E12.2)', ADVANCE='NO') t
        DO i = 1, n
           WRITE(11, '(3X, E12.6, 3X, E12.6, 3X, E12.6)', ADVANCE='NO') & 
               p(i)%p%x, p(i)%p%y, p(i)%p%z
        END DO
        WRITE(11, *)
        t_out = 0.0
     END IF

     t = t + dt
  END DO


  IF (rank .EQ. 0) THEN
     CLOSE(11)
  END IF


  ! Calculate and print the elapsed time of the simulation.
  CALL system_clock(count=finish)
  elapsed_time = real(finish - start, kind=real64) / real(rate, kind=real64)
  IF (elapsed_time >= 60) THEN
     print *, "Elapsed time:", floor(elapsed_time / 60), "min", &
         floor((elapsed_time / 60 - floor(elapsed_time / 60)) * 60), &
         "s"
  ELSE
     PRINT *, "Elapsed time:", elapsed_time, "seconds"
  END IF

END PROGRAM tree
