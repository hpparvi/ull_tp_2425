PROGRAM ex3
  USE geometry !Importing the geometry module
  USE particle !Importing the particle module
  USE barnes   !Importing the barnes module, where the Barnes-Hut algorithm is defined
  USE calcs    !Importing the calcs module, where subroutines to update properties of the particles are defined
  USE data     !Importing the data module, where subroutines to read and save data are defined
  USE iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  USE mpi_types !Importing the mpi_type module, where subroutines to use MPI types are defined
  USE mpi_env   !Importing the mpi_env module, where subroutines to work on MPI are defined
  USE mpi_f08   !Importing the mpi library to use MPI 
  IMPLICIT NONE

  INTEGER(INT64) :: i, n !Loop indexing variable and number of bodies
  REAL(REAL64) :: dt, t_end, t, dt_out, t_out, theta !Time step, final time, current time, output time step, output time and parameter that determines the accuracy of the simulation
  TYPE(vector3d) :: rji !Vector from one particle to another
  TYPE(particle3d), ALLOCATABLE :: p(:)  !Particles
  TYPE(CELL), POINTER :: head, temp_cell !Pointers for the root cell and temporary cells in the tree structure

  INTEGER(INT64) :: start, finish, rate !Variables to store the initial, final and frequency values of the system clock counter
  REAL(REAL64) :: elapsed_time !Variable to calculate and store the elapsed time

  INTEGER :: comsize, rank, ierr !Total number of processes, processes IDs and variable for error handling in MPI operations
  INTEGER :: i_start, i_end, n_local, n_extra !First and last index of the range of particles to process, number of particles per process and extra particles
  INTEGER, ALLOCATABLE :: n_nodes(:), displacements(:) !Number of particles per process and displacement array
  TYPE(MPI_Datatype) :: MPI_PARTICLE3D !MPI 3D particle

  !Initialize the MPI environment
  CALL Initialize_MPI(comsize, rank, ierr)

  !Create an MPI datatype for particle3D type
  CALL Create_MPI_PARTICLE3D(MPI_PARTICLE3D, ierr)

  IF (rank .EQ. 0) THEN 
     CALL system_clock(count_rate = rate) !Get the system clock frequency and stored it in the rate variable
     CALL system_clock(count = start)     !Get the current clock counter state and stored it in the start variable

     CALL read_data(n, p, dt, t_end, dt_out, t_out, theta) !Read the input file

     PRINT*, "Performing the simulation with ", comsize, " processes..." !Output message to notify the user that the simulation is being performed
  END IF

  !Send simulation parameters from the input file to all processes
  CALL send_data(n, dt, dt_out, t_end, theta, ierr)

  IF (rank .NE. 0) THEN
     ALLOCATE(p(n)) !Allocate memory for particles array for all processes
  END IF

  !Send particle data to all processes
  CALL send_particles(p, MPI_PARTICLE3D, ierr)
  
  !Initialize head node for the hierarchical tree
  ALLOCATE(head)
  
  CALL Calculate_ranges(head, p) !Calculate the ranges of the particles
  head%type = 0 !Set the initial type of the head cell to 0 (no particles)
  CALL Nullify_Pointers(head) !Nullify all pointers in the head cell

  !Create the initial tree of cells
  DO i = 1, n
     CALL Find_Cell(head, temp_cell, p(i)) !Find the cell for the particle
     CALL Place_Cell(temp_cell, p(i), i)   !Place the particle in the cell
  END DO
  
  CALL Borrar_empty_leaves(head) !Delete subcells with no particles
  CALL Calculate_masses(head, p) !Calculate the masses and centers of mass for the cells

  CALL particles_nodes(n, comsize, rank, ierr, i_start, i_end, n_local, n_extra, n_nodes, displacements) !Determine how particles are distributed among processes
  CALL set_acceleration(p(i_start:i_end)) !Each process resets the acceleration of its particles to zero
  CALL Calculate_forces(head, p, rji, theta, i_start, i_end, rank) !Each process calculates forces between its particles based on the tree

  t = 0.0     !Initialize time
  t_out = 0.0 !Initialize output time

  !Open the output file
  IF (rank .EQ. 0) THEN
     CALL open_output(n)
  END IF
 
  !Main loop to update properties of particles until final time is reached
  DO WHILE (t .LE. t_end)
     CALL velocity(dt, p(i_start:i_end)) !Each process updates the velocities of its particles
     CALL position(dt, p(i_start:i_end)) !Each process updates the positions of its particles
     CALL gather_particles(p, i_start, i_end, n_nodes, displacements, MPI_PARTICLE3D, rank, ierr) !Gather particle data from all processes
     
     !Delete and reinitialize the tree due to updated positions
     CALL Borrar_tree(head) !Delete the old tree structure
     CALL Calculate_ranges(head, p) !Recalculate the ranges of the particles
     head%type = 0  !Set the type of the head cell to 0 (no particles)
     CALL Nullify_Pointers(head) !Nullify all pointers in the head cell
     
     !Rebuild the tree with updated particle positions
     DO i = 1, n
        CALL Find_Cell(head, temp_cell, p(i)) !Find the cell for the particle
        CALL Place_Cell(temp_cell, p(i), i)   !Place the particle in the cell
     END DO
     
     CALL Borrar_empty_leaves(head) !Delete subcells with no particles
     CALL Calculate_masses(head, p) !Recalculate masses and centers of mass for the cells
     
     CALL set_acceleration(p(i_start:i_end)) !Each process resets the acceleration of its particles to zero
     CALL Calculate_forces(head, p, rji, theta, i_start, i_end, rank) !Each process recalculates forces between its particles
     
     CALL velocity(dt, p(i_start:i_end)) !Each process updates the velocities of its particles
     CALL gather_particles(p, i_start, i_end, n_nodes, displacements, MPI_PARTICLE3D, rank, ierr) !Gather particle data from all processes

     t_out = t_out + dt !Update output time

     IF (rank .EQ. 0) THEN 
        IF (t_out .GE. dt_out) THEN
           CALL save_data(n, p, t) !Write current time and position of each particle to the output file
           t_out = 0.0 !Reset output time
        END IF
     END IF

     t = t + dt !Increment the current time by the time step
     
  END DO

  IF (rank .EQ. 0) THEN
     CALL close_output !Close the output file

     CALL system_clock(count = finish)    !Get the current clock counter value and stored it in the finish variable
     elapsed_time = (finish - start)/rate !Calculates the elapsed time in seconds
     
     !Output the elapsed time
     IF (elapsed_time .GT. 60) THEN
        PRINT *, "Time spent on performing the simulation: ", floor(elapsed_time/60), "min", &
             floor((elapsed_time / 60 - floor(elapsed_time / 60)) * 60), &
             "s" !Output in minutes and seconds if the elapsed_time > 60 s
     ELSE
        PRINT *, "Time spent on performing the simulation: ", elapsed_time, "s" !Ouput in seconds if the elapsed_time < 60 s
     END IF

     PRINT*, "Data calculated stored in file: ", output !End of program message
   END IF

   !Deallocate MPI datatypes
   CALL Deallocate_MPI_Types(MPI_PARTICLE3D, ierr)
  
  !Finalize the MPI environment
  CALL Finalize_MPI(ierr)
  
END PROGRAM ex3
