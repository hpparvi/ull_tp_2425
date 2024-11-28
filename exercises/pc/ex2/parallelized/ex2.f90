PROGRAM ex2
  USE geometry !Importing the geometry module
  USE particle !Importing the particle module
  USE barnes   !Importing the barnes module, where the Barnes-Hut algorithm is defined
  USE calcs    !Importing the calcs module, where subroutines to update properties of the particles are defined
  USE data     !Importing the data module, where subroutines to read and save data are defined
  USE iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  !$ USE omp_lib  !Importing the omp library to use OpenMP
  IMPLICIT NONE

  INTEGER(INT64) :: n !Number of bodies
  REAL(REAL64) :: dt, t_end, t, dt_out, t_out, theta !Time step, final time, current time, output time step, output time and parameter that determines the accuracy of the simulation
  TYPE(vector3d) :: rji !Vector from one particle to another
  TYPE(particle3d), allocatable :: p(:) !Particles
  TYPE(CELL), POINTER :: head, temp_cell !Pointers for the root cell and temporary cells in the tree structure

  INTEGER(INT64) :: start, finish, rate !Variables to store the initial, final and frequency values of the system clock counter
  REAL(REAL64) :: elapsed_time !Variable to calculate and store the elapsed time
  
  CALL system_clock(count_rate = rate) !Get the system clock frequency and stored it in the rate variable
  CALL system_clock(count = start) !Get the current clock counter state and stored it in the start variable

  !Read the input file
  CALL read_data(n, p, dt, t_end, t, dt_out, t_out, theta)
  
  !Initializing head node for the hierarchical tree
  ALLOCATE(head)
  
  CALL Calculate_ranges(head, p) !Calculate the ranges of the particles
  head%type = 0 !Set the initial type of the head cell to 0 (no particles)
  CALL Nullify_Pointers(head) !Nullify all pointers in the head cell

  !Creating the initial tree of cells
  DO i = 1, n
     CALL Find_Cell(head, temp_cell, p(i)) !Find the cell for the particle
     CALL Place_Cell(temp_cell, p(i), i)   !Place the particle in the cell
  END DO
  
  CALL Borrar_empty_leaves(head) !Delete subcells with no particles
  CALL Calculate_masses(head, p) !Calculate the masses and centers of mass for the cells

  CALL set_acceleration(n, p) !Reset accelerations of particles to zero
  CALL Calculate_forces(head, n, p, rji, theta) !Calculate forces between particles based on the tree

  t = 0.0 !Initialize time
  t_out = 0.0 !Initialize output time

  !Open the output file
  CALL open_output(n)
 
  !Main loop to update properties of particles until final time is reached
  DO WHILE (t .LE. t_end)
     CALL velocity(n, dt, p) !Update the velocities of particles
     CALL position(n, dt, p) !Update the positions of particles
     
     !Delete and reinitialize the tree due to updated positions
     CALL Borrar_tree(head) !Delete the old tree structure
     CALL Calculate_ranges(head, p) !Recalculate the ranges of the particles
     head%type = 0  !Set the type of the head cell to 0 (no particles)
     CALL Nullify_Pointers(head) !Nullify all pointers in the head cell
     
     !Rebuild the tree with updated particle positions
     DO i = 1, n
        CALL Find_Cell(head, temp_cell, p(i))!Find the cell for the particle
        CALL Place_Cell(temp_cell, p(i), i)  !Place the particle in the cell
     END DO
     
     CALL Borrar_empty_leaves(head) !Delete subcells with no particles
     CALL Calculate_masses(head, p) !Recalculate masses and centers of mass for the cells
     
     CALL set_acceleration(n, p) !Reset accelerations of particles to zero
     CALL Calculate_forces(head, n, p, rji, theta) !Recalculate forces between particles

     CALL velocity(n, dt, p) !Update the velocities of particles
     t_out = t_out + dt      !Update output time

     IF (t_out .GE. dt_out) THEN
        CALL save_data(n, p, t) !Write current time and position of each particle to the output file
        t_out = 0.0 !Reset output time
        
     END IF

     t = t + dt !Increment the current time by the time step
     
  END DO

  CALL close_output !Close the output file

  CALL system_clock(count = finish)    !Get the current clock counter value and stored it in the finish variable
  elapsed_time = (finish - start)/rate !Calculates the elapsed time in seconds
  
  IF (elapsed_time .GT. 60) THEN
     PRINT *, "Time spent on performing the simulation: ", floor(elapsed_time/60), "min", &
         floor((elapsed_time / 60 - floor(elapsed_time / 60)) * 60), &
         "s" !Output the elapsed time
  ELSE
     PRINT *, "Time spent on performing the simulation: ", elapsed_time, "s"
  END IF

  PRINT*, "Data calculated stored in file: ", output !End of program message
  
END PROGRAM ex2
