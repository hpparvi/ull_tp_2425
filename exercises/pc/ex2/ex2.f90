PROGRAM ex2
  USE geometry !Importing the geometry module
  USE particle !Importing the particle module
  USE barnes   !Importing the barnes module, where the Barnes-Hut algorithm is defined
  USE calcs    !Importing the calcs module, where subroutines to update properties of the particles are defined
  USE iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  IMPLICIT NONE

  INTEGER(INT64) :: n !Number of bodies
  REAL(REAL64) :: dt, t_end, t, dt_out, t_out !Time step, final time, current time, output time step and output time
  TYPE(vector3d) :: rji !Vector from one particle to another
  TYPE(particle3d), allocatable :: p(:) !Particles
  REAL(REAL64) :: dummy !Temporary variable for reading input data
  TYPE(CELL), POINTER :: head, temp_cell !Pointers for the root cell and temporary cells in the tree structure
  
  CHARACTER(len = 50) :: input  !Input file name for initial data
  CHARACTER(len = 50) :: output !Output file name for final results
  CHARACTER(len = 10) :: temp_str !Temporary character variable for the header of the output file
  INTEGER :: io_status !Variable to check the status of I/O operations

  !Read the input file name from the user
  PRINT*, "Insert the input file name: "
  READ*, input

  !Assigning output file name
  IF (input == 'input_angels.dat') THEN
     output = 'output_angels.dat'
     
  ELSE
     output = 'output.dat'
     
  END IF

  !Open the input file and check for errors
  OPEN(12, file = input, status = 'old', action = 'read', iostat = io_status)
  IF (io_status /= 0) THEN
      PRINT*, "Error opening file: ", input !Print error message if file opening fails
      STOP
  END IF
   
  !Read time parameters and number of bodies from the input file and check for errors
  READ(12, *, iostat = io_status) !Read blank line
  IF (io_status /= 0) THEN
      PRINT*, "Error reading from file: ", input !Print error message if reading fails
      STOP
  END IF
  
  READ(12, *, iostat = io_status) dt, dt_out, t_end, n !Read time step, output time step, final time and number of bodies from the input file
  IF (io_status /= 0) THEN
      PRINT*, "Error reading parameters (dt, dt_out, t_end, n) from file" !Print error message if reading fails
      STOP
  END IF

  PRINT*, "Time step: ", dt !Output the time step
  PRINT*, "Output time step: ", dt_out !Output the output time step
  PRINT*, "Final time: ", t_end   !Output the final time
  PRINT*, "Number of bodies: ", n !Output the number of bodies
  PRINT *, "" !Blank line

  !Allocate memory for particles array
  ALLOCATE(p(n))

  READ(12, *, iostat = io_status) !Read blank line
  IF (io_status /= 0) THEN
      PRINT*, "Error reading from file: ", input !Print error message if reading fails
      STOP
  END IF

  READ(12, *, iostat = io_status) !Read blank line
  IF (io_status /= 0) THEN
      PRINT*, "Error reading from file: ", input !Print error message if reading fails
      STOP
  END IF

  !Read particle properties from the input file
  DO i = 1, n
     READ(12, *, iostat = io_status) dummy, p(i)%m, p(i)%p, p(i)%v !Read the mass, position and velocity of each particle

     IF (io_status /= 0) THEN
        PRINT*, "Error reading data for body ", i !Print error message if reading fails
        STOP
     END IF
     
     PRINT*, "Body ", i !Output the current body number
     PRINT*, "Mass: ", p(i)%m !Output the mass of the particle
     PRINT*, "Position (x, y, z): ", p(i)%p !Output the position of the particle
     PRINT*, "Velocity (x, y, z): ", p(i)%v !Output the velocity of the particle
     PRINT *, "" !Blank line
  END DO

  CLOSE(12) !Close the input file

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
  CALL Calculate_forces(head, n, p, rji) !Calculate forces between particles based on the tree

  t = 0.0 !Initialize time
  t_out = 0.0 !Initialize output time

  !Open the output file and check for errors
  OPEN(13, file = output, status = 'replace', action = 'write', iostat = io_status)
  IF (io_status /= 0) THEN
      PRINT*, "Error opening file: ", output !Print error message if file opening fails
      STOP
  END IF

  !Write header for the output file
  WRITE(13, "(A12)", ADVANCE = 'no') "Time"  !Column header for time 

  !Column headers for positions of the particles
  DO i = 1, n
     WRITE(temp_str, "(I0)") i !Convert integer i into a character string
     !Write the position in x, y, z coordinates for particle i
     WRITE(13, "(A12)", ADVANCE = 'no') "x" // TRIM(temp_str) // " "
     WRITE(13, "(A12)", ADVANCE = 'no') "y" // TRIM(temp_str) // " "
     WRITE(13, "(A12)", ADVANCE = 'no') "z" // TRIM(temp_str) // " "
  END DO

  WRITE(13, *) !New line for writing the results
  
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
     CALL Calculate_forces(head, n, p, rji) !Recalculate forces between particles

     CALL velocity(n, dt, p) !Update the velocities of particles
     t_out = t_out + dt      !Update output time

     IF (t_out .GE. dt_out) THEN
        WRITE(13, "(F12.2)", ADVANCE = 'no') t !Write the current time to the output file
        
        DO i = 1, n
           WRITE(13, "(3ES12.4)", ADVANCE = 'no') p(i)%p%x, p(i)%p%y, p(i)%p%z !Write the current position of each particle to the output file
        END DO

        WRITE(13, *) !New line for the next set of results
        
        t_out = 0.0 !Reset output time
     END IF

     t = t + dt !Increment the current time by the time step
     
  END DO

  CLOSE(13) !Close the output file

  PRINT*, "Data calculated stored in file: ", output !End of program message
  
END PROGRAM ex2
