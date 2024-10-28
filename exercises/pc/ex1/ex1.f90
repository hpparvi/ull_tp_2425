PROGRAM ex1
  USE geometry !Importing the geometry module
  USE particle !Importing the particle module
  USE calcs    !Importing the calcs module, where subroutines to update properties of the particles are defined
  USE iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  IMPLICIT NONE
  
  INTEGER(INT64) :: n !Number of bodies
  REAL(REAL64) :: dt, t_end, t, dt_out, t_out !Time step, final time, current time, output time step and output time
  TYPE(vector3d) :: rji !Vector from one particle to another
  TYPE(particle3d), allocatable :: p(:) !Particles

  CHARACTER(len = 50) :: input  !Input file name for initial data
  CHARACTER(len = 50) :: output !Output file name for final results
  INTEGER :: io_status !Variable to check the status of I/O operations

  input = 'initial_data.dat' !Assigning input file name
  output = 'final_data.dat'  !Assigning output file name

  !Open the input file and check for errors
  OPEN(12, file = input, status = 'old', action = 'read', iostat = io_status)
  IF (io_status /= 0) THEN
      PRINT*, "Error opening file: ", input !Print error message if file opening fails
      STOP
  END IF

  !Read time parameters and number of bodies and check for errors
  READ(12, *, iostat = io_status) !Read blank line
  IF (io_status /= 0) THEN
      PRINT*, "Error reading from file: ", input !Print error message if reading fails
      STOP
  END IF
  
  READ(12, *, iostat = io_status) dt, dt_out, t_end, n !Read time step, output time step, final time and number of bodies
  IF (io_status /= 0) THEN
      PRINT*, "Error reading parameters (dt, dt_out, t_end, n) from file" !Print error message if reading fails
      STOP
  END IF

  PRINT*, "Time step: ", dt !Output the time step
  PRINT*, "Output time step: ", dt_out !Output the output time step
  PRINT*, "Final time: ", t_end !Output the final time
  PRINT*, "Number of bodies: ", n !Output the number of bodies
  PRINT *, "" !Blank line

  !Allocate  memory for particles array
  ALLOCATE(p(n))

  !Read particle properties from the input file
  DO i = 1, n
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

     READ(12, *, iostat = io_status) p(i)%m !Read the mass of each particle
     IF (io_status /= 0) THEN
         PRINT*, "Error reading mass for body ", i !Print error message if reading fails
         STOP
     END IF
     PRINT*, "Body ", i !Output the current body number
     PRINT*, "Mass: ", p(i)%m !Output the mass of the particle

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

     READ(12, *, iostat = io_status) p(i)%p !Read the position of each particle
     IF (io_status /= 0) THEN
         PRINT*, "Error reading position for body ", i !Print error message if reading fails
         STOP
     END IF
     PRINT*, "Position (x, y, z): ", p(i)%p !Output the position of the particle

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

     READ(12, *, iostat = io_status) p(i)%v !Read the velocity of each particle
     IF (io_status /= 0) THEN
         PRINT*, "Error reading velocity for body ", i !Print error message if reading fails
         STOP
     END IF
     PRINT*, "Velocity (x, y, z): ", p(i)%v !Output the velocity of the particle
     PRINT *, "" !Blank line
  END DO

  CLOSE(12) !Close the input file

  CALL set_acceleration(n, p) !Call subroutine that sets the initial accelerations of particles to zero
  CALL acceleration(n, dt, p, rji) !Call subroutine to update acceleration between particles

  t = 0.0 !Initialize time
  t_out = 0.0 !Initialize output time

  !Open the output file and check for errors
  OPEN(13, file = output, status = 'replace', action = 'write', iostat = io_status)
  IF (io_status /= 0) THEN
      PRINT*, "Error opening file: ", output !Print error message if file opening fails
      STOP
  END IF

  !Write the header for the output file
  WRITE(13, "(A5, 6A12, A17)") "Body", "x", "y", "z", "vx", "vy", "vz", "Time"

  !Main loop to update properties of particles until final time is reached
  DO WHILE (t .LE. t_end)
     
     CALL velocity(n, dt, p) !Call subroutine to update the velocities of particles
     CALL position(n, dt, p) !Call subroutine to update the positions of particles
     CALL set_acceleration(n, p) !Call subroutine to reset accelerations of particles to zero
     CALL acceleration(n, dt, p, rji) !Call subroutine to update acceleration between particles based on current positions
     CALL velocity(n, dt, p) !Call subroutine to update the velocities of particles again with new accelerations
     
     t_out = t_out + dt !Increment the output time by the time step
     
     IF (t_out .GE. dt_out) THEN !Output results
        DO i = 1, n
           WRITE(13, "(I3, 7X, 6ES12.4, 1F12.2)") i, p(i)%p, p(i)%v, t !Write the current position and velocity of each particle to the output file
        END DO
        
        t_out = 0.0 !Reset output time
     END IF

     t = t + dt !Increment the current time by the time step
     
  END DO

  CLOSE(13) !Close the output file

  PRINT*, "Data calculated stored in file: ", output !End of program message

END PROGRAM ex1
