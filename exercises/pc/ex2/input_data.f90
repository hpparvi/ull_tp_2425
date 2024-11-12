PROGRAM input_data
  USE particle !Importing the particle module 
  USE iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  IMPLICIT NONE
  
  INTEGER(INT64) :: i, n !Number of bodies
  INTEGER(INT32) :: values(1:8), k
  INTEGER(INT32), dimension(:), allocatable :: seed
  TYPE(particle3d), allocatable :: p(:) !Particles
  REAL(REAL64) :: rx, ry, rz, dt, t_end, dt_out !Time step, final time, output time step

  CHARACTER(len = 50) :: input  !Input file name for initial data
  INTEGER :: io_status !Variable to check the status of I/O operations

  input = 'input.dat' !Assigning input file name
  
  CALL date_and_time(values = values)
  CALL random_seed(size = k)
  
  ALLOCATE(seed(1:k))
  
  SEED(:) = values(8)
  CALL random_seed(put = seed)
  
  PRINT*, "Insert the number of bodies: "
  READ*, n

  PRINT*, "Insert the time step: "
  READ*, dt

  PRINT*, "Insert the output time step: "
  READ*, dt_out

  PRINT*, "Insert the final time: "
  READ*, t_end

  ALLOCATE(p(n))
  
  p%m = 1.0/n

  !Open the input file and check for errors
  OPEN(12, file = input, status = 'replace', action = 'write', iostat = io_status)
  IF (io_status /= 0) THEN
      PRINT*, "Error opening file: ", input !Print error message if file opening fails
      STOP
   END IF

  !Write header for the input file
  WRITE(12, "(A9, 6X, 2A16, 5X, A18)") "Time step", "Output time step", "Final time", "Number of bodies"
  !Write the time step, output time step, final time and number of bodies to the input file
  WRITE(12, "(ES12.6, ES15.6, 7X, ES15.6, I14)") dt, dt_out, t_end, n

  WRITE(12, *) !Write blank line
  
  !Write header for the input file
  WRITE(12, "(A5, 7A12)") "Body", "Mass", "x", "y", "z", "vx", "vy", "vz"
  
  DO i = 1, n
     CALL random_number(rx)
     p(i)%p%x = rx
     
     DO
        CALL random_number(ry)
        p(i)%p%y = ry
        IF ((rx**2 + ry**2) .LE. 1) EXIT
     END DO
     
     DO
        CALL random_number(rz)
        p(i)%p%z = rz
        IF ((rx**2 + ry**2 + rz**2) .LE. 1) EXIT
     END DO
     
     WRITE(12, "(I3, 7X, 7ES12.4)") i, p(i)%m, rx, ry, rz, 0.0, 0.0, 0.0 !Write the mass and the current position and velocity of each particle to the input file
  END DO

  CLOSE(12)!Close the input file

  PRINT*, "Initial data stored in file: ", input !End of program message
  
END PROGRAM input_data
