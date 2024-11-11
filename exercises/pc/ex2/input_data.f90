PROGRAM input_data
  USE iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  IMPLICIT NONE
  
  INTEGER(INT64) :: i, n
  INTEGER :: values(1:8), k
  INTEGER, dimension(:), allocatable :: seed
  REAL(REAL64) :: mass, rx, ry, rz

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
  
  mass = 1.0/n

  !Open the input file and check for errors
  OPEN(12, file = input, status = 'replace', action = 'write', iostat = io_status)
  IF (io_status /= 0) THEN
      PRINT*, "Error opening file: ", input !Print error message if file opening fails
      STOP
   END IF

  !Write the header for the input file
  WRITE(12, "(A5, 7A12)") "Body", "Mass", "x", "y", "z", "vx", "vy", "vz"
  
  DO i = 1, n
     CALL random_number(rx)
     
     DO
        CALL random_number(ry)
        IF ((rx**2 + ry**2) .LE. 1) EXIT
     END DO
     
     DO
        CALL random_number(rz)
        IF ((rx**2 + ry**2 + rz**2) .LE. 1) EXIT
     END DO
     
     WRITE(12, "(I3, 7X, 7ES12.4)") i, mass, rx, ry, rz, 0.0, 0.0, 0.0 !Write the mass and the current position and velocity of each particle to the input file
  END DO

  CLOSE(12)!Close the input file

  PRINT*, "Data calculated stored in file: ", input !End of program message
  
END PROGRAM input_data
