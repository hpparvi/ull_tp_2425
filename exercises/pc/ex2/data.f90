MODULE data
  USE geometry !Importing the geometry module
  USE particle !Importing the particle module
  USE barnes   !Importing the barnes module, where the Barnes-Hut algorithm is defined
  USE calcs    !Importing the calcs module, where subroutines to update properties of the particles are defined
  USE iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  IMPLICIT NONE

  CHARACTER(len = 100) :: input, identifier !Input file name for initial data and identifier of the input file
  CHARACTER(len = 100) :: output !Output file name for final results
  INTEGER :: io_status !Variable to check the status of I/O operations

CONTAINS
  
  !Subroutine to read data from an input file
  SUBROUTINE read_data(n, p, dt, t_end, t, dt_out, t_out, theta)
    INTEGER(INT64), INTENT(INOUT) :: n !Number of bodies
    REAL(REAL64), INTENT(INOUT) :: dt, t_end, dt_out, t_out, theta !Time step, final time, output time step, output time and parameter that determines the accuracy of the simulation       
    INTEGER(INT64) :: i !Loop indexing variable
    TYPE(particle3d), allocatable, INTENT(INOUT) :: p(:) !Particles

    REAL(REAL64) :: dummy !Temporary variable for reading input data
    INTEGER(INT64) :: input_pos, dot_pos !Variables for 'input' and '.' positions

    !Read the input file name from the user
    PRINT*, "Insert the input file name: "
    READ*, input

    !Checking if input file name contains 'input'
    input_pos = INDEX(input, 'input')
    dot_pos = INDEX(input, '.')

    IF (input_pos /= 0 .AND. dot_pos .GT. input_pos) THEN
       identifier = input(input_pos + 5 : dot_pos - 1) !Extracting the identifier after 'input' and before '.'
       output = 'output' // TRIM(identifier) // '.dat' !Output file name with the identifier

    ELSE
       output = 'output.dat' !Output file name in other case
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

    READ(12, *, iostat = io_status) dt, dt_out, t_end, n, theta !Read time step, output time step, final time, number of bodies and theta from the input file
    IF (io_status /= 0) THEN
       PRINT*, "Error reading parameters (dt, dt_out, t_end, n, theta) from file" !Print error message if reading fails
       STOP
    END IF

    PRINT*, "Time step: ", dt !Output the time step
    PRINT*, "Output time step: ", dt_out !Output the output time step
    PRINT*, "Final time: ", t_end   !Output the final time
    PRINT*, "Number of bodies: ", n !Output the number of bodies
    PRINT*, "Theta: ", theta !Output the parameter that determines the accuracy of the simulation
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

  END SUBROUTINE read_data

  !Subroutine to open an output file
  SUBROUTINE open_output(n)
    INTEGER(INT64), INTENT(INOUT) :: n !Number of bodies
    INTEGER(INT64) :: i !Loop indexing variable

    CHARACTER(len = 10) :: temp_str !Temporary character variable for the header of the output file

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

  END SUBROUTINE open_output

  !Subroutine to save data to an output file
  SUBROUTINE save_data(n, p, t)
    INTEGER(INT64), INTENT(INOUT) :: n !Number of bodies
    REAL(REAL64), INTENT(INOUT) :: t   !Current time
    INTEGER(INT64) :: i !Loop indexing variable
    TYPE(particle3d), allocatable, INTENT(INOUT) :: p(:) !Particles

    WRITE(13, "(F12.2)", ADVANCE = 'no') t !Write the current time to the output file

    DO i = 1, n
       WRITE(13, "(3ES12.4)", ADVANCE = 'no') p(i)%p%x, p(i)%p%y, p(i)%p%z !Write the current position of each particle to the output file
    END DO

    WRITE(13, *) !New line for the next set of results

  END SUBROUTINE save_data
    
  !Subroutine to close an output file
  SUBROUTINE close_output
    CLOSE(13) !Close the output file
        
  END SUBROUTINE close_output
  
END MODULE data

