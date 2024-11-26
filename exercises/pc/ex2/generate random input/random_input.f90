PROGRAM random_input
  USE geometry !Importing the geometry module
  USE particle !Importing the particle module 
  USE iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  IMPLICIT NONE
  
  INTEGER(INT64) :: i, n !Integer variable for looping indexing  and number of bodies
  INTEGER(INT32) :: values(1:8), k !Variables for storing data and time info and holding the size of the random seed
  INTEGER(INT32), allocatable :: seed(:) !Random seed
  TYPE(particle3d), allocatable :: p(:) !Particles
  REAL(REAL64) :: rx, ry, rz, dt, t_end, dt_out, theta !Random coordinates for the position of particles (x, y, z), time step, final time, output time step and parameter that
                                                       !determines the accuracy of the simulation

  CHARACTER(len = 50) :: input  !Input file name for initial data
  INTEGER :: io_status !Variable to check the status of I/O operations

  !Read the input file name from the user
  PRINT*, "Insert the input file name: "
  READ*, input
  
  CALL date_and_time(values = values) !Get the current date and time
  CALL random_seed(size = k) !Get the size of the random seed

  !Allocate memory for the random seed
  ALLOCATE(seed(1:k))
  
  SEED(:) = values(8) !Initialize the random seed using the current date and time
  CALL random_seed(put = seed) !Set the random seed

  !Read the number of bodies from the user
  PRINT*, "Insert the number of bodies: "
  READ*, n

  !Read the time step from the user
  PRINT*, "Insert the time step: "
  READ*, dt
  
  !Read the output time step from the user
  PRINT*, "Insert the output time step: "
  READ*, dt_out

  !Read the final time from the user
  PRINT*, "Insert the final time: "
  READ*, t_end

  !Read theta from the user
  PRINT*, "Insert theta: "
  READ*, theta
  
  !Allocate memory for particles array
  ALLOCATE(p(n))

  !Initialize the mass of each particle
  p%m = 1.0/n

  !Open the input file and check for errors
  OPEN(12, file = input, status = 'replace', action = 'write', iostat = io_status)
  IF (io_status /= 0) THEN
      PRINT*, "Error opening file: ", input !Print error message if file opening fails
      STOP
   END IF

  WRITE(12, "(A9, 6X, 2A16, 5X, A18, 5X, A12)") "Time step", "Output time step", "Final time", "Number of bodies", "Theta" !Write header to the input file
  WRITE(12, "(ES12.6, ES15.6, 7X, ES15.6, I14, 15X, ES12.6)") dt, dt_out, t_end, n, theta !Write the time step, output time step, final time, number of bodies and theta to the input file

  WRITE(12, *) !Write a blank line
  
  !Write header for particles properties to the input file
  WRITE(12, "(A5, 7A12)") "Body", "Mass", "x", "y", "z", "vx", "vy", "vz"
  
  DO i = 1, n
     CALL random_number(rx) !Generate a random number for the x-coordinate
     p(i)%p%x = rx
     
     DO
        CALL random_number(ry) !Generate a random number for the y-coordinate
        p(i)%p%y = ry

        IF ((NORM(vector3d(rx, ry, 0))**2) .LE. 1) EXIT !Ensure the generated coordinates are within a circle in the xy-plane
     END DO
     
     DO
        CALL random_number(rz) !Generate a random number for the z-coordinate
        p(i)%p%z = rz

        IF ((NORM(vector3d(rx, ry, rz)))**2 .LE. 1) EXIT !Ensure the generated coordinates are within a sphere
     END DO
     

     WRITE(12, "(I3, 7X, 7ES12.4)") i, p(i)%m, rx, ry, rz, 0.0, 0.0, 0.0 !Write the index, mass, position and velocity of each particle to the input file
  END DO

  CLOSE(12) !Close the input file

  PRINT*, "Initial data stored in file: ", input !End of program message
  
END PROGRAM random_input
