PROGRAM ex1
  USE GEOMETRY
  USE PARTICLE
  IMPLICIT NONE

  ! Loop indices
  INTEGER :: i,j,k
  
  ! Number of particles
  INTEGER :: n
  
  ! Timestep, end time, time loop index,
  ! times at which to print
  DOUBLE PRECISION :: dt, t_end, t, dt_out, t_out
  INTEGER :: time_counter, total_timesteps

  ! Vector related quantities (squared, cubed)
  ! I removed rs because it was not used
  DOUBLE PRECISION :: r2, r3 

  ! We will have n instances of the "particle" type
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: particles 
  
  ! Acceleration arrays
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a

  ! Difference vector
  TYPE(vector3d) :: rji

  ! Read the necessary inputs
  INTEGER :: openstatus_input, openstatus_output, readstatus
  CHARACTER(13) :: datafile="data_read.txt"
  CHARACTER(11) :: resultfile='results.txt'
  
  OPEN (UNIT=1, FILE=datafile, STATUS="old", &
       ACTION="read", POSITION="rewind", &
       IOSTAT=openstatus_input)
  IF (openstatus_input > 0) stop "Cannot open file."
  
  READ (1, *, IOSTAT=readstatus) dt
  READ (1, *, IOSTAT=readstatus) dt_out
  READ (1, *, IOSTAT=readstatus) t_end
  READ (1, '(i5)', IOSTAT=readstatus) n

  ! Make sure it was read correctly and sensibly
  PRINT '(A, F8.4)', "This is the selected time step:", dt
  PRINT '(A, F8.4)', "This is the selected time step for outputs:", dt_out
  PRINT '(A, F8.4)', "This is the selected final time:", t_end
  PRINT '(A, I3)', "This is the selected number of particles:", n
  PRINT*, "" ! Blank space

  ! Calculate the necessary timesteps to reach end
  total_timesteps = t_end/dt
  
  ALLOCATE(particles(n))
  ALLOCATE(a(n)) 

  ! Assign the masses & initial conditions
  DO i = 1, n
     READ (1, *, IOSTAT=readstatus) particles(i)%m, & 
          particles(i)%p, particles(i)%v
     PRINT '(A, I2)', "This is the mass for particle", i
     PRINT '(F3.1)', particles(i)%m

     PRINT '(A, I2)', "This is the initial position for particle", i
     PRINT '(F7.3)', particles(i)%p

     PRINT '(A, I2)', "This is the initial velocity for particle", i
     PRINT '(F7.3)', particles(i)%v
     PRINT*, "" ! Blank space
  END DO

  CLOSE (UNIT=1) ! Close the file
  

  ! Make this big block into a function or subroutine
  ! Set all accelerations at 0 initially
  DO i = 1, n
     a(i)%x = 0.
     a(i)%y = 0.
     a(i)%z = 0.
  END DO
  
  ! For each particle
  DO i = 1,n 
     ! And each of its neighbors...
     DO j = i+1,n 
	! Determine difference vector
	rji = particles(j)%p - particles(i)%p
	! The cube of the distance
        r3 = magnitude(rji)**3

	! Compute the accelerations for i, j. Adding to previous value
	a(i) = a(i) + particles(j)%m * rji / r3
	a(j) = a(j) - particles(i)%m * rji / r3 
    END DO
  END DO

  ! Now, compute the velocities and positions after one timestep
  t_out = 0.0 

  OPEN (UNIT=2, FILE=resultfile, STATUS="replace", &
       ACTION="write", POSITION="rewind", &
       IOSTAT=openstatus_output)
  IF (openstatus_output > 0) stop "Cannot open file."

  
  ! For all needed times
  DO time_counter = 0, total_timesteps
  ! DO t = 0.0, t_end, dt 
     ! Compute velocities and positions for 1st time in the timestep
     particles%v = particles%v + a * (dt/2) 
     particles%p = particles%p + particles%v * dt 

     ! Set all accelerations to 0 again and recompute
     DO i = 1, n
        a(i)%x = 0.
        a(i)%y = 0.
        a(i)%z = 0.
     END DO

     
     ! Same as before
     DO i = 1,n 
	DO j = i+1,n
           rji = particles(j)%p - particles(i)%p
    
           r3 = magnitude(rji)**3

           a(i) = a(i) + particles(j)%m * rji / r3 
	   a(j) = a(j) - particles(i)%m * rji / r3 
	END DO 
     END DO

     
     particles%v = particles%v + a * (dt/2)

     
     ! Update the output time (t_out just keeps track of the progress
     ! until print)
     t_out = t_out + dt 
     ! If t_out is bigger than the increments at which we want output:
     IF (t_out >= dt_out) THEN
        WRITE (2, '(A, F5.2)') "t=", dt*time_counter
	! For each of the particles
	DO i = 1,n 
           ! Print ALL the positions if it is time to do so
           WRITE (2, '(F7.3, F7.3, F7.3)') particles(i)%p%x, &
                particles(i)%p%y, particles(i)%p%z
           WRITE (2, '(A)') ""
	END DO
        t_out = 0.0
     END IF

  END DO

  CLOSE (UNIT=2)

  
  
END PROGRAM ex1
