PROGRAM ex1
  USE, INTRINSIC :: iso_fortran_env
  USE GEOMETRY
  USE PARTICLE
  IMPLICIT NONE

  ! Loop indices
  INTEGER :: i,j,k
  
  ! Number of particles
  INTEGER :: n
  
  ! Timestep, end time, time loop indices,
  ! times at which to print
  REAL(real64) :: dt, t_end, t, dt_out, t_out
  INTEGER :: time_counter, total_timesteps

  ! Vector related quantities (magnitude cubed)
  ! I removed rs because it was not used
  REAL(real64) :: r3 

  ! We will have n instances of the "particle" type
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: particles 
  
  ! Acceleration arrays
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a

  ! Difference vector
  TYPE(vector3d) :: rji

  ! To read the necessary inputs from a file
  INTEGER :: openstatus_input, openstatus_output, readstatus
  CHARACTER(20) :: datafile
  CHARACTER(14) :: resultfile='output_ex1.txt'

  ! To test time of execution
  INTEGER :: count_begin, count_end, count_rate

  PRINT '(A)', "Type in the name of the input file with &
       the initial conditions. Must be in the working directory."
  READ (*,*) datafile

  ! Open the file (already existing)
  OPEN (UNIT=1, FILE=datafile, STATUS="old", &
       ACTION="read", POSITION="rewind", &
       IOSTAT=openstatus_input)
  IF (openstatus_input > 0) stop "Cannot open file. Make sure your &
       input data file is named with < 20 characters."

  ! Read the inputs
  READ (1, *, IOSTAT=readstatus) dt
  READ (1, *, IOSTAT=readstatus) dt_out
  READ (1, *, IOSTAT=readstatus) t_end
  READ (1, '(i5)', IOSTAT=readstatus) n

  ! Make sure it was read correctly and sensibly
  PRINT '(A, F8.4)', "This is the selected time step:", dt
  PRINT '(A, F8.4)', "This is the selected time step for outputs:", dt_out
  PRINT '(A, F8.4)', "This is the selected final time:", t_end
  PRINT '(A, I4)', "This is the selected number of particles:", n
  PRINT*, "" ! Blank space

  ! Calculate the necessary timesteps to reach end
  total_timesteps = t_end/dt

  ! Allocate arrays once the dimension is known
  ALLOCATE(particles(n))
  ALLOCATE(a(n))

  ! Assign the masses & initial conditions
  DO i = 1, n
     READ (1, *, IOSTAT=readstatus) particles(i)%m, & 
          particles(i)%p, particles(i)%v
     PRINT '(A, I2)', "This is the mass for particle", i
     PRINT '(F3.1)', particles(i)%m

     PRINT '(A, I4)', "This is the initial position for particle", i
     PRINT '(F7.3)', particles(i)%p
     PRINT '(A, I4)', "This is the initial velocity for particle", i
     PRINT '(F7.3)', particles(i)%v
     PRINT*, "" ! Blank space
  END DO

  CLOSE (UNIT=1) ! Close the file

  ! Start measuring time here
  CALL SYSTEM_CLOCK(count_begin, count_rate)
  ! (If you want:)
  ! PRINT '(A, I8)', "The count rate is ", count_rate
  

  ! Set all accelerations at 0 initially
  a = vector3d(0., 0., 0.)
  
  ! For each particle
  CALL calculate_accelerations(particles, n, a)

  ! Now, compute the velocities and positions after one timestep
  t_out = 0.0 

  ! Open the file for the outputs
  OPEN (UNIT=2, FILE=resultfile, STATUS="replace", &
       ACTION="write", POSITION="rewind", &
       IOSTAT=openstatus_output)
  IF (openstatus_output > 0) stop "Cannot open file."

  
  ! For all needed times
  DO time_counter = 0, total_timesteps
     ! Compute velocities and positions for 1st time in the timestep
     particles%v = particles%v + a * (dt/2) 
     particles%p = particles%p + particles%v * dt 

     ! Set all accelerations to 0 again and recompute
     a = vector3d(0., 0., 0.)

     
     ! Same as before
     CALL calculate_accelerations(particles, n, a)

     ! Update velocity once again
     particles%v = particles%v + a * (dt/2)

     
     ! Update the output time (t_out just keeps track of the progress
     ! until print)
     t_out = t_out + dt 
     ! If t_out is bigger than the increments at which we want output:
     IF (t_out >= dt_out) THEN
        WRITE (2, '(F9.2)', ADVANCE='no') dt*time_counter
	! For each of the particles
	DO i = 1,n 
           ! Print ALL the positions if it is time to do so
           WRITE (2, '(F9.3, F9.3, F9.3)', ADVANCE='no') particles(i)%p%x, &
                particles(i)%p%y, particles(i)%p%z
	END DO
        WRITE (2, '(A)') "" ! Just to advance to the next line
        t_out = 0.0
     END IF

  END DO

  CLOSE (UNIT=2)


  ! End measuring time here (is taking into account writing to file)
  CALL SYSTEM_CLOCK(count_end)

  PRINT '(A, F6.3, A)', "The elapsed time is ", &
       REAL(count_end-count_begin) / REAL(count_rate), &
       " seconds"


CONTAINS

  ! Subroutine that calculates accelerations
  SUBROUTINE calculate_accelerations(particles, n, a)
    TYPE(vector3d),  DIMENSION(:) :: a
    INTEGER :: i,j,n
    TYPE(particle3d), DIMENSION(:) :: particles
    TYPE(vector3d) :: rji
    REAL(real64) :: r3

    DO i = 1,n 
      DO j = i+1,n
        rji = particles(j)%p - particles(i)%p

        r3 = magnitude(rji)**3

        a(i) = a(i) + particles(j)%m * rji / r3 

        a(j) = a(j) - particles(i)%m * rji / r3
        
      END DO
    END DO

  END SUBROUTINE calculate_accelerations


  
END PROGRAM ex1
