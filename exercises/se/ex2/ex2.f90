PROGRAM ex2
  ! USE STATEMENTS
  USE, INTRINSIC :: iso_fortran_env ! for 64-bit reals
  USE geometry
  USE particle
  USE barnes_hut
  IMPLICIT NONE

  ! Loop indices
  INTEGER :: i,j,k

  ! Number of particles
  INTEGER :: n

  ! Time parameters
  REAL(real64) :: dt, t_end, t, dt_out, t_out
  INTEGER :: time_counter, total_timesteps
  
  ! Vector quantities (magnitudes)
  REAL(real64) :: rs, r2, r3

  ! Threshold that defines if a set of particles will be considered as their CM
  REAL(real64), PARAMETER :: theta = 1

  ! Particles
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: particles

  ! Accelerations
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a

  ! Difference vector
  TYPE(vector3d) :: rji

  ! To read the necessary inputs from a file
  INTEGER :: openstatus_input, openstatus_output, readstatus
  CHARACTER(13) :: datafile="input.txt"
  CHARACTER(11) :: resultfile='output.txt'

  TYPE(cell), POINTER :: head, temp_cell

  ! Open the file (already existing)
  OPEN (UNIT=1, FILE=datafile, STATUS="old", &
       ACTION="read", POSITION="rewind", &
       IOSTAT=openstatus_input)
  IF (openstatus_input > 0) stop "Cannot open file."

  ! Read the inputs
  READ (1, *, IOSTAT=readstatus) dt
  READ (1, *, IOSTAT=readstatus) dt_out
  READ (1, *, IOSTAT=readstatus) t_end
  READ (1, '(i5)', IOSTAT=readstatus) n

  ! Make sure it was read correctly and sensibly
  PRINT '(A, F8.4)', "This is the selected time step:", dt
  PRINT '(A, F8.4)', "This is the selected time step for outputs:", dt_out
  PRINT '(A, F8.4)', "This is the selected final time:", t_end
  PRINT '(A, I3)', "This is the selected number of particles:", n
  PRINT*, "" ! Blank space  ! Read the data from the file

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

     PRINT '(A, I2)', "This is the initial position for particle", i
     PRINT '(F7.3)', particles(i)%p
     PRINT '(A, I2)', "This is the initial velocity for particle", i
     PRINT '(F7.3)', particles(i)%v
     PRINT*, "" ! Blank space
  END DO

  CLOSE (UNIT=1) ! Close the file

  ! Initialize the head node
  ALLOCATE(head)

  CALL Calculate_ranges(head, particles) ! space occupied by the head node
  head%type = 0 ! no particle (initialization)
  CALL Nullify_Pointers(head) ! remove all pointers
  

  ! Create the initial tree
  DO i = 1,n ! For all particles
     
     ! Locate the cell where this particle should go
     CALL Find_Cell(head, temp_cell, particles(i))
     
     ! Place the particle inside said cell
     CALL Place_Cell(temp_cell, particles(i), i) ! i is the ID of the particle
                                         ! (called n or pos in the subroutine)
     
  END DO

  
  ! Remove subcells with no particles inside
  CALL Delete_empty_leaves(head)

  ! Calculate the masses (recursively)
  CALL Calculate_masses(head, particles)

  ! Get the initial accelerations (recursive)
  DO i = 1, n
     a(i)%x = 0.
     a(i)%y = 0. 
     a(i)%z = 0.
  END DO
  
  CALL Calculate_forces(head, particles, a)


  ! Open output file
  OPEN (UNIT=2, FILE=resultfile, STATUS="replace", &
       ACTION="write", POSITION="rewind", &
       IOSTAT=openstatus_output)
  IF (openstatus_output > 0) stop "Cannot open file."

  

  ! Main loop
  t_out = 0.0
  DO time_counter = 0, total_timesteps
  ! DO t = 0.0, t_end, dt
     particles%v = particles%v + a * (dt/2) 
     particles%p = particles%p + particles%v * dt 

     ! Delete the tree because the positions have changed
     CALL Delete_tree(head)

     ! Re-generate the tree
     CALL Calculate_ranges(head, particles)
     head%type = 0
     CALL Nullify_Pointers(head)

     ! This updates the tree
     DO i = 1,n
        CALL Find_Cell(head,temp_cell, particles(i))
        CALL Place_Cell(temp_cell, particles(i), i)
     END DO

     ! Remove the leaves with no particles
     CALL Delete_empty_leaves(head)
     
     CALL Calculate_masses(head, particles)

     ! Set all accelerations at 0
     DO i = 1, n
        a(i)%x = 0.
        a(i)%y = 0. 
        a(i)%z = 0.
     END DO
     
     CALL Calculate_forces(head, particles, a)
     particles%v = particles%v + a * (dt/2)
     
     t_out = t_out + dt
     ! If t_out is bigger than the increments at which we want output:
     IF (t_out >= dt_out) THEN
        WRITE (2, '(F10.2)', ADVANCE='no') dt*time_counter
	! For each of the particles
	DO i = 1,n 
           ! Store ALL the positions if it is time to do so
           WRITE (2, '(F15.3, F15.3, F15.3)', ADVANCE='no') particles(i)%p%x, &
                particles(i)%p%y, particles(i)%p%z
	END DO
        WRITE (2, '(A)') "" ! Just to advance to the next line
        t_out = 0.0
     END IF
     
  END DO
  ! End of main loop

  CLOSE (UNIT=2)


  

END PROGRAM ex2
