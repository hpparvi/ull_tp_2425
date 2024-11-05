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





  

END PROGRAM ex2
