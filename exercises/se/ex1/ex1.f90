PROGRAM ex1
  IMPLICIT NONE

  ! Loop indices
  INTEGER :: i,j,k
  
  ! Number of particles
  INTEGER :: n
  
  ! Timestep, end time, time loop index,
  ! times at which to print
  DOUBLE PRECISION :: dt, t_end, t, dt_out, t_out

  ! Vector related quantities (squared, cubed)
  ! I removed rs because it was not used
  DOUBLE PRECISION :: r2, r3 

  ! We will have n instances of the "particle" type
  TYPE(particle), DIMENSION(:), ALLOCATABLE :: particles 
  



END PROGRAM ex1
