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

  ! Vector related quantities (squared, cubed)
  ! I removed rs because it was not used
  DOUBLE PRECISION :: r2, r3 

  ! We will have n instances of the "particle" type
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: particles 
  
  ! Acceleration arrays
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: a

  ! Difference vector
  TYPE(vector3d) :: rji

  ! Take the necessary inputs
  READ*, dt 
  READ*, dt_out 
  READ*, t_end 
  READ*, n 

  ALLOCATE(particles(n))
  ALLOCATE(a(n,3)) 

  ! Assign the masses & initial conditions
  DO i = 1, n 
	READ*, particles(i)%m, particles(i)%p, particles(i)%v 
  END DO


  ! Make this big block into a function or subroutine
  ! Set all accelerations at 0 initially
  a = 0.0 
  ! For each particle
  DO i = 1,n 
     ! And each of its neighbors...
     DO j = i+1,n 
	! Determine difference vector
        ! WARNING: THIS IS NOT OK, %p IS A POINT, - IS FOR VECTOR & POINT
	rji = particles(j)%p - particles(i)%p
	! The cube of the distance
        r3 = magnitude(rji)**3

	! Compute the accelerations for i, j. Adding to previous value
	a(i,:) = a(i,:) + particles(j)%m * rji / r3 
	a(j,:) = a(j,:) - particles(i)%m * rji / r3 
    END DO
  END DO

  ! Now, compute the velocities and positions after one timestep
  t_out = 0.0 

  ! For all needed times
  DO t = 0.0, t_end, dt 
     ! Compute velocities and positions for 1st time in the timestep
     particles%v = particles%v + a * dt/2 
     particles%p = particles%p + particles%v * dt 

     ! Set all accelerations to 0 again and recompute
     a = 0.0 
     ! Same as before
     DO i = 1,n 
	DO j = i+1,n
           ! WARNING: THIS IS NOT OK, %p IS A POINT, - IS FOR VECTOR & POINT
           rji = particles(j)%p - particles(i)%p
    
           r3 = magnitude(rji)**3

           a(i,:) = a(i,:) + particles(j)%m * rji / r3 
	   a(j,:) = a(j,:) - particles(i)%m * rji / r3 
	END DO 
     END DO

     
     particles%v = particles%v + a * dt/2

     ! Update the output time (t_out just keeps track of the progress
     ! until print)
     t_out = t_out + dt 
     ! If t_out is bigger than the increments at which we want output:
     IF (t_out >= dt_out) THEN 
	! For each of the particles
	DO i = 1,n 
           ! Print ALL the positions if it is time to do so
	   PRINT*, particles(i)%p 
	END DO 
	t_out = 0.0 
     END IF

  END DO


  
  
END PROGRAM ex1
