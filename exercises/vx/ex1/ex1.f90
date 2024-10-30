PROGRAM leapfrog       !This program does the calculations for solving the N-body problem.
  use iso_fortran_env  !This module ensures all variables are defined as 64-bit.
  use geometry         !This module defines 3D vector and point operations for vector3d and point3d types.
  use particle         !This module defines the particle structure with position, velocity, and mass.
  IMPLICIT NONE

  INTEGER(int64) :: i, j, k
  INTEGER(int64) :: n                                 !Number of bodies.
  REAL(real64) :: dt, t_end, t, dt_out, t_out         !Time step, total simulation time, output time step.
  REAL(real64) :: rs, r2, r3                          !Position, squared distance between two points and cubed distance.
  type(vector3d) :: rji                               !Distance between two bodies.
  
  type(particle3d), allocatable :: p(:)               !Particles of the problem.
  type(vector3d), allocatable :: a(:)                 !Their accelerations.


  !Input and output files.
  character(len=30) :: data                          
  character(len=30) :: orbits

  data = 'data_input.dat'                            
  orbits = 'data_output.dat'

  
  !Read the needed parameters.
  OPEN(10, file=data, status='old', action='read')

  READ(10, *) dt
  print *, "The read value of the time step (dt) is ", dt

  READ(10, *) dt_out
  print *, "The read value of the output time step (dt_out) is", dt_out

  READ(10, *) t_end
  print *, "The read value of the total simulation time (t_end) is ", t_end

  READ(10, *) n
  print *, "The read value of the number of bodies (n) is ", n

  !Allocate memory for particles.
  ALLOCATE(p(n))
  ALLOCATE(a(n))

  !Read positions, velocities, and masses of particles.
  DO i = 1, n
     READ(10, *) p(i)%p
     print *, "The read value of the position of particle", i, "is:", p(i)%p

     READ(10, *) p(i)%v
     print *, "The read value of the velocity of particle", i, "is:", p(i)%v

     READ(10, *) p(i)%m
     print *, "The read value of the mass of particle", i, "is:", p(i)%m
  END DO

  CLOSE(10)


  !Initialize acceleration to zero.
  DO i = 1, n
     a(i) = vector3d(0.0, 0.0, 0.0)
  END DO

  
  OPEN(11, file=orbits, status='old', action='write')  !Open the file where output data will be written.


  !Main code.
  DO i = 1, n
     DO j = i + 1, n
        rji = distance(p(i)%p, p(j)%p)
        r2 = rji%x**2 + rji%y**2 + rji%z**2
        r3 = r2 * sqrt(r2)
        a(i) = a(i) + ((p(j)%m * rji) / r3)
        a(j) = a(j) - ((p(i)%m * rji) / r3)
     END DO
  END DO

  t = 0.0                                                  
  t_out = 0.0
  DO WHILE (t <= t_end)                                    
     DO i = 1, n                                           
        p(i)%v = p(i)%v + (a(i) * (dt / 2))
        p(i)%p = p(i)%p + (p(i)%v * dt)
        a(i) = vector3d(0.0, 0.0, 0.0)
     END DO

     DO i = 1, n
        DO j = i + 1, n
           rji = distance(p(i)%p, p(j)%p)
           r2 = rji%x**2 + rji%y**2 + rji%z**2
           r3 = r2 * sqrt(r2)
           a(i) = a(i) + ((p(j)%m * rji) / r3)
           a(j) = a(j) - ((p(i)%m * rji) / r3)
        END DO
     END DO

     DO i = 1, n                                            
        p(i)%v = p(i)%v + (a(i) * (dt / 2))
     END DO
   
     t_out = t_out + dt
     IF (t_out >= dt_out) THEN
        DO i = 1, n
           WRITE(11, '(E12.6, 3X, E12.6, 3X, E12.6)') p(i)%p%x, p(i)%p%y, p(i)%p%z  !Write output data.
        END DO

        t_out = 0.0
     END IF
     t = t + dt
  END DO

  CLOSE(11)
  
END PROGRAM leapfrog
