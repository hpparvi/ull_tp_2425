PROGRAM ex1
  USE particle
  IMPLICIT NONE
  REAL(real64) :: dt, t_end, dt_out, t_out
  REAL(real64) :: t
  TYPE(particle3d), ALLOCATABLE :: particles(:)
  INTEGER :: rc, stat, i, j
  INTEGER :: n = 0
  CHARACTER(len=*), PARAMETER :: filename = 'initial_conditions.dat', outname = 'result.dat'
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: aa
  TYPE(vector3d) :: rji
  REAL(real64) :: r2

  

  OPEN (file = filename, action = 'read', status = 'old', unit = 3, iostat = rc)
  IF (rc/=0) WRITE (*,*) 'Cannot open file ' , filename

  DO i = 1, 3, 1 ! skips the first 3 lines as those are not particles
     READ(3, *)
  END DO
  
  DO ! this loop counts the amount of particles
     READ(3, *, iostat = stat)
     IF (stat/=0) EXIT
     n = n + 1
  END DO
  
  REWIND(3)
  
  READ (3, *) dt
  READ (3, *) t_end
  READ (3, *) dt_out

  ALLOCATE(particles(n))

  DO i = 1, n
     READ (3, *) particles(i)%m, particles(i)%p%xx, particles(i)%p%yy, particles(i)%p%zz,&
          &particles(i)%v%xx, particles(i)%v%yy, particles(i)%v%zz
  END DO   
  CLOSE(3)

  ALLOCATE(aa(n))
  aa = vector3d(0.,0.,0.)
  
  
  CALL calculate_accelerations(particles, aa)

  t_out = 0.0

  t = 0.0

  OPEN (file = outname, action = 'write', status = 'replace', unit = 4, iostat = rc)
  IF (rc/=0) WRITE (*,*) 'Cannot open file ' , outname

  DO WHILE (t .LT. t_end)

     particles%v = particles%v + aa * dt/2
     particles%p = particles%p + particles%v * dt
     aa = vector3d(0., 0., 0.)
     CALL calculate_accelerations(particles, aa)

     particles%v = particles%v + aa * dt/2
     t_out = t_out + dt

     IF (t_out >= dt_out) THEN
        DO i = 1, n
           WRITE(4, *) particles(i)%p
        END DO
     t_out = 0.0
     END IF

     t = t + dt
   END DO

  CLOSE(4)
  
  

CONTAINS
  
  SUBROUTINE calculate_accelerations(bodies, accelerations)
    TYPE(particle3d), DIMENSION(:), INTENT(in) :: bodies
    TYPE(vector3d), DIMENSION(:), INTENT(inout) :: accelerations

    DO i = 1, n
       DO j = i+1, n
          rji =  normalize(bodies(j)%p - bodies(i)%p)
          r2 = distance(bodies(j)%p, bodies(i)%p)
          accelerations(i) = accelerations(i) + bodies(j)%m * rji/r2
          accelerations(j) = accelerations(j) - bodies(i)%m * rji/r2
       END DO
    END DO

  END SUBROUTINE calculate_accelerations


END PROGRAM ex1
