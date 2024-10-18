PROGRAM ex1
  USE particle !only using the module particle because it's using the module geometry as well
  IMPLICIT NONE
  INTEGER :: i, j, k, n !n is the number of bodies
  REAL :: dt, t_end, t, dt_out, t_out, rs, r2, r3 !time step, final time, time, final time step, out time, ...
  TYPE(vector3d) :: rji
  TYPE(particle3d), allocatable :: p(:)
  
  PRINT*, "Introduce the time step"
  READ*, dt
  PRINT*, "Introduce the ... time step"
  READ*, dt_out
  PRINT*, "Introduce the final time"
  READ*, t_end
  PRINT*, "Introduce the number of bodies for the simulation"
  READ*, n 

  !Assigning memory for particles' properties
  ALLOCATE(p(n))

  DO i = 1, n
     PRINT*, "Insert particle's mass, initial position and velocity"
     READ*, p(i)%m, p(i)%p, p(i)%v !Reading particles' masses, positions and velocities
     
     p(i)%a = vector3d(0.0, 0.0, 0.0) !Setting the initial particles' accelerations to zero
    
  END DO
  
  
  DO i = 1, n
     DO j = i + 1, n
        rji = distance(p(j)%p, p(i)%p)
        r2 = (norm(rji))**2
        r3 = r2 * SQRT(r2)
        p(i)%a = p(i)%a + ((p(j)%m * rji)/r3)
        p(j)%a = p(j)%a - ((p(i)%m * rji)/r3)
     END DO
  END DO

  t = 0.0
  t_out = 0.0
  
  DO WHILE (t .LE. t_end)
     DO i = 1, n
     p(i)%v = p(i)%v + (p(i)%a * (dt * 0.5))
     p(i)%p = p(i)%p + (p(i)%v *  dt)
     p(i)%a = vector3d(0.0, 0.0, 0.0)
     END DO
  
     DO i = 1, n
        DO j = i + 1, n
            rji = distance(p(j)%p, p(i)%p)
            r2 = (norm(rji))**2
            r3 = r2 * SQRT(r2)
            p(i)%a = p(i)%a + ((p(j)%m * rji)/r3)
            p(j)%a = p(j)%a - ((p(i)%m * rji)/r3)
        END DO
     END DO

     DO i = 1, n
        p(i)%v = p(i)%v + (p(i)%a * (dt * 0.5))
     END DO
     
     t_out = t_out + dt 
     
     IF (t_out .GE. dt_out) THEN
        DO i = 1, n
           PRINT*, p(i)%p
        END DO
        
        t_out = 0.0
        
     END IF

     t = t + dt
     
  END DO
END PROGRAM ex1
