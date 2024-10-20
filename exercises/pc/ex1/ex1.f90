PROGRAM ex1
  USE calcs !only using the module calcs because it's using the geometry and particle modules as well
  IMPLICIT NONE
  INTEGER :: n !n is the number of bodies
  REAL :: dt, t_end, t, dt_out, t_out !time step, final time, time, results printing time step, results printing time, position squared and cubed
  TYPE(vector3d) :: rji !position between two bodies
  TYPE(particle3d), allocatable :: p(:)

  CHARACTER(len = 20) :: input
  CHARACTER(len = 20) :: output
  INTEGER :: io_status

  input = 'initial_data.dat'
  output = 'final_data.dat'

  OPEN(12, file = input, status = 'old', action = 'read', iostat = io_status)
  
  READ(12, *, IOSTAT = io_status) dt, dt_out, t_end, n
  PRINT*, "Time step: ", dt
  PRINT*, "Results printing time interval: ", dt_out
  PRINT*, "Final time: ", t_end
  PRINT*, "Number of bodies: ", n

  !Assigning memory for particles' properties
  ALLOCATE(p(n))

  DO i = 1, n
     READ(12, *, IOSTAT = io_status) p(i)%m, p(i)%p, p(i)%v !reading particles' masses, positions and velocities
     PRINT*, "Body ", i
     PRINT*, "Mass: ", p(i)%m
     PRINT*, "Position (x, y, z): ", p(i)%p
     PRINT*, "Velocity (x, y, z): ", p(i)%v    
  END DO

  CALL set_acceleration(n, p) !calling subroutine that sets the initial particles' accelerations to zero

  CLOSE(12)

  OPEN(13, file = output, status = 'old', action = 'write', iostat = io_status)
  WRITE(13, '(A)') "Position (x, y, z)"
  
  CALL acceleration(n, dt, p, rji) !Calling subroutine to calculate particles' acceleration 

  t = 0.0
  t_out = 0.0

  DO WHILE (t .LE. t_end)
     
     CALL velocity(n, dt, p)
     CALL position(n, dt, p)
     CALL set_acceleration(n, p)
     CALL acceleration(n, dt, p, rji)
     CALL velocity(n, dt, p)
     
     t_out = t_out + dt 
     
     IF (t_out .GE. dt_out) THEN
        DO i = 1, n
           WRITE(13, '(3F12.6)') p(i)%p
        END DO
        
        t_out = 0.0
        
     END IF

     t = t + dt
     
  END DO

  CLOSE(13)

END PROGRAM ex1
