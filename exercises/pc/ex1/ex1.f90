PROGRAM ex1
  USE calcs !only using the module calcs because it's using the geometry and particle modules as well
  IMPLICIT NONE
  
  INTEGER(INT64) :: n !n is the number of bodies
  REAL(REAL64) :: dt, t_end, t, dt_out, t_out !time step, final time, time, results printing time step, results printing time, position squared and cubed
  TYPE(vector3d) :: rji !position between two bodies
  TYPE(particle3d), allocatable :: p(:)

  CHARACTER(len = 50) :: input
  CHARACTER(len = 50) :: output
  INTEGER :: io_status

  input = 'initial_data.dat'
  output = 'final_data.dat'

  OPEN(12, file = input, status = 'old', action = 'read', iostat = io_status)
  READ(12, *, iostat  = io_status)
  READ(12, *, iostat  = io_status) dt, dt_out, t_end, n
  PRINT*, "Time step: ", dt
  PRINT*, "Results printing time interval: ", dt_out
  PRINT*, "Final time: ", t_end
  PRINT*, "Number of bodies: ", n

  !Assigning memory for particles' properties
  ALLOCATE(p(n))
 
  DO i = 1, n
     READ(12, *, iostat  = io_status)
     READ(12, *, iostat  = io_status)
 
     READ(12, *, iostat  = io_status) p(i)%m !reading particles' masses
     PRINT*, "Body ", i
     PRINT*, "Mass: ", p(i)%m
 
     READ(12, *, iostat  = io_status)
     READ(12, *, iostat  = io_status)
  
     READ(12, *, iostat  = io_status) p(i)%p !reading particles' positions
     PRINT*, "Position (x, y, z): ", p(i)%p
  
     READ(12, *, iostat  = io_status)
     READ(12, *, iostat  = io_status)
  
     READ(12, *, iostat  = io_status) p(i)%V !reading particles' velocities
     PRINT*, "Velocity (x, y, z): ", p(i)%v    
  END DO

  CLOSE(12)

  CALL set_acceleration(n, p) !calling subroutine that sets the initial particles' accelerations to zero
  CALL acceleration(n, dt, p, rji) !Calling subroutine to calculate particles' acceleration

  t = 0.0
  t_out = 0.0

  OPEN(13, file = output, status = 'replace', action = 'write', iostat = io_status)
  WRITE(13, "(A5, 6A12, A17)") "Body", "x", "y", "z", "vx", "vy", "vz", "Time"

  DO WHILE (t .LE. t_end)
     
     CALL velocity(n, dt, p)
     CALL position(n, dt, p)
     CALL set_acceleration(n, p)
     CALL acceleration(n, dt, p, rji)
     CALL velocity(n, dt, p)
     
     t_out = t_out + dt 
     
     IF (t_out .GE. dt_out) THEN
        DO i = 1, n
           WRITE(13, "(I3, 7X, 6ES12.4, 1F12.2)") i, p(i)%p, p(i)%v, t

           
        END DO
        
        t_out = 0.0
        
     END IF

     t = t + dt
     
  END DO

  CLOSE(13)

  PRINT*, "Data calculated stored in file final_data.dat"

END PROGRAM ex1
