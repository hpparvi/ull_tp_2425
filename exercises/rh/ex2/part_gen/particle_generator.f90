program particle_generator
  use iso_fortran_env
  IMPLICIT NONE
  INTEGER :: I,N
  INTEGER :: values(1:8), k
  INTEGER, dimension(:), allocatable :: seed
  REAL(real64) :: mass,rx,ry,rz
  REAL(real64) :: dt, t_end, t, dt_out! time variables
  CHARACTER(len=*), PARAMETER :: filename = 'particle_position.dat'
  INTEGER :: rc 
  
  OPEN (file = filename, action = 'write', status = 'replace', unit = 4, iostat = rc) 
    IF (rc/=0) WRITE (*,*) 'Cannot open file ' , filename
  
  PRINT*, "dt?"
  READ*, dt
  WRITE(4, *) dt
  
  PRINT*, "dt_out?"
  READ*, dt_out
  WRITE(4, *) dt_out
  
  PRINT*, "t_end?"
  READ*, t_end
  WRITE(4, *) t_end
  
  CALL date_and_time(values=values)
  CALL random_seed(size=k)
  ALLOCATE(seed(1:k))

  SEED(:) = values(8)
  CALL random_seed(put=seed)
  PRINT*, "Number of bodies?"
  READ*, N
  WRITE(4, *) N
  mass = 1.0 / N
  DO I= 1,N
    CALL random_number(rx)
    DO
      CALL random_number(ry)
        IF ((rx**2 + ry**2) .LE. 1) EXIT
    END DO
    DO
      CALL random_number(rz)
      IF ((rx**2 + ry**2 + rz**2) .LE. 1) EXIT
    END DO

    WRITE(4, *) mass, rx,ry,rz, 0,0,0
    PRINT*, "dist", SQRT(rx**2 + ry**2 + rz**2)
  END DO

end program particle_generator
