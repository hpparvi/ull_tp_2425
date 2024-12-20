PROGRAM particulas
  use iso_fortran_env
  use geometry
  use particle
  IMPLICIT NONE
  INTEGER :: I, N
  INTEGER :: values(1:8), k
  INTEGER, dimension(:), allocatable :: seed
  INTEGER :: stat !for the error check when opening/closing files
  REAL(REAL64) :: mass, rx, ry, rz, dt, dt_out, t_end
  CHARACTER(LEN=*), parameter :: filename = 'initial_conditions.dat'

  dt = 0.001
  dt_out = 0.1
  t_end = 1.0

  
  CALL date_and_time(values=values)
  CALL random_seed(size=k)
  ALLOCATE(seed(1:k))
  SEED(:) = values(8)
  CALL random_seed(put=seed)
  PRINT*, "Number of bodies?"
  READ*, N
  mass = 1.0 / N
  OPEN (file = filename, action = 'write', status = 'replace', unit = 3, iostat = stat) !opens input file
  IF (stat/=0) WRITE (*,*) 'Cannot open file ' , filename

  WRITE(3, '(F6.3)') dt
  WRITE(3, '(F6.3)') dt_out
  WRITE(3, '(F11.3)') t_end
  
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
     WRITE(3,'(F6.3, 3F11.8,3F11.8)') mass, rx,ry,rz, 0.,0.,0.
     ! PRINT*, "dist", SQRT(rx**2 + ry**2 + rz**2)
  END DO
  CLOSE(3)
END PROGRAM particulas
