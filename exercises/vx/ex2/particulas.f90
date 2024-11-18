program particulas

  use iso_fortran_env
  IMPLICIT NONE
  
  INTEGER(int64) :: I,N
  INTEGER :: values(1:8), k
  INTEGER, dimension(:), allocatable :: seed
  REAL(real64) :: mass,rx,ry,rz
  character(len=30) :: data

  data = 'data_input.dat'

  open(10, file=data, status='replace', action='write')
  
  CALL date_and_time(values=values)
  CALL random_seed(size=k)
  ALLOCATE(seed(1:k))
  SEED(:) = values(8)
  CALL random_seed(put=seed)

  
  PRINT*, "Number of bodies?"
  READ*, N
  mass = 1.0 / N   !Todas las particulas tienen la misma masa?

  
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
     
     WRITE(10,'(F6.3, 3F11.8,3I2)') mass, rx,ry,rz, 0,0,0
     ! PRINT*, "dist", SQRT(rx**2 + ry**2 + rz**2)
  END DO

  close(10)
  
END program particulas
