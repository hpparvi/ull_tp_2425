program particulas

  use iso_fortran_env
  IMPLICIT NONE
  
  INTEGER(int64) :: I,N
  INTEGER :: values(1:8), k
  INTEGER, dimension(:), allocatable :: seed
  REAL(real64) :: mass,rx,ry,rz
  real(real64) :: vx, vy, vz
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
  WRITE(10, *) 0.001    ! dt
  WRITE(10, *) 0.1      ! dt_out
  WRITE(10, *) 100.0    ! t_end
  WRITE(10, *) n        ! Número de cuerpos
  
  DO I= 1,N
     CALL random_number(rx)
     DO
        rx = rx * 100
        CALL random_number(ry)
        ry = ry * 100
        IF ((rx**2 + ry**2) .LE. 1) EXIT
     END DO
     
     DO
        rz = rz * 100
        CALL random_number(rz)
        IF ((rx**2 + ry**2 + rz**2) .LE. 1) EXIT
     END DO

     CALL random_number(vz)
     vz = vz/10
     DO
        CALL random_number(vy)
        vy = vy/10
        IF ((vz**2 + vy**2) .LE. 1) EXIT
     END DO

     DO
        CALL random_number(vx)
        vx = vx/10
        IF ((vx**2 + vy**2 + vz**2) .LE. 1) EXIT
     END DO
     
     WRITE(10,'(F6.3, 6F11.8)') mass, rx,ry,rz, vx,vy,vz
     ! PRINT*, "dist", SQRT(rx**2 + ry**2 + rz**2)
  END DO

  close(10)
  
END program particulas
