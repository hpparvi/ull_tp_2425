PROGRAM init_cond
  IMPLICIT NONE
  
  INTEGER :: I,N
  INTEGER :: values(1:8), k
  INTEGER, dimension(:), allocatable :: seed
  REAL :: mass,rx,ry,rz
  
  CALL date_and_time(values=values)
  CALL random_seed(size=k)
  
  ALLOCATE(seed(1:k))
  SEED(:) = values(8)
  CALL random_seed(put=seed)
  
  PRINT*, "Number of bodies?"
  READ*, N
  mass = 1.0 / N
  DO I= 1,N
     DO
        CALL random_number(rx)
        CALL random_number(ry)
        IF ((rx**2 + ry**2) .LE. 1) EXIT
     END DO
     DO
        CALL random_number(rz)
        IF ((rx**2 + ry**2 + rz**2) .LE. 1) EXIT
     END DO
     ! Now we have created a point with mass, r coordinates (inside sphere of radius 1)
     ! and we set initial velocities as 0
     WRITE(*,’(F6.3, 3F11.8,3I2)’) mass, rx,ry,rz, 0,0,0
  END DO
