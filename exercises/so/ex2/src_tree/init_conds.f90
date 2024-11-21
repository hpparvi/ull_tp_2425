MODULE init_conds_mod
  IMPLICIT NONE
CONTAINS
  SUBROUTINE create_init_conds(output_name, N)

  IMPLICIT NONE
  
  integer, intent(in) :: N
  CHARACTER(len=100) , intent(in):: output_name  !  = "data/random_bodies.txt"

  INTEGER :: I
  INTEGER :: values(1:8), k
  INTEGER, dimension(:), allocatable :: seed
  REAL :: mass,rx,ry,rz
  INTEGER:: ios


  CALL date_and_time(values=values)
  CALL random_seed(size=k)
  
  ALLOCATE(seed(1:k))
  SEED(:) = values(8)
  CALL random_seed(put=seed)

  print *, "Creating random initial conditions with ", N, "bodies."
  
  mass = 1.0 / N

  ! Open the output file
  open(unit=20, file=output_name, status="replace", action="write", iostat=ios)
     if (ios /= 0) then
        print *, "Error: Could not open file ", output_name
        stop
     end if
     
  DO I= 1,N
     DO
        ! uniformly distributed x and y values around (0,1)
        CALL random_number(rx)
        CALL random_number(ry)
        ! to re-centre around (-1,1)
        rx = 2*rx -1 
        ry = 2*ry -1 
        IF ((rx**2 + ry**2) .LE. 1) EXIT
     END DO
     DO
        CALL random_number(rz)
        rz = 2*rz - 1
        IF ((rx**2 + ry**2 + rz**2) .LE. 1) EXIT
     END DO
     ! Now we have created a point with mass, r coordinates (inside sphere of radius 1)
     ! and we set initial velocities as 0
     WRITE(20,'(F6.3, 3F20.15, 3I2)') mass, rx,ry,rz, 0,0,0
  END DO
  close(20)
  PRINT*, "Initial conditions created."
END SUBROUTINE CREATE_INIT_CONDs
END MODULE init_conds_mod
