MODULE init_conds_mod
  IMPLICIT NONE
CONTAINS
  SUBROUTINE create_init_conds(output_name, N, R)

  IMPLICIT NONE
  
  integer, intent(in) :: N
  CHARACTER(len=100) , intent(in):: output_name  !  = "data/random_bodies.txt"

  INTEGER :: I
  INTEGER :: values(1:8), k
  INTEGER, dimension(:), allocatable :: seed
  REAL :: mass,rx,ry,rz, R
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
        ! to re-centre around (-R,R)
        rx = 2*R*rx - R 
        ry = 2*R*ry - R 
        IF ((rx**2 + ry**2) .LE. R**2) EXIT
     END DO
     DO
        CALL random_number(rz)
        rz = 2*R*rz - R
        IF ((rx**2 + ry**2 + rz**2) .LE. R**2) EXIT
     END DO
     ! Now we have created a point with mass, r coordinates (inside sphere of radius R)
     ! and we set initial velocities as 0
     WRITE(20,'(F6.3, 3F20.15, 3I2)') mass, rx,ry,rz, 0,0,0
  END DO
  close(20)
  PRINT*, "Initial conditions created."
END SUBROUTINE CREATE_INIT_CONDS
END MODULE init_conds_mod
