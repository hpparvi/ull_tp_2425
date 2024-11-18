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
  
  !to count the number of arguments. If there is one, use as the number of  particles.
  !If there is  none, ask for the number of particles to be created
  INTEGER:: num_args, ios
  CHARACTER(len=100) :: n_string

  CHARACTER(len=100) :: dummy_string


  CALL date_and_time(values=values)
  CALL random_seed(size=k)
  
  ALLOCATE(seed(1:k))
  SEED(:) = values(8)
  CALL random_seed(put=seed)

  print *, "Creating random initial conditions with ", N, "bodies."

 ! num_args = command_argument_count()
 ! if (num_args<1) then
 !    print *, "No number of bodies provided to generate initial conditions."
 !    print *, "You can skip this step by stating the number of bodies when  exectuing:"
 !    print *, "  ./program <number_bodies>"
 !    PRINT*, "You can also change the output file name by executing ./program <number_bodies> <output_name.txt>"
 !    PRINT *,""
 !    PRINT*, "Type the number of bodies desired."
 !    READ*, N
 ! else
 !    call get_command_argument(1, n_string)
 !    read(n_string, *, iostat=ios) N

 !    PRINT*, N, " bodies distributed randomly through a sphere."
 ! end if

  !if (num_args>1) then
 !    CALL get_command_argument(2, dummy_string)
 !!    print *, "Changing output file name from ", output_name, " to ", dummy_string
 !    CALL get_command_argument(2, output_name)
 ! end if
  
  mass = 1.0 / N

  ! Open the output file
  open(unit=20, file=output_name, status="replace", action="write", iostat=ios)
     if (ios /= 0) then
        print *, "Error: Could not open file ", output_name
        stop
     end if
     
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
     WRITE(20,'(F6.3, 3F20.15, 3I2)') mass, rx,ry,rz, 0,0,0
  END DO
  close(20)
  PRINT*, "Initial conditions created."
END SUBROUTINE CREATE_INIT_CONDs
END MODULE init_conds_mod
