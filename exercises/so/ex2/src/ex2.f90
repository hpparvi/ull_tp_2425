
PROGRAM leapfrog

  use geometry
  use particles
  use init_conds_mod
  use parameter_reader

  IMPLICIT NONE

  INTEGER :: i,j,k, i_t, ios
  INTEGER :: n, n_t!, N_bodies
  INTEGER :: num_args       ! to count the number of given arguments
  character(len=100) :: par_filename, dummy_line     ! to import the stars data and parameters data
  !character(len=100)  :: !input_file, output_file
  real :: t, t_out !dt, t_end, dt_out
  !logical :: create_bodies
  character(len=100) :: fmt
  
  type(particle), dimension(:), allocatable :: partics
  type(vector3d)                            :: rji_v
  real                                      :: r3_2


  !Expected usage:
  ! - command line: executable <ic_data_filename> <time_param_filename>
  
  !call create_init_conds()
  ! Get the number of command-line arguments
    num_args = command_argument_count()

    if (num_args < 1) then
        print *, "Error: No parameters file provided. Usage: ./program <par_file.par>"
        stop
    end if

    ! Get the input file name from the command-line arguments
    call get_command_argument(1, par_filename)

    ! Remove any trailing spaces from the filename
    par_filename = trim(adjustl(par_filename))

    call read_parameters(par_filename)
    print*, "Reading parameters from ", par_filename
    print*, "..."
    print*, "Parameters loaded:"
    print*, "  -  dt = ", dt
    print*, "  -  dt_out = ", dt_out
    print*, "  -  t_end = ", t_end
    print*, "  -  input_file = ", input_file
    print*, "  -  create_bodies = ", create_bodies
    if (create_bodies .eqv. .true.) then
    print*, "  -  N_bodies = ", N_bodies
 end if
    print*, "  -  output_file = ", output_file

    input_file = trim(input_file)
 

    if (create_bodies .eqv. .true.) then
       call create_init_conds(input_file, N_bodies)
    end if
    
    ! Open the file with initial data
    open(unit=10, file=input_file, status="old", action="read", form="formatted", iostat=ios)
    if (ios /= 0) then
        print *, "Error: Could not open file ", input_file
        stop
     end if
     
    ! Count the number of stars
    n = 0
    do
       read(10, '(A)', iostat=ios) dummy_line  ! reading each line as a whole string
       if (ios /= 0) then
          if (ios == 1) then
          print*, "Error reading ", trim(input_file)
          end if
          if (ios == -1) then
             exit
          end if
          
       end if
       
        n = n + 1
     end do
    rewind(10) !so that pointer is at the beginnin for reading and importing data later

    print*, " Number of bodies = ", n

    ! Allocate arrays based on the number of stars
    allocate(partics(n))
     
    do i = 1, n
        read(10, *, iostat=ios) partics(i)%m, partics(i)%p%x, partics(i)%p%y, partics(i)%p%z, partics(i)%v
        if (ios /= 0) then
            print *, "Error reading file, stopped at star:", i
            exit
         end if
     end do
     print *, "Values saved to particle type variables."

    ! Close the file
     close(10)
     

  n_t = t_end/dt

! initialize accelerations
  do i = 1,n
     partics(i)%a = vector3d(0,0,0)
  end do
! initialize values
  DO i = 1,n
     DO j = i+1,n
        ! Particles method
        rji_v = vecpp(partics(i)%p, partics(j)%p)
        r3_2 = norm(rji_v)**3 + epsilon **2
        partics(i)%a = partics(i)%a + (partics(j)%m * rji_v / r3_2)
        partics(j)%a = partics(j)%a - (partics(i)%m * rji_v / r3_2)
     END DO
  END DO
  
  t_out = 0.0
  
  ! Open the output file
  open(unit=20, file=output_file, status="replace", action="write", iostat=ios)
     if (ios /= 0) then
        print *, "Error: Could not open file ", output_file
        stop
     end if
     
! start all calculations
  DO i_t = 0, n_t
     t = i_t * dt
     partics(:)%v = partics(:)%v + (partics(:)%a * (dt/2))
     partics(:)%p = partics(:)%p + (partics(:)%v * dt)
     partics(:)%a = vector3d(0,0,0)
    
     DO i = 1,n
        DO j = i+1,n
           rji_v = vecpp(partics(i)%p, partics(j)%p)
           r3_2 = norm(rji_v)**3
           partics(i)%a = partics(i)%a + (partics(j)%m * rji_v / r3_2)
           partics(j)%a = partics(j)%a - (partics(i)%m * rji_v / r3_2)
        END DO
     END DO
     ! Store all values of positions:
     ! t p1x, p1y, ... pNz
     write(fmt, "(A,I0,A,I0,A)") "(F10.2, ", 3*N, "F12.5)"
     write(20, fmt) t, (partics(k)%p, k=1,N)
     partics(:)%v = partics(:)%v + (partics(:)%a * (dt/2))
     t_out = t_out + dt
     IF (t_out >= dt_out) THEN
        print*, "Time = ", t, "", "Iteration = ", i_t
        DO i = 1,n
           !print*, "Particle ", i
           !print*, partics(i)%p
        END DO
        t_out = 0.0
     END IF
     
  END DO

  close(20)
  
END PROGRAM leapfrog
