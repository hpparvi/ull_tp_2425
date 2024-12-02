
PROGRAM main

  use geometry
  use particles
  use init_conds_mod
  use parameter_reader
  use barnes_hut

  !$ use omp_lib

  IMPLICIT NONE

  INTEGER :: i, k, i_t, ios
  INTEGER :: n, n_t
  INTEGER :: num_args       ! to count the number of given arguments
  character(len=100) :: par_filename, dummy_line     ! to import the stars data and parameters data
  real :: t, t_out !dt, t_end, dt_out
  real :: t_ini, t_fin, t2, t3_i, t3_f, t_dummy
  integer :: it_percent=0, it_percent_previous=0
  character(len=100) :: fmt
  logical :: parallelized = .false.
  
  type(particle), dimension(:), allocatable :: partics

  TYPE (CELL), POINTER :: head   ! for the head cell of the tree and

  !Expected usage:
  ! - command line: executable <custom.par>
  ! The file <custom.par> defines the custom variables (as number of bodies in random distribution
  ! or filename with initial conditions, filename for the output, timestep, total time, gravitational
  ! softening scale and so on...

  ! ------------------------------------------------------------
  ! ------------------------------------------------------------
  ! Block 1: loading all data from the parameters file specified
  ! ------------------------------------------------------------
  ! ------------------------------------------------------------

  call cpu_time(t_ini)

  !$ parallelized = .true. 
  
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

  ! Print info of the parameters specified
    call read_parameters(par_filename)
    print*, "Reading parameters from: ", trim(par_filename)
    print*, "..."
    print*, "Parameters loaded:"
    print '(A10, A10, A10, A10, A10, A16, A10, A10)', "dt", "dt_out","t_end","epsilon", "theta", "create_bodies", "N_bodies", "radius"
    print '(F10.2,F10.2,F10.2,F10.2,F10.2,L16,I10, F10.2)', dt, dt_out, t_end, epsilon, theta, create_bodies, N_bodies, radius
    print*, ""
    print*, "  -  input_file = ", input_file
    print*, "  -  output_file = ", trim(output_file)
    print*, ""
    print'(A,L1)', " - Parallelization: ", parallelized
    if (parallelized .eqv. .true.) then
       print*, "     -  N_threads = ", N_threads
    end if
    
    print*, ""

    ! Set number of threads for parallel mode
    !$ call omp_set_num_threads(N_threads)
    ! Used when going through the tree for calculations for each particle
    
    input_file = trim(input_file)
 
    ! if a new file with the positions of N_bodies is desired,
    ! the create_init_conds routine creates it with the desired name
    if (create_bodies .eqv. .true.) then
       call create_init_conds(input_file, N_bodies, radius)
    end if
    
    ! Open the file with initial data
    open(unit=10, file=input_file, status="old", action="read", form="formatted", iostat=ios)
    if (ios /= 0) then
        print *, "Error: Could not open file ", input_file
        stop
     end if
     
  ! Count the number of stars to check that it matches the desired one
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
    print '(A, I0)', " Number of bodies read: ", n

  ! Allocate arrays based on the number of stars
    allocate(partics(n))
  ! Read and save the input data in the particle elements
    do i = 1, n
        read(10, *, iostat=ios) partics(i)%m, partics(i)%p%x, partics(i)%p%y, partics(i)%p%z, partics(i)%v
        if (ios /= 0) then
            print *, "Error reading file, stopped at star:", i
            exit
         end if
      end do
    print*,""
    print *, "Input data saved to particle type variables."
    close(10)

! Specify the number of time-steps that are going to be required   
  n_t = INT(t_end/dt)

! initialize accelerations
  do i = 1,n
     partics(i)%a = vector3d(0,0,0)
  end do

! Open the output file
  open(unit=20, file=output_file, status="replace", action="write", iostat=ios)
     if (ios /= 0) then
        print *, "Error: Could not open file ", output_file
        stop
     end if
     
call cpu_time(t2)     

print*, ""
print*, "-------------------------------------------------------------"
print*, "Beginning of calculations"
print*, "-------------------------------------------------------------"

  ! ------------------------------------------------------------
  ! ------------------------------------------------------------
  ! Block 2: CALCULATIONS. First step, creation of tree
  ! ------------------------------------------------------------
  ! ------------------------------------------------------------

! Initialize head node
  ALLOCATE(head)
  ! Creation of tree
  CALL Make_Tree(head, partics, n, theta, epsilon)
  
!---------------------------------------
! Main loop
!---------------------------------------
  print*,""
  print*,"..."
  print*,""
  t_out = 0.0
  DO i_t = 0, n_t
     it_percent = i_t*100/n_t
     if (it_percent .ne. it_percent_previous) then
       write(*, '(A, I3, A)', advance="no") " Progress: ", it_percent, "%"
       write(*, '(A)', ADVANCE="NO") CHAR(13)
       it_percent_previous = it_percent
     end if
  
     t = i_t*dt
     partics(:)%v = partics(:)%v + (partics(:)%a * (dt/2))
     partics(:)%p = partics(:)%p + (partics(:)%v * dt)
     ! positions have changed so the tree has to be rebuilt
     CALL Delete_Tree(head)
     CALL Make_Tree(head, partics, n, theta, epsilon)
     
     partics(:)%v = partics(:)%v + (partics(:)%a * (dt/2))
     t_out = t_out + dt

     ! Store all values of positions:
     ! t p1x, p1y, ... pNz
     write(fmt, "(A,I0,A,I0,A)") "(F10.2, ", 3*N, "F12.5)"
     write(20, fmt) t, (partics(k)%p, k=1,N)
          
  END DO

  close(20)

  print*, ""
print*, "-------------------------------------------------------------"
print*, "Calculations finished. ", n_t, " iterations completed." 
print*, "-------------------------------------------------------------"

print*, trim(output_file), " created."

call cpu_time(t_fin)

print '("Total time = ",f0.3," seconds.")',t_fin-t_ini
print '("  - Time spent in input setup = ",f0.3," seconds.")',t2-t_ini
print '("  - Time spent in main loop = ",f0.3," seconds.")',t_fin-t2
  
END PROGRAM main
