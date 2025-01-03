
PROGRAM main

  use geometry
  use particles
  use init_conds_mod
  use parameter_reader
  use barnes_hut

  !$ use omp_lib
  use mpi_f08
  use mpi_particles

  IMPLICIT NONE

  INTEGER :: i, k, i_t, ios
  INTEGER :: n, n_t
  INTEGER :: num_args       ! to count the number of given arguments
  character(len=100) :: par_filename, dummy_line     ! to import the stars data and parameters data
  real :: t, t_out !dt, t_end, dt_out
  !integer :: clock_ini, clock_fin, clock_rate, clock_dummy, clock_dummy_2
  real :: t_total, t_ini, t_fin, t2, t_tree=0, t_del_tree=0, t_parts=0, t_write=0, t_bcast = 0, t_dummy = 0, t_dummy_2=0, t_gather = 0
  real, dimension(7) :: t_calcs = [0,0,0,0,0,0,0]
  integer :: it_percent=0, it_percent_previous=0
  character(len=100) :: fmt
  
  type(particle), dimension(:), allocatable :: partics

  TYPE (CELL), POINTER :: head   ! for the head cell of the tree and

  ! for MPI usage
  integer :: ierr, rank, size
  type(MPI_Datatype) :: MPI_particle, MPI_vector3d, MPI_point3d
  integer :: master_rank = 0, i_send_min, i_send_max
  integer, allocatable :: n_send(:), displacements(:)




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

  ! -----------------MPI initialization ------------------------------------------
  
  CALL MPI_INIT(ierr)
  ! Obtain the rank and size
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

  ! The whole reading and settig up parameters is done only by master process
  if (rank .eq. master_rank) then
     t_ini = MPI_Wtime()
     !call system_clock(clock_ini, clock_rate)
     ! call cpu_time(t_ini)

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
     print*, "  -  input_file = ", trim(input_file)
     print*, "  -  output_file = ", trim(output_file)
     print*, ""
     print'(A,L1)', " - Parallelization: ", .true.
     print*, "     -  N_threads = ", size
     print*, ""

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

     t2 = MPI_Wtime()

     print*, ' Broadcasting data to all processes...'

  end if

  ! Create mpi_particle datatype
  call create_mpi_datatypes(MPI_point3d, MPI_vector3d, MPI_particle)

  ! Sharing of data (number of particles and timesteps, total time
  ! calculation parameters and particle information)
  
  CALL MPI_Bcast(n,   1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(n_t, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(t,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(epsilon, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(theta,  1,  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  if (rank.ne.master_rank) ALLOCATE(partics(n))
  CALL MPI_Bcast(partics, n, MPI_particle, 0, MPI_COMM_WORLD, ierr)

  ! number of particles to send to each process:
  allocate(n_send(size))
  n_send(:) = floor(real(n) / real(size))
  ! remaining bodies sent to last node
  n_send(size) = n - (size-1) * floor(real(n)/real(size))
  ! indexes of minimum and maximum value sent
  if (rank.eq.master_rank) then
     i_send_min = 1
  else
     i_send_min = SUM(n_send(:rank)) +1
  end if
  i_send_max = SUM(n_send(:rank+1))
  ! displacements
  allocate(displacements(size))
  displacements(1) = 0
  do i = 2, size
     displacements(i) = sum(n_send(:i-1))
  end do     

  if (rank.eq.master_rank) then
     t_bcast = MPI_Wtime()
  end if
  
! -------------------MPI initialization end ------------------------------------------


 
  ! ------------------------------------------------------------
  ! ------------------------------------------------------------
  ! Block 2: CALCULATIONS. First step, creation of tree
  ! 
  ! Creation of tree will be done on all nodes, and force
  ! calculations will be the only thing parallelized.
  ! Master node will divide particle information and write output.
  ! ------------------------------------------------------------
  ! ------------------------------------------------------------
  if (rank.eq.master_rank) then
     print*, ""
     print*, "-------------------------------------------------------------"
     print*, "Beginning of calculations"
     print*, "-------------------------------------------------------------"
  end if

! Initialize head node
  ALLOCATE(head)
  ! Creation of tree
  CALL Make_Tree(head, partics, n, theta, epsilon, t_calcs, n_send, master_rank, rank)
  if (rank.eq.master_rank) then
     t_tree = MPI_Wtime() - t_bcast
  end if
  ! After tree calculations each rank has updated its cells accelerations. Now it has to 
  !call system_clock(clock_dummy_2)
  !call system_clock(clock_dummy)
  
  !t_tree = real(clock_dummy_2-clock_dummy) / real(clock_rate)
  
!---------------------------------------
! Main loop
  !---------------------------------------
  if (rank.eq.master_rank) then
     print*,""
     print*,"..."
     print*,""
     t_out = 0.0
  end if
  DO i_t = 0, n_t
     if (rank.eq.master_rank) then
        it_percent = i_t*100/n_t
        if (it_percent .ne. it_percent_previous) then
           write(*, '(A, I3, A)', advance="no") " Progress: ", it_percent, "%"
           write(*, '(A)', ADVANCE="NO") CHAR(13)
           it_percent_previous = it_percent
        end if
     end if

     t = i_t*dt
     
     t_dummy = MPI_Wtime()
     partics(i_send_min:i_send_max)%v = partics(i_send_min:i_send_max)%v + (partics(i_send_min:i_send_max)%a * (dt/2))
     partics(i_send_min:i_send_max)%p = partics(i_send_min:i_send_max)%p + (partics(i_send_min:i_send_max)%v * dt)
     t_dummy_2 = MPI_Wtime()

     t_parts = t_parts + t_dummy_2-t_dummy

     ! positions have changed so the tree has to be rebuilt. First we update all particle information to all processes
     CALL MPI_Allgatherv(MPI_IN_PLACE, n_send(rank + 1), MPI_particle, &
          & partics, n_send, displacements, MPI_particle, MPI_COMM_WORLD, ierr)
     t_dummy = MPI_Wtime()
     t_gather = t_gather + t_dummy - t_dummy_2
     t_dummy_2 = t_dummy
     
     CALL Delete_Tree(head)
     !call system_clock(clock_dummy)
     t_dummy = MPI_Wtime()
     t_del_tree = t_del_tree + t_dummy - t_dummy_2
     
     CALL Make_Tree(head, partics, n, theta, epsilon, t_calcs, n_send, master_rank, rank)
     !call system_clock(clock_dummy_2)
     t_dummy_2 = MPI_Wtime()
     t_tree = t_tree + t_dummy_2 - t_dummy
     
     partics(:)%v = partics(:)%v + (partics(:)%a * (dt/2))
     !call system_clock(clock_dummy)
     t_dummy = MPI_Wtime()
     t_parts = t_parts + t_dummy - t_dummy_2
     
     t_out = t_out + dt

     if (t_out >= dt_out) then
        if (rank.eq.master_rank) then
           ! Store all values of positions:
           ! t p1x, p1y, ... pNz
           write(fmt, "(A,I0,A,I0,A)") "(F10.2, ", 3*N, "F12.5)"
           write(20, fmt) t, (partics(k)%p, k=1,N)
           t_dummy_2 = MPI_Wtime()
           t_write = t_write + t_dummy_2 - t_dummy
           t_out = 0.0
        end if
     end if
  
  END DO
  
  if (rank.eq.master_rank) then

     close(20)

     print*, ""
     print*, "-------------------------------------------------------------"
     print*, "Calculations finished. ", n_t, " iterations completed." 
     print*, "-------------------------------------------------------------"

     print*, trim(output_file), " created."

     t_fin = MPI_Wtime()
     !call system_clock(clock_fin)
     t_total =t_fin - t_ini

     print '("Total wall time = ",f0.3," seconds.")', t_fin-t_ini
     print '("  - Time spent in input setup = ",f0.3," seconds.")',t2-t_ini
     print '("     - Time broadcasting initial values = ",f0.3," seconds.")',t_bcast-t2
     print '("  - Time spent in main loop = ",f0.3," seconds.")',t_fin-t2
     print '("     - Makig trees = ",f0.3," seconds.")', t_tree
     print '("        - Calculating ranges = ",f0.3," seconds.")', t_calcs(1)
     print '("        - Nullifying pointers = ",f0.3," seconds.")', t_calcs(2)
     print '("        - Finding and placing cells = ",f0.3," seconds.")', t_calcs(3)
     print '("        - Deleting empty leaves = ",f0.3," seconds.")', t_calcs(4)
     print '("        - Calculating masses = ",f0.3," seconds.")', t_calcs(5)
     print '("        - Initializing accelerations = ",f0.3," seconds.")', t_calcs(6)
     print '("        - Calculating forces = ",f0.3," seconds.")', t_calcs(7)
     print '("        - Gathering process information = ",f0.3," seconds.")', t_gather
     print '("     - Deleting trees = ",f0.3," seconds.")', t_del_tree
     print '("     - Updating particles = ",f0.3," seconds.")', t_parts
     print '("     - Writing results = ",f0.3," seconds.")', t_write

     ! Print with comma-separated values
     print*, "All time values for copying into python and plotting:"
     WRITE(*, "(F6.2,',')", advance='NO') t_total
     WRITE(*, "(F6.2,',')", advance='NO') t_tree
     WRITE(*, "(F6.2,',')", advance='NO') t_calcs(1)
     WRITE(*, "(F6.2,',')", advance='NO') t_calcs(2)
     WRITE(*, "(F6.2,',')", advance='NO') t_calcs(3)
     WRITE(*, "(F6.2,',')", advance='NO') t_calcs(4)
     WRITE(*, "(F6.2,',')", advance='NO') t_calcs(5)
     WRITE(*, "(F6.2,',')", advance='NO') t_calcs(6)
     WRITE(*, "(F6.2,',')", advance='NO') t_calcs(7)
     WRITE(*, "(F6.2,',')", advance='NO') t_del_tree
     WRITE(*, "(F6.2,',')", advance='NO') t_parts
     WRITE(*, "(F6.2,',')", advance='NO') t_write
     WRITE(*, "(F6.2)", advance='YES') t_gather
  end if

  CALL MPI_Finalize()


  
END PROGRAM main
