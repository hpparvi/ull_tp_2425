program e3
  ! import required libraries/modules  
  use, intrinsic ::  iso_fortran_env
  use geometry
  use particle
  use tree_algorithm
  
  ! included MPI library and module
  use mpi_f08
  use types_mpi
  IMPLICIT NONE

  INTEGER :: i,j,k,n,rc
  REAL (real64) :: dt, t_end, t, dt_out, t_out
  REAL(real64), PARAMETER :: theta = 1
  
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: particles 
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: acc !acceleration
  CHARACTER(len=*), PARAMETER :: filename = 'init_files/example.dat', outname = 'output.dat' ! i.c. input/output files names
  TYPE (CELL), POINTER :: head, temp_cell ! create cell (as pointer)
  
  ! MPI variables for parallelise the program
  INTEGER :: master_rank = 0, rank, pr, ierr
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: mpistatus
  !INTEGER :: my_n, 
  
  ! To measure computational time
  INTEGER :: start_time, end_time, count_rate
  REAL :: elapsed_time
  
  ! Let's initialise MPI 
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr) ! process ID number 
  call MPI_COMM_SIZE(MPI_COMM_WORLD, pr, ierr) ! number of processes contained in a communicator
  
  ! Call subroutine to create MPI types used in the simulation
  call create_mpi_types()
  
  ! Start timer (only on root process)
  IF (rank == master_rank) THEN
    call system_clock(start_time, count_rate)
  END IF
  
  ! Only master reads init file
  if (rank == master_rank) then  
    ! open the input file
    OPEN (file = filename, action = 'read', status = 'old', unit = 3, iostat = rc)
    IF (rc/=0) WRITE (*,*) 'Cannot open file ' , filename 
   
    ! save the initial conditions
    READ (3, *) dt
    READ (3, *) dt_out
    READ (3, *) t_end
    READ (3, *) n
    
    ! If there are more processes than particles, abort execution
    if (n.lt.pr) then
      print *, 'Parallelization error: More processes than particles.'
      call mpi_abort(mpi_comm_world, 13)
      call mpi_finalize(ierr)
      stop
    end if
   
    ALLOCATE(particles(n))
    
    DO i = 1, n
      READ (3, *) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
            & particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
    END DO
    
    CLOSE(UNIT=3)
  
  end if 
 
  ! broadcast the initial conditions (number of particles 
  ! and time variables)
  CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(dt,1,MPI_REAL,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(dt_out,1,MPI_REAL,0,MPI_COMM_WORLD,error)
  CALL MPI_BCAST(t_end,1,MPI_REAL,0,MPI_COMM_WORLD,error)
  
  if (rank.nq.master_rank) ALLOCATE(particles(n))
  
  CALL MPI_BCAST(particles,n,mpi_particle3d,0,MPI_COMM_WORLD,error) 
  
  !!--------------------------------------------------------
  !! Calculate block to each node 
  !! this would be can be divided by number of processes()
  
  my_n = n/pr ! number of particles to each node 
  my_start = (my_n * my_rank) + 1 
  my_end = my_start + my_n - 1
  
  !!-----------------------------------------
  
  
  !! Initialise head node 
  ALLOCATE(head)
  CALL Calculate_ranges(head, particles) 
  head%type = 0 ! no particle
  CALL Nullify_Pointers(head) ! null all the pointers first just in case
  

  ! Create initial tree (must be calculated in all the nodes)
  DO i = 1,n
    CALL Find_Cell(head,temp_cell,particles(i)) 
    CALL Place_Cell(temp_cell,particles(i),i)

  END DO
  
  CALL Borrar_empty_leaves(head)
  CALL Calculate_masses(head, particles)

  ! Allocate and calculate initial accelerations
  allocate(acc)
  acc = vector3d(0.0,0.0,0.0)
    
    
  CALL Calculate_forces(head,particles,n,acc)
    
  ! open the output file (master node only)
  if (rank.eq.master_rank) then  
    OPEN (file = outname, action = 'write', status = 'replace', unit = 4, iostat = rc) 
      IF (rc/=0) WRITE (*,*) 'Cannot open file ' , outname
  
    ! the first record of the output file is the initial positions of the particles
    WRITE(4, *) t, particles%p
  end if
  
  !! Main loop 
  !!!!!!!!!!!!!!!!!!
    t_out = 0.0
    DO  WHILE (t <= t_end)

      particles%v = particles%v + acc * (dt/2.)
      particles%p = particles%p + particles%v * dt
      
      CALL Borrar_tree(head) ! remove previous tree
      CALL Calculate_ranges(head, particles) ! calculate head range again
      head%type = 0 
      CALL Nullify_Pointers(head) 
      
      DO i = 1,n
        CALL Find_Cell(head,temp_cell,particles(i))
        CALL Place_Cell(temp_cell,particles(i),i)
      END DO
      
      CALL Borrar_empty_leaves(head)
      CALL Calculate_masses(head, particles)
      
      acc = vector3d(0.0,0.0,0.0)
 
      CALL Calculate_forces(head,particles,n,pr,acc)
      
      particles%v = particles%v + acc * (dt/2.)
      
      t_out = t_out + dt
      IF (t_out >= dt_out) THEN
        IF (rank.eq.master_rank) THEN
          WRITE(4, *) t, particles%p ! time and positions in one row (one particle position after another)
          t_out = 0.0
        END IF
      END IF
      
      t = t + dt
      
    END DO 
  
  ! Close output file (but again only in the master node)
  IF (rank.eq.master_rank) CLOSE(UNIT=4)
  
  ! End timer (only on root process)
  if (rank.eq.master_rank) then
  
    call system_clock(end_time)
    
    ! Calculate elapsed time
    elapsed_time = REAL(end_time-start_time) / REAL(count_rate)
    
    ! Output the elapsed time
    PRINT *, "The elapsed time is ", elapsed_time, " seconds"
  end if
  
  ! Cleanup routine for new mpi types previously created
  CALL free_mpi_types()

  CALL MPI_FINALIZE(ierr)

end program e3
