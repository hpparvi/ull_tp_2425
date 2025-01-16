program e3
  ! import required libraries/modules  
  use, intrinsic ::  iso_fortran_env
  use geometry
  use particle
  use tree_algorithm
  
  ! included MPI library and MPI new datypes module
  use mpi_f08
  use types_mpi
  IMPLICIT NONE

  INTEGER :: i,j,k,n,rc
  REAL (real64) :: dt, t_end, dt_out, t=0., t_out=0.
  REAL(real64), PARAMETER :: theta = 1
  
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: particles, part_chunk 
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: acc, acc_chunk
  
  CHARACTER(len=*), PARAMETER :: filename = 'init_files/initial_conditions_10b.dat', outname = 'output.dat' ! i.c. input/output files names
  TYPE (CELL), POINTER :: head, temp_cell ! create cell (as pointer)
  
  ! MPI variables for parallelise the program
  INTEGER :: master_rank = 0, rank, pr, ierr
  INTEGER :: n_start, n_end
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: mpistatus
  INTEGER, DIMENSION(:), ALLOCATABLE :: n_counts, n_disp
  
  ! To measure computational time
  INTEGER :: start_time, end_time, count_rate
  REAL :: elapsed_time
  
  ! Let's initialise MPI 
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr) ! process ID number 
  call MPI_COMM_SIZE(MPI_COMM_WORLD, pr, ierr) ! number of processes contained in a communicator
  
  ! Allocate arrays that will save how many particles will carry each chunk
  ! to calculate in every process, and the displacement too
  ALLOCATE(n_counts(pr))
  ALLOCATE(n_disp(pr)) 
  
  ! Call subroutine to create MPI types used in the simulation
  call create_mpi_types()
  
  ! Start timer (only on root process)
  IF (rank == master_rank) THEN
    call system_clock(start_time, count_rate)
  END IF
  
  ! We are going to read the initial condition file in the master node,
  ! and then broadcast the info to slaves nodes
  ! Read init file
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
    ! Read particles info
    DO i = 1, n
      READ (3, *) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
            & particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
    END DO
    
    CLOSE(UNIT=3)
  
  end if 
 
  ! broadcast the initial conditions (number of particles 
  ! and time variables)
  CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(dt,1,MPI_Double_precision,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(dt_out,1,MPI_Double_precision,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(t_end,1,MPI_Double_precision,0,MPI_COMM_WORLD,ierr)
  
  ! All nodes needs particles info to construct the tree
  ! using Barnes-Hut algorithm
  if (rank.ne.master_rank) ALLOCATE(particles(n))
  
  CALL MPI_BCAST(particles,n,mpi_particle3d,0,MPI_COMM_WORLD,ierr) 
  
  !!! Now we proceed to use Barnes-Hut algorithm to implement N body sim 
  
  ! First, initialise head node (done in all processes)
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
  

  ! Calculate how many elements will have each node  
  n_counts(:) = FLOOR(REAL(n)/REAL(pr))
  
  ! Last node will have more particles if n is not
  ! divisible by the number of process (just a few more,
  ! that correspond to the remainder of the division)
  n_counts(pr) = n_counts(pr) + MOD(n,pr) 
  
  ! Allocate particles and acceleration chunk that corresponds
  ! to subsets that will be used for calculation in each node
  ALLOCATE(part_chunk(n_counts(rank+1)))
  ALLOCATE(acc_chunk(n_counts(rank+1))) 
  
  
  ! Calculate where start and end each chunk
  n_start =  rank*FLOOR(REAL(n)/REAL(pr)) + 1
  if ((rank+1).ne.pr) then
    n_end = (rank+1)*FLOOR(REAL(n)/REAL(pr))
  else
    n_end = n
  end if
  
  ! Calculate number of displacements for each chunk
  ! that will be parallel used 
  DO i = 1, pr
     IF (i.EQ.1) THEN
        n_disp(i) = 0
     ELSE
        n_disp(i) = SUM(n_counts(1:(i-1)))
     END IF
  END DO
  
  ! Allocate and calculate initial accelerations
  ALLOCATE(acc(n))
  acc = vector3d(0.0,0.0,0.0)   
  CALL Calculate_forces(head,particles,acc,n_start,n_end)
  ! Save the acceleration calculated in each node
  acc_chunk = acc(n_start:n_end)
  
  ! Make sure that all the process are done
  CALL mpi_barrier(mpi_comm_world)
  
  ! open the output file (master node only)
  if (rank.eq.master_rank) then  
    OPEN (file = outname, action = 'write', status = 'replace', unit = 4, iostat = rc) 
      IF (rc/=0) WRITE (*,*) 'Cannot open file ' , outname
  
    ! the first record of the output file is the initial positions of the particles
    WRITE(4, *) t, particles%p
  end if
  
  !! Main loop 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    DO  WHILE (t <= t_end)
      ! Sends every node their chunk of particles
      CALL mpi_scatterv(particles, n_counts, n_disp, mpi_particle3d, &
           & part_chunk, n_counts(rank+1), mpi_particle3d, master_rank, MPI_comm_world)
      
      part_chunk%v = part_chunk%v + acc_chunk * (dt/2.)
      part_chunk%p = part_chunk%p + part_chunk%v * dt
      
      ! Then gathers particles again and shares it with everyone
      CALL mpi_allgatherv(part_chunk, n_counts(rank+1), mpi_particle3d,&
           & particles, n_counts, n_disp, mpi_particle3d, MPI_COMM_WORLD)
    
      
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
      CALL Calculate_forces(head,particles,acc,n_start,n_end)
      acc_chunk = acc(n_start:n_end)
      
      
      ! Divide the particles set again
      CALL mpi_scatterv(particles, n_counts, n_disp, mpi_particle3d, &
           & part_chunk, n_counts(rank+1), mpi_particle3d, master_rank, MPI_comm_world)
           
      part_chunk%v = part_chunk%v + acc_chunk * (dt/2.)
      
      ! Then gather particles info just like before
      CALL mpi_allgatherv(part_chunk, n_counts(rank+1), mpi_particle3d,&
          & particles, n_counts, n_disp, mpi_particle3d, MPI_COMM_WORLD)
       
      t_out = t_out + dt
      IF (t_out >= dt_out) THEN
        IF (rank.eq.master_rank) THEN
          WRITE(4, *) t, particles%p ! time and positions in one row (one particle position after another)
          ! make sure every node is in the same step
        END IF
        t_out = 0.0
      END IF
      ! Make sure all nodes are in the same step
      CALL mpi_barrier(mpi_comm_world)
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
