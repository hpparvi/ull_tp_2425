PROGRAM ex3
  ! USE STATEMENTS
  USE mpi_f08
  USE, INTRINSIC :: iso_fortran_env ! for 64-bit reals
  USE geometry
  USE particle
  USE barnes_hut
  IMPLICIT NONE

  ! For initializing MPI and setting new MPI type (vector3d)
  INTEGER :: rank, comsize, ierr
  TYPE(MPI_datatype) :: MPI_vector3d, MPI_point3d, MPI_particle3d
  
  ! For distributing work within MPI
  INTEGER :: npart_per_node, npart_last_node
  INTEGER, DIMENSION(:), ALLOCATABLE :: counts_nodes, displacements
  INTEGER :: first, last
  !todo Include these in the README
  
  ! (Considering more particles than nodes)
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: vel_msg_per, vel_msg_last
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a_msg_per, a_msg_last
  
  ! Loop indices
  INTEGER :: i,j,k

  ! Number of particles
  INTEGER :: n

  ! Time parameters
  REAL(real64) :: dt, t_end, dt_out, t_out
  INTEGER :: time_counter, total_timesteps
  
  ! Particles
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: particles, particles_node

  ! Accelerations; velocities to be received from worker processes (?)
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a, a_node

  ! To read the necessary inputs from a file
  INTEGER :: openstatus_input, openstatus_output, readstatus
  CHARACTER(20) :: datafile
  CHARACTER(14) :: resultfile="output.txt"

  ! To test time of execution
  INTEGER :: count_begin, count_end, count_rate

  ! For tree cells
  TYPE(cell), POINTER :: head, temp_cell
  

  ! Call this first
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, comsize, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

  ! Create new structures for MPI
  !todo Is it possible to just call the particle, since it internally calls the other 2? 
  CALL create_MPI_vector(MPI_vector3d)
  CALL create_MPI_point(MPI_point3d)
  CALL create_MPI_particle(MPI_particle3d)

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Get necessary info from input file (only on 1 thread)  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Begin one thread
  IF (rank==0) THEN
     PRINT '(A)', "Type in the name of the input file with &
          the initial conditions. Must be in the working directory."
     READ (*,*) datafile

     ! Open the file (already existing)
     OPEN (UNIT=1, FILE=datafile, STATUS="old", &
          ACTION="read", POSITION="rewind", &
          IOSTAT=openstatus_input)
     IF (openstatus_input > 0) stop "Cannot open file. Make sure your &
          input data file is named with < 20 characters."

     ! Read the inputs
     READ (1, *, IOSTAT=readstatus) dt
     READ (1, *, IOSTAT=readstatus) dt_out
     READ (1, *, IOSTAT=readstatus) t_end
     READ (1, '(i5)', IOSTAT=readstatus) n

     ! Make sure it was read correctly and sensibly
     PRINT '(A, F8.4)', "This is the selected time step:", dt
     PRINT '(A, F8.4)', "This is the selected time step for outputs:", dt_out
     PRINT '(A, F8.4)', "This is the selected final time:", t_end
     PRINT '(A, I4)', "This is the selected number of particles:", n
     PRINT*, "" ! Blank space  ! Read the data from the file

     ! Calculate the necessary timesteps to reach end
     total_timesteps = t_end/dt

     ! Allocate arrays once the dimension is known
     ALLOCATE(particles(n))
     ALLOCATE(a(n))

     ! Assign the masses & initial conditions
     DO i = 1, n
        READ (1, *, IOSTAT=readstatus) particles(i)%m, & 
             particles(i)%p, particles(i)%v
        PRINT '(A, I4)', "This is the mass for particle", i
        PRINT '(F3.1)', particles(i)%m

        PRINT '(A, I4)', "This is the initial position for particle", i
        PRINT '(F7.3)', particles(i)%p
        PRINT '(A, I4)', "This is the initial velocity for particle", i
        PRINT '(F7.3)', particles(i)%v
        PRINT*, "" ! Blank space
     END DO

     CLOSE (UNIT=1) ! Close the file
     
  END IF

  CALL MPI_barrier(MPI_COMM_WORLD) ! Avoid premature abort due to worker
  ! processes reaching next step

  ! Send the particle info to all nodes
  !todo If this does not work, bcast n, allocate the particles to all nodes and then bcast particles
  CALL MPI_Bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  IF (rank /= 0) ALLOCATE(particles(n))
  CALL MPI_Bcast(particles, n, MPI_particle3d, 0, MPI_COMM_WORLD, ierr)
   
   
  ! We now know how many particles we have and also how many
  ! MPI processes/nodes. Get the necessary split-up:
  ! (This is done on all processes; very fast)
  IF (n/comsize > 0) THEN  ! CASE: there are more particles than nodes
     npart_per_node = n/comsize ! This is equivalent to FLOOR() when dividing integers (tested)
     npart_last_node = (n - comsize*npart_per_node) + npart_per_node
  ELSE IF (n/comsize < 1 .AND. rank==0) THEN
     PRINT '(A, I3, A, I3, A)', "The number of processes (", comsize, ") exceeds the number &
          of particles (", n, "). Please try again with fewer nodes."
     !todo Make sure this works properly
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

     
  ! Start measuring time here (only 1 node)
  IF (rank==0) THEN
     CALL SYSTEM_CLOCK(count_begin, count_rate)
  END IF

  ! Then, create the tree in all nodes, since pointers cannot be messaged
  ! Initialize the head node
  ALLOCATE(head)

  CALL Calculate_ranges(head, particles) ! space occupied by the head node
  head%type = 0 ! no particle (initialization)
  CALL Nullify_Pointers(head) ! remove all pointers


  ! Create the initial tree
  DO i = 1,n ! For all particles

     ! Locate the cell where this particle should go
     CALL Find_Cell(head, temp_cell, particles(i))

     ! Place the particle inside said cell
     CALL Place_Cell(temp_cell, particles(i), i) ! i : ID of particle
     ! (called n or pos in the subroutine)

  END DO


  ! Remove subcells with no particles inside
  CALL Delete_empty_leaves(head)

  ! Calculate the masses (recursively)
  CALL Calculate_masses(head, particles)



  ! =========================================================================
  ! Now the work sharing begins. First, get all the necessary ints and arrays
  ! =========================================================================

  ! Get the beginning and end of the particles to be assigned to each node
  ! (via the already created integers npart_per_node and npart_last_node)
  !todo Make sure this is correct
  IF (rank == 0) THEN
     first = 1
     last  = npart_per_node
     PRINT*, "Rank:", rank, "First and last:", first, last
  ELSE IF (rank < comsize-1) THEN
     first = rank*npart_per_node + 1 
     last  = first + npart_per_node -1
     PRINT*, "Rank:", rank, "First and last:", first, last
  ELSE
     first = rank*npart_per_node + 1
     last = n
     PRINT*, "Rank:", rank, "First and last:", first, last
  END IF

  ! Arrays for MPI messaging
  ALLOCATE(counts_nodes(comsize))
  ALLOCATE(displacements(comsize))

  ! Set the counts and displacements depending on whether nodes < particles
  IF (n/comsize > 0) THEN
     counts_nodes(:comsize-1) = npart_per_node
     counts_nodes(comsize)    = npart_last_node
     DO i=0, comsize-1
        displacements(i+1) = i*npart_per_node
     END DO
  ELSE
     PRINT '(A)', "The number of particles is smaller than the number of &
          & nodes. Please, give it another go."
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  ! Allocate the smaller arrays within each node
  ALLOCATE(particles_node(counts_nodes(rank+1)))
  ALLOCATE(a_node(counts_nodes(rank+1)))



  ! ================================================================
  ! Now calculate the initial forces, already worksharing
  ! ================================================================
  
  ! Set accelerations to 0
  IF (rank /= 0) ALLOCATE(a(n)) ! all non-master nodes (already done)
  a = vector3d(0., 0., 0.)

  ! Get the initial forces
  
  CALL Calculate_forces(head, particles, a, first, last)
  a_node = a(first:last)

  ! Broadcast the time variables so all processes have them
  CALL MPI_Bcast(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(t_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(time_counter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(total_timesteps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(dt_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)


  ! ==========================================
  ! Time loop and saving the results to a file
  ! ==========================================

  ! Open the output file
  IF (rank == 0) THEN
     ! Open output file
     OPEN (UNIT=2, FILE=resultfile, STATUS="replace", &
          ACTION="write", POSITION="rewind", &
          IOSTAT=openstatus_output)
     IF (openstatus_output > 0) stop "Cannot open file."
  END IF

  ! Main loop
  t_out = 0.0
  DO time_counter = 0, total_timesteps

     ! ==========================================
     ! Sharing workload of leapfrog
     ! ==========================================
     ! * : info received
     CALL MPI_scatterv(particles, counts_nodes, displacements, MPI_particle3d, & ! *
          particles_node, counts_nodes(rank+1), MPI_particle3d, & ! Info received
          0, MPI_COMM_WORLD) ! From where and through which channel
     
     particles_node%v = particles_node%v + a_node * (dt/2) 
     particles_node%p = particles_node%p + particles_node%v * dt

     CALL MPI_allgatherv(particles_node, counts_nodes(rank+1), MPI_particle3d, & ! *
          particles, counts_nodes, displacements, MPI_particle3d, & ! Info received
          MPI_COMM_WORLD) ! Channel
    



     
     ! ========================================
     ! Re-do the whole tree in all processes...
     ! ========================================
     CALL Delete_tree(head)

     CALL Calculate_ranges(head, particles)
     head%type = 0
     CALL Nullify_Pointers(head)

     DO i = 1,n
        CALL Find_Cell(head, temp_cell, particles(i))
        CALL Place_Cell(temp_cell, particles(i), i)
     END DO

     CALL Delete_empty_leaves(head)

     CALL Calculate_masses(head, particles)

     a = vector3d(0., 0., 0.)



     ! ==================================================
     ! Then, calculate the forces by sharing the workload
     ! ==================================================
     CALL Calculate_forces(head, particles, a, first, last)
     a_node = a(first:last) ! The only updated a elements are
     ! those between the particles assigned to the process

     
     ! ========================
     ! Leapfrog time once again
     ! ========================
     CALL MPI_scatterv(particles, counts_nodes, displacements, MPI_particle3d, & ! *
          particles_node, counts_nodes(rank+1), MPI_particle3d, & ! Info received
          0, MPI_COMM_WORLD) ! From where and through which channel
     
     particles_node%v = particles_node%v + a_node * (dt/2) 

     CALL MPI_allgatherv(particles_node, counts_nodes(rank+1), MPI_particle3d, & ! *
          particles, counts_nodes, displacements, MPI_particle3d, & ! Info received
          MPI_COMM_WORLD) ! Channel



     ! ======================================
     ! Save the results to file if it is time
     ! ======================================
     t_out = t_out + dt
     ! If t_out is bigger than the increments at which we want output:
     ! Write only with master node to avoid disaster
     IF (t_out >= dt_out) THEN
        IF (rank==0) THEN
           WRITE (2, '(F10.2)', ADVANCE='no') dt*time_counter
           ! For each of the particles
           DO i = 1,n 
              ! Store ALL the positions if it is time to do so
              WRITE (2, '(F15.5, F15.5, F15.5)', ADVANCE='no') particles(i)%p%x, &
                   particles(i)%p%y, particles(i)%p%z
           END DO
           WRITE (2, '(A)') "" ! Just to advance to the next line
        END IF

        !CALL MPI_barrier(MPI_COMM_WORLD)
        t_out = 0.0
     END IF

  END DO
  ! ======= End of time loop =======

  

  CALL MPI_barrier(MPI_COMM_WORLD) ! Typically not necessary, but here it
  ! prevents the early call of the system clock before all processes are
  ! done with the calculations.
  
  ! End measuring time here (is taking into account writing to file)
  IF (rank==0) THEN
     CALL SYSTEM_CLOCK(count_end)

     PRINT '(A, F6.3, A)', "The elapsed time is ", &
          REAL(count_end-count_begin) / REAL(count_rate), &
          " seconds"

     CLOSE (UNIT=2)
  END IF


       
  
  ! Finish/close MPI
  CALL mpi_finalize(ierr)



CONTAINS

  SUBROUTINE create_MPI_vector(MPI_vector3d)
    TYPE(MPI_datatype), INTENT(out) :: MPI_vector3d
    INTEGER(MPI_ADDRESS_KIND), DIMENSION(3) :: offsets_vector
    TYPE(vector3d) :: example_vector
    INTEGER :: ierr
    
    example_vector = vector3d(1.,1.,1.)
    offsets_vector = [0*sizeof(example_vector%x), sizeof(example_vector%x), 2*sizeof(example_vector%x)]
    CALL MPI_Type_create_struct(3, [1, 1, 1], offsets_vector, &
         [MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION], &
         MPI_vector3d, ierr)
    CALL MPI_Type_commit(MPI_vector3d, ierr)

  END SUBROUTINE create_MPI_vector

  
  SUBROUTINE create_MPI_point(MPI_point3d)
    TYPE(MPI_datatype), INTENT(out) :: MPI_point3d
    INTEGER(MPI_ADDRESS_KIND), DIMENSION(3) :: offsets_point
    TYPE(point3d) :: example_point
    INTEGER :: ierr
    
    example_point = point3d(1.,1.,1.)
    offsets_point = [0*sizeof(example_point%x), sizeof(example_point%x), 2*sizeof(example_point%x)]
    CALL MPI_Type_create_struct(3, [1, 1, 1], offsets_point, &
         [MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION], &
         MPI_point3d, ierr)
    CALL MPI_Type_commit(MPI_point3d, ierr)

  END SUBROUTINE create_MPI_point


  SUBROUTINE create_MPI_particle(MPI_particle3d)
    TYPE(MPI_datatype), INTENT(out) :: MPI_particle3d
    TYPE(MPI_datatype) :: MPI_vector3d, MPI_point3d ! To make particle's (r,v)
    INTEGER(MPI_ADDRESS_KIND), DIMENSION(3) :: offsets_particle
    TYPE(particle3d) :: example_particle
    INTEGER :: ierr

    CALL create_MPI_vector(MPI_vector3d)
    CALL create_MPI_point(MPI_point3d)
    
    example_particle = particle3d(point3d(1., 1., 1.), vector3d(1., 1., 1.), 1)
    offsets_particle = [0*sizeof(example_particle%m), sizeof(example_particle%m), &
         sizeof(example_particle%m) + sizeof(example_particle%p)]
    
    CALL MPI_Type_create_struct(3, [1, 1, 1], offsets_particle, &
         [MPI_DOUBLE_PRECISION, MPI_point3d, MPI_vector3d], &
         MPI_particle3d, ierr)
    CALL MPI_Type_commit(MPI_particle3d, ierr)

  END SUBROUTINE create_MPI_particle
  

END PROGRAM ex3
