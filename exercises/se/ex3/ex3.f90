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
  INTEGER(MPI_ADDRESS_KIND), DIMENSION(3) :: offsets_vector
  TYPE(vector3d) :: example_vector
  TYPE(MPI_datatype) :: MPI_vector3d, MPI_point3d, MPI_particle3d
  
  ! For distributing work within MPI
  INTEGER :: npart_per_node, npart_last_node
  INTEGER, DIMENSION(:), ALLOCATABLE :: counts_nodes, displacements
  
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
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: particles

  ! Accelerations; velocities to be received from worker processes (?)
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a

  ! To read the necessary inputs from a file
  INTEGER :: openstatus_input, openstatus_output, readstatus
  CHARACTER(20) :: datafile
  CHARACTER(14) :: resultfile="output_ex2.txt"

  ! To test time of execution
  INTEGER :: count_begin, count_end, count_rate

  ! For tree cells
  TYPE(cell), POINTER :: head, temp_cell
  

  ! Call this first
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD, comsize, ierr)
  CALL mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

  ! Create new structure for MPI
  example_vector = vector3d(1.,1.,1.)
  offsets_vector = [0*sizeof(example_vector%x), sizeof(example_vector%x), 2*sizeof(example_vector%x)]
  CALL MPI_Type_create_struct(3, [1, 1, 1], offsets_vector, &
       [MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION], &
       MPI_vector3d, ierr)
  CALL MPI_Type_commit(MPI_vector3d, ierr)

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Get necessary info from input file (only on 1 thread)  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Begin one thread
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
  ALLOCATE(vels(n))

  ! Assign the masses & initial conditions
  DO i = 1, n
     READ (1, *, IOSTAT=readstatus) particles(i)%m, & 
          particles(i)%p, particles(i)%v
     PRINT '(A, I2)', "This is the mass for particle", i
     PRINT '(F3.1)', particles(i)%m

     PRINT '(A, I2)', "This is the initial position for particle", i
     PRINT '(F7.3)', particles(i)%p
     PRINT '(A, I2)', "This is the initial velocity for particle", i
     PRINT '(F7.3)', particles(i)%v
     PRINT*, "" ! Blank space
  END DO

  CLOSE (UNIT=1) ! Close the file

     
  ! We now know how many particles we have and also how many
  ! MPI processes/nodes. Let us decide how to divide the work:
  
  n_workers = comsize - 1 ! Exclude the main process for simplicity

  IF (n/n_workers > 0) THEN  ! CASE: there are more particles than nodes
     npart_per_node = n/n_workers
     npart_last_node = (n - n_workers*npart_per_node) + npart_per_node
  END IF

     
  ! Start measuring time here
  CALL SYSTEM_CLOCK(count_begin, count_rate)
  ! PRINT '(A, I8)', "The count rate is ", count_rate

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

  ! Set accelerations to 0
  a = vector3d(0., 0., 0.)

  ! THIS CAN BE PARALLELIZED
  ! Get the initial accelerations (recursive)
  CALL Calculate_forces(head, particles, a)

  CALL MPI_Bcast(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(t_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(t_end, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_Bcast(dt_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  IF (rank == 0) THEN
     
     ! Open output file
     OPEN (UNIT=2, FILE=resultfile, STATUS="replace", &
          ACTION="write", POSITION="rewind", &
          IOSTAT=openstatus_output)
     IF (openstatus_output > 0) stop "Cannot open file."


     ! Main loop
     t_out = 0.0
     DO time_counter = 0, total_timesteps
        ! DO t = 0.0, t_end, dt
        particles%v = particles%v + a * (dt/2) 
        particles%p = particles%p + particles%v * dt

        ! Delete the tree because the positions have changed
        CALL Delete_tree(head)

        ! Re-generate the tree
        CALL Calculate_ranges(head, particles)
        head%type = 0
        CALL Nullify_Pointers(head)

        ! This updates the tree
        DO i = 1,n
           CALL Find_Cell(head, temp_cell, particles(i))
           CALL Place_Cell(temp_cell, particles(i), i)
        END DO

        ! Remove the leaves with no particles
        CALL Delete_empty_leaves(head)

        CALL Calculate_masses(head, particles)

        ! Set all accelerations at 0
        a = vector3d(0., 0., 0.)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Parallelization implemented here (and
        ! inside calculate_forces):

        ! Send messages INSIDE CALCULATE FORCES and receive them
        ! If nonviable to do so in both places like with OpenMP,
        ! then maybe redefine calculate_forces to only do with 1 particle
        CALL Calculate_forces(head, particles, a)

        ! What I need to send here is: the velocity and a of each particle
        ! need to split up the array in np pieces, where np is the number of
        ! processes/cores being used (do this manually by taking equal parts
        ! + a remainder)

        ! Scatter the information
        ALLOCATE(counts_nodes(n_workers))
        ALLOCATE(displacements(n_workers))

        ! Set the counts and displacements depending on whether nodes < particles
        IF (n/comsize > 0) THEN
           counts_nodes(:n_workers-1) = npart_per_node
           counts_nodes(n_workers)    = npart_last_node
           DO i=0, n_workers-1
              displacements(i) = i*npart_per_node
           END DO
        ELSE
           PRINT '(A)', "The number of particles is smaller than the number of &
                & nodes. Please, give it another go."
        END IF

        CALL MPI_scatterv(particles%v, counts_nodes, displacements, MPI_vector3d, &
             vel_msg_r, counts_nodes(rank), MPI_vector3d, 0, MPI_COMM_WORLD)
        CALL MPI_scatterv(a, counts_nodes, displacements, MPI_vector3d, &
             a_msg_r, counts_nodes(rank), MPI_vector3d, 0, MPI_COMM_WORLD)

        ! The workers are doing their calculations in between this

        CALL MPI_gather(...)
        
        ! particles%v = particles%v + a * (dt/2)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        t_out = t_out + dt
        ! If t_out is bigger than the increments at which we want output:
        IF (t_out >= dt_out) THEN
           WRITE (2, '(F10.2)', ADVANCE='no') dt*time_counter
           ! For each of the particles
           DO i = 1,n 
              ! Store ALL the positions if it is time to do so
              WRITE (2, '(F15.5, F15.5, F15.5)', ADVANCE='no') particles(i)%p%x, &
                   particles(i)%p%y, particles(i)%p%z
           END DO
           WRITE (2, '(A)') "" ! Just to advance to the next line
           t_out = 0.0
        END IF

     END DO
     ! End of main loop


     ! End measuring time here (is taking into account writing to file)
     CALL SYSTEM_CLOCK(count_end)

     PRINT '(A, F6.3, A)', "The elapsed time is ", &
          REAL(count_end-count_begin) / REAL(count_rate), &
          " seconds"

     CLOSE (UNIT=2)


  ! ===================
  !   WORKER PROCESSES
  ! ===================
  ELSE
     IF (n/comsize > 0) THEN
        

     ELSE
        
        
     END IF
           

  END IF
   


  
  

  ! Finish/close MPI
  CALL mpi_finalize(ierr)



  

END PROGRAM ex3
