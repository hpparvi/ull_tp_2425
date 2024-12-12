PROGRAM ex3
  USE, INTRINSIC :: iso_fortran_env
  USE mpi_f08
  USE geometry
  USE particle
  USE tree
  USE create_mpi_types
  IMPLICIT NONE

  TYPE(particle3d), ALLOCATABLE :: parts(:), parts_sub(:)
  TYPE(MPI_Datatype) :: MPI_particle3d, mpi_vector3d, mpi_point3d
  CHARACTER(len=*), PARAMETER :: filename = 'initial_conditions.dat', outname = 'output.dat'
  INTEGER :: i, j, stat, n=0, master=0
  REAL(real64) :: dt, dt_out, t_end, theta = 1., t_out=0., t=0.
  INTEGER :: start_time, end_time, rate
  INTEGER :: rank, size, ierr, istart, iend
  REAL :: elapsed_time
  TYPE(vector3d), ALLOCATABLE :: aa(:), aa_sub(:)
  TYPE(cell), POINTER :: head, temp_cell
  INTEGER, ALLOCATABLE :: sendcounts(:), displs(:)


  ! Initialize mpi
  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, size, ierr)
  CALL mpi_comm_rank(mpi_comm_world, rank, ierr)

  ! these will be used for sending and receiving the messages
  ALLOCATE(sendcounts(size), displs(size))
  
  ! build the derived types (these subroutines are in create_mpi_types)
  CALL build_derived_particle(MPI_particle3d)
  CALL build_derived_vector(MPI_vector3d)
  CALL build_derived_point(MPI_point3d)

  !Start timer
  if (rank.eq.master) then
     CALL system_CLOCK(count_rate = rate)
     CALL system_CLOCK(count = start_time)
  end if

  ! read the initial data with the master node
  IF (rank.EQ.master) THEN
     OPEN (file = filename, action = 'read', status = 'old', unit = 3, iostat = stat) !opens input file
     IF (stat/=0) WRITE (*,*) 'Cannot open file ' , filename

     DO i = 1, 3, 1 ! skips the first 3 lines as those are not particles
        READ(3, *)
     END DO

     DO ! this loop counts the amount of particles
        READ(3, *, iostat = stat)
        IF (stat/=0) EXIT
        n = n + 1
     END DO

     REWIND(3) ! goes back to the beginning of the file

     READ (3, *) dt
     READ (3, *) dt_out
     READ (3, *) t_end

     ALLOCATE(parts(n)) !allocates the particles now that it knows how many there are

     DO i = 1, n !reads initial conditions for all particles
        READ (3, *) parts(i)
     END DO
     CLOSE(3)
  END IF

  ! and send the relevant data to all other nodes
  CALL mpi_bcast(n, 1, MPI_INTEGER, master, mpi_comm_world)

  IF (rank.NE.master) ALLOCATE(parts(n))
  call mpi_bcast(dt_out, 1, MPI_Double_precision, master, mpi_comm_world)
  call mpi_bcast(t_end, 1, MPI_Double_precision, master, mpi_comm_world)
  CALL mpi_bcast(dt, 1, MPI_Double_precision, master, mpi_comm_world)
  CALL mpi_bcast(parts, n, MPI_particle3d, master, mpi_comm_world)

  !! Initializing the tree. we do it with all nodes because pointers
  ALLOCATE(head)
  CALL calculate_ranges(parts, head)
  head%TYPE = 0
  CALL nullify_pointers(head)

  !! Creating the starting tree
  DO i = 1, n
     CALL find_cell(head, temp_cell, parts(i))
     CALL place_cell(temp_cell, parts(i), i)
  END DO

  CALL delete_empty_leaves(head)
  CALL calculate_masses(head)

  ! calculate how much of each array each node will be working on
  sendcounts(:) = FLOOR(REAL(n)/REAL(size))
  sendcounts(size) = FLOOR(REAL(n)/REAL(size)) + MOD(REAL(n), REAL(size))

  
  ALLOCATE(aa(n))
  aa = vector3d(0.,0.,0.)

  ! give each node the start and end id of the particles it works on
  IF (rank.EQ.0) THEN
     istart=1
  ELSE
     istart=(SUM(sendcounts(1:rank)))+1
  END IF
  iend = SUM(sendcounts(1:(rank+1)))

  DO i = 1, size
     IF (i.EQ.1) THEN
        displs(i) = 0
     ELSE
        displs(i) = SUM(sendcounts(1:(i-1)))
     END IF
  END DO


  ! and allocate the subarrays it'll work on
  ALLOCATE(aa_sub(sendcounts(rank+1)))
  ALLOCATE(parts_sub(sendcounts(rank+1)))

  ! each node calculates the accelerations for particles in (istart, iend)
  CALL calculate_forces(head, aa, parts, theta, istart, iend)

  ! and saves it to its subarray
  aa_sub = aa(istart:iend)

  ! wait til everyone is done
  CALL mpi_barrier(mpi_comm_world)

  ! opens the output file ONLY ONCE
  if (rank.eq.master) then
     OPEN (file = outname, action = 'write', status = 'replace', unit =&
          & 4, iostat = stat)
     IF (stat/=0) WRITE (*,*) 'Cannot open file ', outname
  end if

  do while (t .le. t_end) ! main loop
     ! sends everyone their chunk of particles
     CALL mpi_scatterv(parts, sendcounts, displs, MPI_particle3d, &
          & parts_sub, sendcounts(rank+1), MPI_particle3d, master, MPI_comm_world)
     ! leapfrog part 1
     parts_sub%v = parts_sub%v + aa_sub * (dt/2.)
     parts_sub%p = parts_sub%p + parts_sub%v * dt
     ! and then gathers the updated particles and shares it with everyone
     CALL mpi_allgatherv(parts_sub, sendcounts(rank+1), MPI_particle3d,&
          & parts, sendcounts, displs, MPI_particle3d, MPI_COMM_WORLD)

     ! redo the tree (again, pointers)
     call delete_tree(head)

     call calculate_ranges(parts, head)
     head%type = 0
     call nullify_pointers(head)

     do i = 1, n
        call find_cell(head, temp_cell, parts(i))
        call place_cell(temp_cell, parts(i), i)
     end do

     call delete_empty_leaves(head)
     call calculate_masses(head)

     ! and it's back to calculating the forces
     aa = vector3d(0.,0.,0.)

     call calculate_forces(head, aa, parts, theta, istart, iend)
     aa_sub = aa(istart:iend)

     ! re-divide the particles
     CALL mpi_scatterv(parts, sendcounts, displs, MPI_particle3d, &
          & parts_sub, sendcounts(rank+1), MPI_particle3d, master, MPI_comm_world)
     ! leapfrog part 2
     parts_sub%v = parts_sub%v + aa_sub * (dt/2.)
     ! and gather/share again
     CALL mpi_allgatherv(parts_sub, sendcounts(rank+1), MPI_particle3d,&
          & parts, sendcounts, displs, MPI_particle3d, MPI_COMM_WORLD)

     ! update output timer
     t_out = t_out + dt
     if (t_out.ge.dt_out) then
        if (rank.eq.master) then
           ! write ONLY ONCE with the master node
           WRITE(4, fmt='(F11.3)', advance='no') t
           DO i = 1, n
              WRITE(4, fmt='(4ES13.4)', advance='no') parts(i)%p
           END DO
           WRITE(4,*) ''
        end if
        ! wait just in case
        call mpi_barrier(mpi_comm_world)
        t_out = 0.0
     end if
     
     ! update general time
     t=t+dt
  end do

  ! close the file (all file operations are handled by the master node)
  if (rank.eq.master) close(4)

  ! make sure everyone is done
  call mpi_barrier(mpi_comm_world)
  ! and then check the time
  if (rank.eq.master) then
     CALL system_CLOCK(count = end_time)
     elapsed_time = REAL(end_time-start_time)/REAL(rate)
     print *, 'Elapsed time: ', elapsed_time, 'seconds'
  end if

  ! free the data types and finalize mpi
  CALL mpi_type_free(mpi_vector3d)
  CALL mpi_type_free(mpi_point3d)
  CALL mpi_type_free(mpi_particle3d)
  CALL mpi_finalize(ierr)

END PROGRAM ex3
