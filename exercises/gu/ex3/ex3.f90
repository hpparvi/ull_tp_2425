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


  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, size, ierr)
  CALL mpi_comm_rank(mpi_comm_world, rank, ierr)

  ALLOCATE(sendcounts(size), displs(size))

  CALL build_derived_particle(MPI_particle3d)
  CALL build_derived_vector(MPI_vector3d)
  CALL build_derived_point(MPI_point3d)

  !Start timer
  CALL system_CLOCK(count_rate = rate)
  CALL system_CLOCK(count = start_time)

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
     ALLOCATE(aa(n))

     DO i = 1, n !reads initial conditions for all particles
        READ (3, *) parts(i)
     END DO
     CLOSE(3)
  END IF

  CALL mpi_bcast(n, 1, MPI_INTEGER, master, mpi_comm_world)

  IF (rank.NE.master) ALLOCATE(parts(n))
  call mpi_bcast(dt_out, 1, MPI_Double_precision, master, mpi_comm_world)
  call mpi_bcast(t_end, 1, MPI_Double_precision, master, mpi_comm_world)
  CALL mpi_bcast(dt, 1, MPI_Double_precision, master, mpi_comm_world)
  CALL mpi_bcast(parts, n, MPI_particle3d, master, mpi_comm_world)

  !! Initializing the head node
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

  sendcounts(:) = FLOOR(REAL(n)/REAL(size))
  sendcounts(size) = FLOOR(REAL(n)/REAL(size)) + MOD(REAL(n), REAL(size))

  IF (rank.NE.master) ALLOCATE(aa(n))
  aa = vector3d(0.,0.,0.)

  IF (rank.EQ.0) THEN
     istart=1
  ELSE
     istart=(SUM(sendcounts(1:rank)))+1
  END IF
  iend = SUM(sendcounts(1:(rank+1)))

  ALLOCATE(aa_sub(sendcounts(rank+1)))
  ALLOCATE(parts_sub(sendcounts(rank+1)))

  CALL calculate_forces(head, aa, parts, theta, istart, iend)

  aa_sub = aa(istart:iend)

  CALL mpi_barrier(mpi_comm_world)

  DO i = 1, size
     IF (i.EQ.1) THEN
        displs(i) = 0
     ELSE
        displs(i) = SUM(sendcounts(1:(i-1)))
     END IF
  END DO


!!$  CALL mpi_gatherv(aa_sub, sendcounts(rank+1), MPI_vector3d,&
!!$       &aa, sendcounts, displs, MPI_vector3d, master,&
!!$       &MPI_comm_world)
  !CALL mpi_bcast(aa, n, MPI_vector3d, master, mpi_comm_world)


  if (rank.eq.master) then
     OPEN (file = outname, action = 'write', status = 'replace', unit =&
          & 4, iostat = stat)
     IF (stat/=0) WRITE (*,*) 'Cannot open file ', outname
  end if

  t_out = 0.0
  t=0.

  do while (t .le. t_end)
!!$     CALL mpi_scatterv(aa, sendcounts, displs, MPI_vector3d, &
!!$          &aa_sub, sendcounts(rank+1), MPI_vector3d, master, mpi_Comm_World)
     CALL mpi_scatterv(parts, sendcounts, displs, MPI_particle3d, &
          & parts_sub, sendcounts(rank+1), MPI_particle3d, master, MPI_comm_world)
     parts_sub%v = parts_sub%v + aa_sub * (dt/2.)
     parts_sub%p = parts_sub%p + parts_sub%v * dt
     CALL mpi_gatherv(parts_sub, sendcounts(rank+1), MPI_particle3d,&
          & parts, sendcounts, displs, MPI_particle3d, master, MPI_COMM_WORLD)
     CALL mpi_bcast(parts, n, MPI_particle3d, master, mpi_comm_world)
     
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

     aa = vector3d(0.,0.,0.)

     call calculate_forces(head, aa, parts, theta, istart, iend)
     aa_sub = aa(istart:iend)

     CALL mpi_scatterv(parts, sendcounts, displs, MPI_particle3d, &
          & parts_sub, sendcounts(rank+1), MPI_particle3d, master, MPI_comm_world)
     parts_sub%v = parts_sub%v + aa_sub * (dt/2.)
     CALL mpi_gatherv(parts_sub, sendcounts(rank+1), MPI_particle3d,&
          & parts, sendcounts, displs, MPI_particle3d, master, MPI_COMM_WORLD)
     call mpi_bcast(parts, n, MPI_particle3d, master, mpi_comm_world)
     t_out = t_out + dt
     if (t_out.ge.dt_out) then
        if (rank.eq.master) then
           WRITE(4, fmt='(F11.3)', advance='no') t
           DO i = 1, n
              WRITE(4, fmt='(4ES13.4)', advance='no') parts(i)%p
           END DO
           WRITE(4,*) ''
        end if
        call mpi_barrier(mpi_comm_world)
        t_out = 0.0
     end if
     t=t+dt
  end do

  if (rank.eq.master) close(4)

  call mpi_barrier(mpi_comm_world)
  if (rank.eq.master) then
     CALL system_CLOCK(count = end_time)
     elapsed_time = REAL(end_time-start_time)/REAL(rate)
     print *, 'Elapsed time: ', elapsed_time, 'seconds'
  end if

  CALL mpi_type_free(mpi_vector3d)
  CALL mpi_type_free(mpi_point3d)
  CALL mpi_type_free(mpi_particle3d)
  CALL mpi_finalize(ierr)

END PROGRAM ex3
