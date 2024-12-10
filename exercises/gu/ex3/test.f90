program test
  use, intrinsic :: iso_fortran_env
  use mpi_f08
  use geometry
  use particle
  use tree
  use create_mpi_types
  implicit none

  type(particle3d), allocatable :: parts(:)
  type(MPI_Datatype) :: MPI_particle3d, mpi_vector3d, mpi_point3d
  CHARACTER(len=*), PARAMETER :: filename = 'initial_conditions.dat'
  INTEGER :: i, j, stat, n=0, master=0
  REAL(real64) :: dt, dt_out, t_end, theta = 1., t_out=0., t=0.
  INTEGER :: start_time, end_time, rate
  integer :: rank, size, ierr, istart, iend
  REAL :: elapsed_time
  TYPE(vector3d), ALLOCATABLE :: aa(:), aa_sub(:)
  TYPE(cell), POINTER :: head, temp_cell
  integer, allocatable :: sendcounts(:), displs(:)
  

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, size, ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)

  allocate(sendcounts(size), displs(size))

  call build_derived_particle(MPI_particle3d)
  call build_derived_vector(MPI_vector3d)
  call build_derived_point(MPI_point3d)

  !Start timer
  call system_clock(count_rate = rate)
  call system_clock(count = start_time)

  if (rank.eq.master) then
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
  end if

  call mpi_bcast(n, 1, MPI_INTEGER, master, mpi_comm_world)

  if (rank.ne.master) allocate(parts(n))
  call mpi_bcast(dt, 1, MPI_Double_precision, master, mpi_comm_world)
  call mpi_bcast(parts, n, MPI_particle3d, master, mpi_comm_world)

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

  sendcounts(:) = floor(real(n)/real(size))
  sendcounts(size) = floor(real(n)/real(size)) + mod(real(n), real(size))

  if (rank.ne.master) allocate(aa(n))
  aa = vector3d(0.,0.,0.)
  
  if (rank.eq.master) then
     aa = vector3d(0.,0.,0.)
  end if

  if (rank.eq.0) then
     istart=1
  else
     istart=(sum(sendcounts(1:rank)))+1
  end if
  iend = sum(sendcounts(1:(rank+1)))
  print *, 'rank ', rank, 'istart ', istart, 'iend ', iend

  allocate(aa_sub(istart:iend))

  do i=1,size
     displs(i) = (i-1)
  end do


  CALL mpi_scatterv(aa, sendcounts, displs, MPI_vector3d, &
       &aa_sub, sendcounts(rank+1), MPI_vector3d, master, mpi_comm_world)
  print *, 'got this far'

  print *, aa_sub(istart:iend), aa_sub(1)
  !call calculate_forces(head, aa_sub, parts, theta, istart, iend)
  
  call system_clock(count = end_time)
  elapsed_time = real(end_time-start_time)/real(rate)
  print *, rank, 'Elapsed time: ', elapsed_time, 'seconds'
  
  call mpi_type_free(mpi_vector3d)
  call mpi_type_free(mpi_point3d)
  call mpi_type_free(mpi_particle3d)
  call mpi_finalize(ierr)
  
end program test
