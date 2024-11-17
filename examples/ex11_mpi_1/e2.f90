program e2
    use mpi_f08
    implicit none
    integer :: rank, comsize, ierr

    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, comsize, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

    if (rank==0) then
       print *, 'Main process'
    else
       print *,'Process ', rank, 'of', comsize
    end if
       
    call mpi_finalize(ierr)
 end program e2
