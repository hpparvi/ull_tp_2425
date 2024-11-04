program e1
    use mpi_f08
    implicit none
    integer :: rank, comsize, ierr

    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, comsize, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

    print *,'Process ', rank, 'of', comsize

    call mpi_finalize(ierr)
 end program e1
