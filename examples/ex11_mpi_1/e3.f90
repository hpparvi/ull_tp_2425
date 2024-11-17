program e3
    use mpi_f08
    implicit none
    integer :: rank, comsize, ierr
    integer :: i, smsg, rmsg
    type(mpi_status) :: status
    
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, comsize, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

    if (rank==0) then
       do i = 1, comsize-1
          smsg = i**2
          print *, "Sending ", smsg, "to", i
          call mpi_send(smsg, 1, MPI_INTEGER, i, 10, MPI_COMM_WORLD, ierr)
       end do
    else
       call mpi_recv(rmsg, 1, MPI_INTEGER, 0, 10, MPI_COMM_WORLD, status, ierr)
       print *,'Task ', rank, 'received', rmsg
    end if
       
    call mpi_finalize(ierr)
 end program e3
