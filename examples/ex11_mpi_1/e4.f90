program e4
    use mpi_f08
    implicit none
    integer :: i, rank, comsize, ierr
    type(mpi_status) :: status
    real, dimension(:), allocatable :: res
    real :: smsg, rmsg

    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, comsize, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! Main process
    ! ------------
    if (rank==0) then
       allocate(res(comsize-1))
       
       do i = 1, comsize-1
          smsg = real(i)
          print *, "Sending ", smsg, "to", i
          call mpi_send(smsg, 1, MPI_REAL, i, 10, MPI_COMM_WORLD, ierr)
       end do
       
       do i = 1, comsize-1
          call mpi_recv(res(i), 1, MPI_REAL, i, 10, MPI_COMM_WORLD, status, ierr)
       end do
       print *, res

    ! Worker process
    ! --------------
    else
       call mpi_recv(rmsg, 1, MPI_REAL, 0, 10, MPI_COMM_WORLD, status, ierr)
       print *, "Task", rank, "Received ", rmsg, "from the main process"
       
       rmsg = sqrt(rmsg)
       
       call mpi_send(rmsg, 1, MPI_REAL, 0, 10, MPI_COMM_WORLD, ierr)
       print *, "Task", rank, "sent ", rmsg, "to the main process"
    end if
       
    call mpi_finalize(ierr)
 end program e4
