program e5
    use mpi_f08
    implicit none
    type(mpi_status) :: status
    real, dimension(:,:), allocatable :: data, res
    real, dimension(:), allocatable :: work
    integer :: i, nwork, nrow, ncol, target, rank, comsize, ierr
    real :: smsg, rmsg

    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, comsize, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    
    if (rank==0) then
       nrow = 20
       ncol = 100
       
       allocate(data(nrow, ncol), res(nrow, ncol))
       
       do i = 1, nrow
          target = 1 + mod(i, comsize-1)
          print *, "Sending row", i, "to", target
          call mpi_send(data(i,:), ncol, MPI_REAL, target, 10, MPI_COMM_WORLD, ierr)
       end do
       
       !do i = 1, comsize-1
       !   call mpi_recv(res(i), 1, MPI_REAL, i, 10, MPI_COMM_WORLD, status, ierr)
       !end do
       !print *, res
    else
       !call mpi_recv(rmsg, 1, MPI_REAL, 0, 10, MPI_COMM_WORLD, status, ierr)
       print *, "Task", rank, "Received ", rmsg, "from the main process"
       
       !rmsg = sqrt(rmsg)
       
       !call mpi_send(rmsg, 1, MPI_REAL, 0, 10, MPI_COMM_WORLD, ierr)
       !print *, "Task", rank, "sent ", rmsg, "to the main process"
    end if
       
    call mpi_finalize(ierr)
 end program e5
