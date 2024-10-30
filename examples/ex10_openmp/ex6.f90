program ex5
  !$ use omp_lib
  implicit none
  real, dimension(4,3) :: a
  integer :: nt = 1, tid = 0
  integer :: i, j

  print *, "Before parallel"
  
  !$omp parallel private(nt, tid, i, j)
  !$ nt = omp_get_num_threads()
  !$ tid = omp_get_thread_num()

  !$omp do collapse(2)
  do i=1,4
     do j=1,3
        a(i,j) = i+j
        print '(xA,i3,i3,xA,i3,xA,i3)', "Do loop, (i,j) =", i, j, "thread ", tid, " of ", nt
     end do
  end do
  !$omp end do
  !$omp end parallel

  print *, "After parallel"
  do i=1,4
     print '(3f5.1)', a(i,:)
  end do
end program ex5
