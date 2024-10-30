program ex4
  !$ use omp_lib
  implicit none
  real, dimension(9) :: a
  integer :: nt = 1, tid = 0
  integer :: i

  print *, "Before parallel"
  
  !$omp parallel private(nt, tid, i) shared(a)
  !$ nt = omp_get_num_threads()
  !$ tid = omp_get_thread_num()

  !$omp do
  do i=1,9
     a(i) = i**2 + 3.
     print '(xA,i3,xA,i2,xA,i3)', "Do loop, i =", i, "thread ", tid, " of ", nt
  end do
  !$omp end do
  !$omp end parallel

  print *, "After parallel"
  print '(9f5.1)', a
  
end program ex4
