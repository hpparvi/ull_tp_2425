program ex3
  !$ use omp_lib
  implicit none
  integer :: nt = 1, tid = 0

  print *, "Before parallel"
  
  !$omp parallel private(nt, tid)
  !$ nt = omp_get_num_threads()
  !$ tid = omp_get_thread_num()
  
  print '(A17,i3,A5,i3)', "Inside parallel ", tid, " of ", nt
  
  !$omp end parallel

  print *, "After parallel"
  
end program ex3
