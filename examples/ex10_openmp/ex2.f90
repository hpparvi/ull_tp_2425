program ex2
  use omp_lib
  implicit none
  integer :: nt = 1, tid = 0

  print *, "Before parallel"
  
  !$omp parallel private(nt, tid)
  nt = omp_get_num_threads()
  tid = omp_get_thread_num()
  
  print '(XA,i3,XA,i3)', &
       & "Inside parallel, thread ", &
       & tid, " of ", nt
  !$omp end parallel

  print *, "After parallel"
  
end program ex2
