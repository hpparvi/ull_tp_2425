program ex1
  implicit none
  integer :: tid = 0
  
  print *, "Before parallel"
  
  !$omp parallel

  print *, "Inside parallel"

  !$omp end parallel

  print *, "After parallel"
  
end program ex1
