program ex7
  !$ use omp_lib
  implicit none
  real, dimension(5) :: a, b
  real :: res = 0
  integer :: i

  a = [0.1, 3.4, 7.4, 1.2, 8.3]
  b = [3.8, 0.1, 1.6, 0.1, 0.3]

  !$omp parallel private(i) shared(a, b, res)
  !$omp do reduction(+ : res)
  do i=1,5
     res = res + a(i) * b(i)
  end do
  !$omp end do
  !$omp end parallel

  print *, res
end program ex7
