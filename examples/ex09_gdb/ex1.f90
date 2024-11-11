program ex1
  use t1
  implicit none
  integer :: i, j
  real :: f
  real, dimension(3) :: s = 0.0
  real, dimension(:), pointer :: x
  real, pointer :: r

  r => x(1)

  call mt(r, 0.3)
  
  do i = 1, 6
     x(i) = x(i) + i * 0.2
     if (i>3) deallocate(x)
  end do

  print *, "sum of 1.2 and 4.3 is", bad_sum(1.2, 4.3)
  print *, "Factorial of 5 is", fac(5)
  print *, "Factorial of 0 is", fac(0)

contains

  integer function fac(n) result(f)
    integer, intent(in) :: n
    integer :: i
    do i = 0, n
       f = f * i
    end do
  end function fac
  
end program ex1
