module t1
  implicit none

contains
  subroutine mt(v, x)
    real, pointer :: v
    real :: x
    v = sin(x)
  end subroutine mt

  real function bad_sum(a, b)
    real :: a, b
    bad_sum = a - b
  end function bad_sum
end module t1
