module m1
  implicit none
  real, parameter :: pi = 3.14

contains

  real function foo(a)
    real, intent(in) :: a
    foo = sqrt(a)
  end function foo
    
end module m1

