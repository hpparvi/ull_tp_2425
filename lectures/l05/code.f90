module geometry
  implicit none
  type :: vector2d
     real :: x, y
  end type vector2d
end module geometry

program test
  use geometry
  implicit none
  type(vector2d) :: a, b, c
  real :: s = 2.3
  a = vector2d(0.0, 1.0)
  b = vector2d(1.0, 1.0)
  c = vector2d(s*a%x+b%x, s*a%y+b%y)
end program test

module geometry
  implicit none
  type :: vector2d
     real :: x, y
  end type vector2d

contains
  pure type(vector2d) function sumvv(a, b)
    type(vector2d), intent(in) :: a, b
    sumvv = vector2d(a%x + b%x, a%y + b%y)
  end function sumvv
  
  pure type(vector2d) function mulrv(a, b)
    real, intent(in) :: a
    type(vector2d), intent(in) :: b
    mulrv = vector2d(a*b%x, a*b%y)
  end function mulrv
end module geometry

program test
  use geometry
  implicit none
  type(vector2d) :: a, b, c
  real :: s = 2.3
  a = vector2d(0.0, 1.0)
  b = vector2d(1.0, 1.0)
  c = sumvv(mulrv(s, a), b)
end program test

!-------------------------------------------
module geometry
  implicit none
  type :: vector2d
     real :: x, y
  end type vector2d

  interface operator(+)
     module procedure sumvv
  end interface

  interface operator(*)
     module procedure mulrv
  end interface
  
contains
  pure type(vector2d) function sumvv(a, b)
    type(vector2d), intent(in) :: a, b
    sumvv = vector2d(a%x + b%x, a%y + b%y)
  end function sumvv
  
  pure type(vector2d) function mulrv(a, b)
    real, intent(in) :: a
    type(vector2d), intent(in) :: b
    mulrv = vector2d(a*b%x, a*b%y)
  end function mulrv
end module geometry

program test
  use geometry
  implicit none
  type(vector2d) :: a, b, c
  real :: s = 2.3
  a = vector2d(0.0, 1.0)
  b = vector2d(1.0, 1.0)
  c = s*a + b
end program test
!------------------------------------------

interface operator(+)
   module procedure newsum
end interface

interface operator(.op.)
   module procedure newop
end interface

!------------------------------------------

module mymodule
  implicit none

  ! Type definitions

  ! Operator definitions
  
contains

  ! Functions and subroutines
  
end module mymodule

!--------------------------------------------

module geometry
  implicit none

  type :: vector2d
     real :: x, y
  end type vector2d

  interface operator(+)
     module procedure sumvv
  end interface

  interface operator(-)
     module procedure subvv
  end interface

  interface operator(*)
     module procedure mulrv, muliv
  end interface

  interface operator(/)
     module procedure divvr
  end interface

  interface operator(.dot.)
     module procedure dotvv
  end interface

contains
  pure type(vector2d) function sumvv(a, b)
    type(vector2d), intent(in) :: a, b
    sumvv = vector2d(a%x + b%x, a%y + b%y)
  end function sumvv

  pure type(vector2d) function subvv(a, b)
    type(vector2d), intent(in) :: a, b
    subvv = vector2d(a%x - b%x, a%y - b%y)
  end function subvv

  pure type(vector2d) function mulrv(a, b)
    real, intent(in) :: a
    type(vector2d), intent(in) :: b
    mulrv = vector2d(a*b%x, a*b%y)
  end function mulrv

  pure type(vector2d) function muliv(i, v)
    integer, intent(in) :: i
    type(vector2d), intent(in) :: v
    muliv = mulrv(real(i), v)
  end function muliv

  pure type(vector2d) function divvr(v, r)
    real, intent(in) :: r
    type(vector2d), intent(in) :: v
    divvr = vector2d(v%x/r, v%y/r)
  end function divvr

  pure real function dotvv(a, b)
    type(vector2d), intent(in) :: a, b
    dotvv = a%x*b%x + a%y*b%y
  end function dotvv

end module geometry
