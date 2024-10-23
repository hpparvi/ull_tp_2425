module geometry
  use iso_fortran_env
  implicit none

  type :: vector3d
     real :: x, y, z = real64
  end type vector3d

  type :: point3d 
     real :: x, y, z = real64
  end type point3d

! Operators assignment

  interface operator(+)
     module procedure sumvv, sumvp, sumpv
  end interface
  
  interface operator(-)
     module procedure subvv, subvp, subpv, subpp
  end interface
  
  interface operator(*)
     module procedure mulvr, mulrv, mulvi
  end interface
  
  interface operator(/)
     module procedure divvr
  end interface

contains

! Addition of two vectors
  elemental type(vector3d) function sumvv(vector1, vector2)
    type(vector3d), intent(in) :: vector1, vector2
    sumvv = vector3d(vector1%x + vector2%x, vector1%y + vector2%y, vector1%z + vector2%z)
  end function sumvv
  
  ! Addition of a vector and a point
  elemental type(point3d) function sumvp(vector, point)
    type(vector3d), intent(in) :: vector
    type(point3d), intent(in) :: point
    sumvp = point3d(vector%x + point%x, vector%y + point%y, vector%z + point%z)
  end function sumvp
  
  ! Addition of a point and a vector
  elemental type(point3d) function sumpv(point, vector)
    type(vector3d), intent(in) :: vector
    type(point3d), intent(in) :: point
    sumpv = point3d(vector%x + point%x, vector%y + point%y, vector%z + point%z)
  end function sumpv
  
! Subtraction of two vectors
  elemental type(vector3d) function subvv(vector1, vector2)
    type(vector3d), intent(in) :: vector1, vector2
    subvv = vector3d(vector1%x - vector2%x, vector1%y - vector2%y, vector1%z - vector2%z)
  end function subvv

! Subtraction of a vector with a point
  elemental type(point3d) function subvp(vector, point)
    type(vector3d), intent(in) :: vector
    type(point3d), intent(in) :: point
    subvp = point3d(vector%x - point%x, vector%y - point%y, vector%z - point%z)
  end function subvp
  
  ! Subtraction of a point with a vector
  pure type(point3d) function subpv(point, vector)
    type(vector3d), intent(in) :: vector
    type(point3d), intent(in) :: point
    subpv = point3d(point%x - vector%x, point%y - vector%y, point%z - vector%z)
  end function subpv

! Subtraction of a point with another point
  elemental type(vector3d) function subpp(point1, point2)
    type(point3d), intent(in) :: point1
    type(point3d), intent(in) :: point2
    subpp = vector3d(point1%x - point2%x, point1%y - point2%y, point1%z - point2%z)
  end function subpp

! Multiplication of a vector with a real number function  
  elemental type(vector3d) function mulvr(vector, sreal) 
    real, intent(in) :: sreal
    type(vector3d), intent(in) :: vector
    mulvr = vector3d(vector%x*sreal, vector%y*sreal, vector%z*sreal)
  end function mulvr

! Multiplication of real number with a vector function  
  elemental type(vector3d) function mulrv(sreal, vector) 
    real, intent(in) :: sreal
    type(vector3d), intent(in) :: vector
    mulrv = vector3d(sreal*vector%x, sreal*vector%y, sreal*vector%z)
  end function mulrv

! Multiplication of a vector with a integer number function
  elemental type(vector3d) function mulvi(vector, sint)
    integer, intent(in) :: sint
    type(vector3d), intent(in) :: vector
    mulvi = vector3d(sint*vector%x, sint*vector%y, sint*vector%z)
  end function mulvi
  
! Division of a vector with a real number function
  elemental type(vector3d) function divvr(vector, sreal)
    real, intent(in) :: sreal
    type(vector3d), intent(in) :: vector
    divvr = vector3d(vector%x/sreal, vector%y/sreal, vector%z/sreal)
  end function divvr

! Distance between two points function  
  elemental type(real) function distance(point1, point2)
    type(point3d), intent(in) :: point1, point2
    distance = ((point1%x - point2%x)**2 + (point1%y - point2%y)**2 + (point1%z - point2%z)**2)**0.5
  end function distance
  
! Module of a vector function  
  elemental type(real) function modulev(vector)
    type(vector3d), intent(in) :: vector
    modulev = ((vector%x)**2 + (vector%y)**2 + (vector%z)**2)**0.5
  end function modulev 
  
! Scalar product of two vectors function
  elemental type(real) function dot_prod(vector1, vector2)
    type(vector3d), intent(in) :: vector1, vector2
    dot_prod = (vector1%x*vector2%x + vector1%y*vector2%y + vector1%z*vector2%z)
  end function dot_prod
 
  
! Angle between two vectors function  
  elemental type(real) function angle(vector1, vector2)
    type(vector3d), intent(in) :: vector1, vector2
    angle = acos(dot_prod(vector1,vector2)/(modulev(vector1)*modulev(vector2)))
  end function angle
  
! Vector normalization function
  elemental type(vector3d) function normalize(vector)
    type(vector3d), intent(in) :: vector
    normalize = divvr(vector, modulev(vector))
  end function normalize

! Cross product function  
  elemental type(vector3d) function cross_product(vector1, vector2)
    type(vector3d), intent(in) :: vector1, vector2
    cross_product%x = vector1%y*vector2%z - vector1%z*vector2%y
    cross_product%y = vector1%z*vector2%x - vector1%x*vector2%z
    cross_product%z = vector1%x*vector2%y - vector1%y*vector2%x
  end function cross_product

  
end module geometry




