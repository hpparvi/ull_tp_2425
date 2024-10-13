module geometry

  implicit none

  type :: vector3d
     real :: x, y, z
  end type vector3d

  type :: point3d !What is the difference between point and a vector
  		  ! maybe the vector is like a point but with a direction coordinate assignation
     real :: x, y, z
  end type point3d

! Operators asignation

  interface operator(+)
     module procedure sumvv, sumvp, sumpv
  end interface
  
  interface operator(-)
     module procedure subvv, subvp, subpv
  end interface
  
  interface operator(*)
     module procedure mulvr, mulrv, mulvi
  end interface
  
  interface operator(/)
     module procedure divvr
  end interface

contains

! Addition of two vectors
  pure type(vector3d) function sumvv(vector1, vector2)
    type(vector3d), intent(in) :: vector1, vector2
    sumvv = vector3d(vector1%x + vector2%x, vector1%y + vector2%y, vector1%z + vector2%z)
  end function sumvv
  
  ! Addition of a vector and a point
  pure type(point3d) function sumvp(vector, point)
    type(vector3d), intent(in) :: vector
    type(point3d), intent(in) :: point
    sumvp = point3d(vector%x + point%x, vector%y + point%y, vector%z + point%z)
  end function sumvp
  
  ! Addition of a point and a vector
  pure type(point3d) function sumpv(point, vector)
    type(vector3d), intent(in) :: vector
    type(point3d), intent(in) :: point
    sumpv = point3d(vector%x + point%x, vector%y + point%y, vector%z + point%z)
  end function sumpv
  
! Subtraction of two vectors
  pure type(vector3d) function subvv(vector1, vector2)
    type(vector3d), intent(in) :: vector1, vector2
    subvv = vector3d(vector1%x - vector2%x, vector1%y - vector2%y, vector1%z - vector2%z)
  end function subvv

! Subtraction of a vector with a point
  pure type(point3d) function subvp(vector, point)
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

! Multiplication of a vector with a real number function  
  pure type(vector3d) function mulvr(vector, sreal) 
    real, intent(in) :: sreal
    type(vector3d), intent(in) :: vector
    mulvr = vector3d(sreal*vector%x, sreal*vector%y, sreal*vector%z)
  end function mulvr

! Multiplication of real number with a vector function  
  pure type(vector3d) function mulrv(sreal, vector) 
    real, intent(in) :: sreal
    type(vector3d), intent(in) :: vector
    mulrv = vector3d(sreal*vector%x, sreal*vector%y, sreal*vector%z)
  end function mulrv

! Multiplication of a vector with a integer number function
  pure type(vector3d) function mulvi(vector, sint)
    integer, intent(in) :: sint
    type(vector3d), intent(in) :: vector
    mulvi = vector3d(sint*vector%x, sint*vector%y, sint*vector%z)
  end function mulvi
  
! Division of a vector with a real number function
  pure type(vector3d) function divvr(vector, sreal)
    real, intent(in) :: sreal
    type(vector3d), intent(in) :: vector
    divvr = vector3d(vector%x/sreal, vector%y/sreal, vector%z/sreal)
  end function divvr

! Distance between two points function  
  pure type(real) function distance(point1, point2)
    type(point3d), intent(in) :: point1, point2
    distance = ((point1%x - point2%x)**2 + (point1%y - point2%y)**2 + (point1%z - point2%z)**2)**0.5
  end function distance
  
! Module of a vector function  
  pure type(real) function modulev(vector)
    type(vector3d), intent(in) :: vector
    modulev = ((vector%x)**2 + (vector%y)**2 + (vector%z)**2)**0.5
  end function modulev 
  
! Scalar product of two vectors function
  pure type(real) function mulvv(vector1, vector2)
    type(vector3d), intent(in) :: vector1, vector2
    mulvv = (vector1%x*vector2%x + vector1%y*vector2%y + vector1%z*vector2%z)
  end function mulvv
 
  
! Angle between two vectors function  
  pure type(real) function angle(vector1, vector2)
    type(vector3d), intent(in) :: vector1, vector2
    angle = acos(mulvv(vector1,vector2)/(modulev(vector1)*modulev(vector2)))
  end function angle
  
! Vector normalization function
  pure type(vector3d) function normalize(vector)
    type(vector3d), intent(in) :: vector
    normalize = divvr(vector, modulev(vector))
  end function normalize
  
  pure type(vector3d) function cross_product(vector1, vector2)
    type(vector3d), intent(in) :: vector1, vector2
    cross_product%x = vector1%y*vector2%z - vector1%z*vector2%y
    cross_product%y = vector1%z*vector2%x - vector1%x*vector2%z
    cross_product%z = vector1%x*vector2%y - vector1%y*vector2%x
  end function cross_product

  
end module geometry

program test
  use geometry
  implicit none
  type(vector3d) :: v1, v2, v3
  type(point3d) :: p1, p2
  real :: theta, pi
  real :: dist
  real :: s = 2.
  v1 = vector3d(2.0, 1.0, 2.0)
  v2 = vector3d(0.0, 1.0, 0.0)
  p1 = point3d(0.0, 1.0, 2.0)
  p2 = point3d(2.0, 1.0, 2.0)
  pi= 3.1415927
  v3 = cross_product(v1,v2)
  dist = distance(p1,p2)
  theta = angle(v3,v1)
  print *, v3
  !print *, vector3d(2./3., 1./3., 2./3.)
  print *, theta
  print *, pi/2.
  
end program test


