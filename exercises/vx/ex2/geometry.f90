module geometry        !This module defines 3D vector and point operations for vector3d and point3d types.
  use iso_fortran_env  !This module ensures all variables are defined as 64-bit.
  implicit none

  
  !Definition of types.
  type :: vector3d
     real(real64) :: x, y, z
  end type vector3d

  type :: point3d
     real(real64) :: x, y, z
  end type point3d

  
  !Definition of custom operators.
  interface operator (+)
     module procedure sumvp, sumpv, sumpp, sumvv
  end interface operator (+)

  interface operator (-)
     module procedure subvp, subpv, subpp, subvv
  end interface operator (-)

  interface operator (*)
     module procedure mulrv, mulvr, mulrp, mulpr
  end interface operator (*)

  interface operator (/)
     module procedure divvr, divpr
  end interface operator (/)
  
  

contains
  
  !Definition of sums between points and vectors.
  pure type(point3d) function sumvp(v, p)            
    type(vector3d), intent(in) :: v
    type(point3d), intent(in) :: p
    sumvp = point3d(v%x + p%x, v%y + p%y, v%z + p%z)
  end function sumvp

  pure type(point3d) function sumpv(p, v)             
    type(point3d), intent(in) :: p
    type(vector3d), intent(in) :: v
    sumpv = point3d(p%x + v%x, p%y + v%y, p%z + v%z)
  end function sumpv

  pure type(vector3d) function sumvv(a, b)            
    type(vector3d), intent(in) :: a, b
    sumvv = vector3d(a%x + b%x, a%y + b%y, a%z + b%z)
  end function sumvv

  pure type(vector3d) function sumpp(a, b)
    type(point3d), intent(in) :: a, b
    sumpp = vector3d(a%x + b%x, a%y + b%y, a%z + b%z)
  end function sumpp

  
  !Definition of substractions between points and vectors.
  pure type(vector3d) function subvv(a, b)            
    type(vector3d), intent(in) :: a, b
    subvv = vector3d(a%x - b%x, a%y - b%y, a%z - b%z)
  end function subvv

  pure type(vector3d) function subvp(v, p)
    type(vector3d), intent(in) :: v
    type(point3d), intent(in) :: p
    subvp = vector3d(v%x - p%x, v%y - p%y, v%z - p%z)
  end function subvp

  pure type(vector3d) function subpv(p, v)
    type(point3d), intent(in) :: p
    type(vector3d), intent(in) :: v
    subpv = vector3d(p%x - v%x, p%y - v%y, p%z - v%z)
  end function subpv

  pure type(vector3d) function subpp(a, b)
    type(point3d), intent(in) :: a, b
    subpp = vector3d(a%x - b%x, a%y - b%y, a%z - b%z)
  end function subpp


  !Definition of multiplication between reals and vectors.
  pure type(vector3d) function mulrv(r, v)
    real(real64), intent(in) :: r
    type(vector3d), intent(in) :: v
    mulrv = vector3d(r*v%x, r*v%y, r*v%z)
  end function mulrv

  pure type(vector3d) function mulvr(v, r)
    type(vector3d), intent(in) :: v
    real(real64), intent(in) :: r
    mulvr = vector3d(v%x*r, v%y*r, v%z*r)
  end function mulvr

  !Definition of multiplication between reals and points.
   pure type(point3d) function mulrp(r, p)
    real(real64), intent(in) :: r
    type(point3d), intent(in) :: p
    mulrp = point3d(r*p%x, r*p%y, r*p%z)
  end function mulrp

  pure type(point3d) function mulpr(p, r)
    type(point3d), intent(in) :: p
    real(real64), intent(in) :: r
    mulpr = point3d(p%x*r, p%y*r, p%z*r)
  end function mulpr


  !Definition of division between reals and vectors.
  pure type(vector3d) function divvr(v, r)
    type(vector3d), intent(in) :: v
    real(real64), intent(in) :: r
    if (r == 0) then
       divvr = vector3d(0.0, 0.0, 0.0)
    else
       divvr = vector3d(v%x/r, v%y/r, v%z/r)
    end if
  end function divvr

  pure type(point3d) function divpr(p, r)
    type(point3d), intent(in) :: p
    real(real64), intent(in) :: r
    if (r == 0) then
       divpr = point3d(0.0, 0.0, 0.0)
    else
       divpr = point3d(p%x/r, p%y/r, p%z/r)
    end if
  end function divpr

  !ULTIMATE DEBUGGER PRO 9000 (+18 only)
  function DEBUGGERPRO9000A(p) result(v)
    type(point3d), intent(in) :: p
    type(vector3d) :: v
    v%x = p%x
    v%y = p%y
    v%z = p%z
  end function DEBUGGERPRO9000A

  function DEBUGGERPRO9000B(v) result(p)
    type(vector3d), intent(in) :: v
    type(point3d) :: p
    p%x = v%x
    p%y = v%y
    p%z = v%z
  end function DEBUGGERPRO9000B


  !Definition of functions.
  pure type(vector3d) function distance(a, b)          !Calculates the distance vector between two 3d points.
    type(point3d), intent(in) :: a, b
    distance = vector3d(b%x - a%x, b%y - a%y, b%z - a%z)
  end function distance

  pure real(real64) function modulus(a)                !Calculates the magnitude of a vector.
    type(vector3d), intent(in) :: a
    modulus = sqrt((-a%x)**2 + (-a%y)**2 + (-a%z)**2)
  end function modulus
  
  pure type(vector3d) function normalize(a)            !Normalizes a vector.
    type(vector3d), intent(in) :: a
    normalize = divvr(a, modulus(a))
  end function normalize

  pure real(real64) function dotproduct(a, b)          !Calculates the dot product of two vectors.
    type(vector3d), intent(in) :: a, b
    dotproduct = a%x*b%x + a%y*b%y + a%z*b%z
  end function dotproduct

  pure real(real64) function angle(a, b)               !Calculates the angle between two vectors in radians.
    type(vector3d), intent(in) :: a, b
    angle = acos(dotproduct(normalize(a),  normalize(b)))
  end function angle

  pure type(vector3d) function cross_product(a, b)     !Calculates the cross product of two vectors.
    type(vector3d), intent(in) :: a, b
    cross_product = vector3d(a%y*b%z - a%z*b%y, a%z*b%x - a%x*b%z, a%x*b%y - a%y*b%x)
  end function cross_product
  
  pure type(vector3d) function orthv(a, b)             !Returns a vector orthogonal to two given vectors.
    type(vector3d), intent(in) :: a, b
    orthv = cross_product(a, b)
  end function orthv

  
end module geometry
