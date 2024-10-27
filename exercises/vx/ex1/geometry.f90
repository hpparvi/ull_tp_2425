module geometry
  use iso_fortran_env
  implicit none
  
  
  type :: vector3d
     real(real64) :: x, y, z
  end type vector3d

  type :: point3d
     real(real64) :: x, y, z
  end type point3d

  interface operator (+)
     module procedure sumvp, sumpv, sumpp, sumvv
  end interface operator (+)

  interface operator (-)
     module procedure subvp, subpv, subpp, subvv
  end interface operator (-)

  interface operator (*)
     module procedure mulrv, mulvr
  end interface operator (*)

  interface operator (/)
     module procedure divvr
  end interface operator (/)
  
  

contains
  pure type(point3d) function sumvp(v, p)            !changed point + vector = point
    type(vector3d), intent(in) :: v
    type(point3d), intent(in) :: p
    sumvp = point3d(v%x + p%x, v%y + p%y, v%z + p%z)
  end function sumvp

  pure type(point3d) function sumpv(p, v)             !same change as ^
    type(point3d), intent(in) :: p
    type(vector3d), intent(in) :: v
    sumpv = point3d(p%x + v%x, p%y + v%y, p%z + v%z)
  end function sumpv

  pure type(vector3d) function sumvv(a, b)            !added function
    type(vector3d), intent(in) :: a, b
    sumvv = vector3d(a%x + b%x, a%y + b%y, a%z + b%z)
  end function sumvv

  pure type(vector3d) function sumpp(a, b)
    type(point3d), intent(in) :: a, b
    sumpp = vector3d(a%x + b%x, a%y + b%y, a%z + b%z)
  end function sumpp

  pure type(vector3d) function subvv(a, b)            !added function
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

  pure type(vector3d) function divvr(v, r)
    type(vector3d), intent(in) :: v
    real(real64), intent(in) :: r
    divvr = vector3d(v%x/r, v%y/r, v%z/r)
  end function divvr

  pure type(vector3d) function distance(a, b)
    type(point3d), intent(in) :: a, b
    distance = vector3d(b%x - a%x, b%y - a%y, b%z - a%z)
  end function distance

  pure real(real64) function modulus(a)
    type(vector3d), intent(in) :: a
    modulus = sqrt((-a%x)**2 + (-a%y)**2 + (-a%z)**2)
  end function modulus
  
  pure type(vector3d) function normalize(a)
    type(vector3d), intent(in) :: a
    normalize = divvr(a, modulus(a))
  end function normalize

  pure real(real64) function dotproduct(a, b)
    type(vector3d), intent(in) :: a, b
    dotproduct = a%x*b%x + a%y*b%y + a%z*b%z
  end function dotproduct

  pure real(real64) function angle(a, b)
    type(vector3d), intent(in) :: a, b
    angle = acos(dotproduct(normalize(a),  normalize(b)))
  end function angle

  pure type(vector3d) function cross_product(a, b)
    type(vector3d), intent(in) :: a, b
    cross_product = vector3d(a%y*b%z - a%z*b%y, a%z*b%x - a%x*b%z, a%x*b%y - a%y*b%x)
  end function cross_product
  
!  pure type(vector3d) function orthv(a, b)
!    type(vector3d), intent(in) :: a, b
!    orthv = vector3d
!  end function orthv

  
end module geometry
