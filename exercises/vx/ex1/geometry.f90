module geometry
  implicit none

  type :: vector3d
     real :: x, y, z
  end type vector3d

  type :: point3d
     real :: x, y, z
  end type point3d

  interface operator (.vsp.)
     module procedure sumvp
  end interface operator (.vsp.)

  interface operator (.psv.)
     module procedure sumpv
  end interface operator (.psv.)

  interface operator (.vrp.)
     module procedure subvp
  end interface operator (.vrp.)

  interface operator (.prv.)
     module procedure subpv
  end interface operator (.prv.)

  interface operator (.rpv.)
     module procedure mulrv
  end interface operator (.rpv.)

  interface operator (.vpr.)
     module procedure mulvr
  end interface operator (.vpr.)

  interface operator (.ver.)
     module procedure divvr
  end interface operator (.ver.)

  interface operator (.vsv.)                          !also added
     module procedure sumvv
  end interface operator (.vsv.)

  interface operator (.vrv.)                          !also added
     module procedure subvv
  end interface operator (.vrv.)


  

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

  pure type(vector3d) function subvv(a, b)            !added function
    type(vector3d), intent(in) :: a, b
    subvv = vector3d(a%x - b%x, a%y - b%y, a%z - b%z)
  end function subvv

  pure type(point3d) function subvp(v, p)
    type(vector3d), intent(in) :: v
    type(point3d), intent(in) :: p
    subvp = point3d(v%x - p%x, v%y - p%y, v%z - p%z)
  end function subvp

  pure type(point3d) function subpv(p, v)
    type(point3d), intent(in) :: p
    type(vector3d), intent(in) :: v
    subpv = point3d(p%x - v%x, p%y - v%y, p%z - v%z)
  end function subpv

  pure type(vector3d) function mulrv(r, v)
    real, intent(in) :: r
    type(vector3d), intent(in) :: v
    mulrv = vector3d(r*v%x, r*v%y, r*v%z)
  end function mulrv

  pure type(vector3d) function mulvr(v, r)
    type(vector3d), intent(in) :: v
    real, intent(in) :: r
    mulvr = vector3d(v%x*r, v%y*r, v%z*r)
  end function mulvr

  pure type(vector3d) function divvr(v, r)
    type(vector3d), intent(in) :: v
    real, intent(in) :: r
    divvr = vector3d(v%x/r, v%y/r, v%z/r)
  end function divvr

  pure type(vector3d) function distance(a, b)
    type(point3d), intent(in) :: a, b
    distance = vector3d(b%x - a%x, b%y - a%y, b%z - a%z)
  end function distance

  pure real function modulus(a)
    type(vector3d), intent(in) :: a
    modulus = sqrt(a%x**2 + a%y**2 + a%z**2)
  end function modulus
  
  pure type(vector3d) function normalize(a)
    type(vector3d), intent(in) :: a
    normalize = divvr(a, modulus(a))
  end function normalize

  pure real function dotproduct(a, b)
    type(vector3d), intent(in) :: a, b
    dotproduct = a%x*b%x + a%y*b%y + a%z*b%z
  end function dotproduct

  pure real function angle(a, b)
    type(vector3d), intent(in) :: a, b
    angle = acos(dotproduct(normalize(a),  normalize(b)))
  end function angle

  pure type(vector3d) function cross_product(a, b)
    type(vector3d), intent(in) :: a, b
    cross_product = vector3d(a%y*b%z - a%z*b%y, a%z*b%x - a%x*b%z, a%x*b%y - a%y*b%x)
  end function cross_product
  
  pure type(vector3d) function orthv(a, b)
    type(vector3d), intent(in) :: a, b
    orthv = cross_product(a, b)
  end function orthv

  
end module geometry
