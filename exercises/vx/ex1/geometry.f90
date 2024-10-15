module geometry
  implicit none

  type :: vector3d
     real :: x, y, z
  end type vector3d

  type :: point3d
     real :: x, y, z
  end type point3d

  interface operator (.v+p.)
     module procedure sumvp
  end interface operator (.v+p.)

  interface operator (.p+v.)
     module procedure sumpv
  end interface operator (.p+v.)

  interface operator (.v-p.)
     module procedure subvp
  end interface operator (.v-p.)

  interface operator (.p-v.)
     module procedure subpv
  end interface operator (.p-v.)

  interface operator (.r*v.)
     module procedure mulrv
  end interface operator (.r*v.)

  interface operator (.v*r.)
     module procedure mulvr
  end interface operator (.v*r.)

  interface operator (.v/r.)
     module procedure divvr
  end interface operator (.v/r.)

  

contains
  pure type(vector3d) function sumvp(v, p)
    type(vector3d), intent(in) :: v
    type(point3d), intent(in) :: p
    sumvp = vector3d(v%x + p%x, v%y + p%y, v%z + p%z)
  end function sumvp

  pure type(vector3d) function sumpv(p, v)
    type(point3d), intent(in) :: p
    type(vector3d), intent(in) :: v
    sumpv = vector3d(p%x + v%x, p%y + v%y, p%z + v%z)
  end function sumpv

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

  pure type(vector3d) function mulrv(r, v)
    real, intent(in) :: r
    type(vector3d), intent(in) :: v
    mulrv = vector3d(r*v%x, r*v%y, r*v%z)
  end function mulrv

  pure type(vector3d) function mulvr(v, r)
    type(vector3d), intent(in) :: v
    real, intent(in) :: r
    mulvr = vector3d(v%x*r, v%y*r, v%z*r)
  end function mulrv

  pure type(vector3d) function divvr(v, r)
    type(vector3d), intent(in) :: v
    real, intent(in) :: r
    divvr = vector3d(v%x/r, v%y/r, v%z/r)
  end function divvr

  pure type(real) function distance(a, b)
    type(vector3d), intent(in) :: a, b
    distance = real(sqrt((b%x-a%x)**2+(b%y-a%y)**2+(b%z-a%z)**2))
  end function distance

  pure type(real) function modulus(a)
    type(vector3d), intent(in) :: a
    modulus = real(sqrt(a%x**2 + a%y**2 + a%z**2))
  end function modulus
  
  pure type(vector3d) function normalize(a)
    type(vector3d), intent(in) :: a
    normalize = divvr(a, modulus(a))
  end function normalize

  pure type(real) function dotproduct(a, b)
    type(vector3d), intent(in) :: a, b
    dotproduct = real(a%x*b%x + a%y*b%y + a%z*b%z)
  end function dotproduct

  pure type(real) function angle(a, b)
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
