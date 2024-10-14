module geometry
  ! module that contains
  ! - data types: vector3d, point3d
  ! - functions to sum, substract, multiply and divide among them
  ! - operators associated to those functions
  ! - personal functions:
  !   + norm of a vector
  !   + vector between two points
  !   + dot product of two vectors
  ! - some other functions:
  !   + distance between points
  !   + angle between two vectors (in radians)
  !   + normalize a vector
  !   + cross product of vectors
  !   + orthogonal vector of two given ones
  
  implicit none

  type :: vector3d
     double precision :: x, y, z
  end type vector3d

  type :: point3d
     double precision :: x, y, z
  end type point3d

  interface operator(+)
     module procedure sumvp, sumpv, sumvv
  end interface

  interface operator(-)
     module procedure subvp, subpv, subvv
  end interface

  interface operator(*)
     module procedure mulrv, mulvr
  end interface

  interface operator(/)
     module procedure divvr
  end interface

  interface operator(.dot.)
     module procedure dotvv
  end interface

contains
  
! Adding vectors and points
  elemental type(point3d) function sumvp(av, bp)
    type(vector3d), intent(in) :: av
    type(point3d), intent(in)  :: bp
    sumvp = point3d(av%x + bp%x, av%y + bp%y, av%z + bp%z)
  end function sumvp

  elemental type(point3d) function sumpv(ap, bv)
    type(vector3d), intent(in) :: bv
    type(point3d), intent(in)  :: ap
    sumpv = point3d(ap%x + bv%x, ap%y + bv%y, ap%z + bv%z)
  end function sumpv

  elemental type(vector3d) function sumvv(a, b)
    type(vector3d), intent(in) :: a, b
    sumvv = vector3d(a%x + b%x, a%y + b%y, a%z + b%z)
  end function sumvv
  
! Substracting vectors and points
  elemental type(point3d) function subvp(av, bp)
    type(vector3d), intent(in) :: av
    type(point3d), intent(in)  :: bp
    subvp = point3d(av%x - bp%x, av%y - bp%y, av%z - bp%z)
  end function subvp

  elemental type(point3d) function subpv(ap, bv)
    type(vector3d), intent(in) :: bv
    type(point3d), intent(in)  :: ap
    subpv = point3d(ap%x - bv%x, ap%y - bv%y, ap%z - bv%z)
  end function subpv

  elemental type(vector3d) function subvv(a, b)
    type(vector3d), intent(in) :: a, b
    subvv = vector3d(a%x - b%x, a%y - b%y, a%z - b%z)
  end function subvv
  
! Multiplication of real number and vector
  elemental type(vector3d) function mulrv(a, b)
    double precision, intent(in) :: a
    type(vector3d), intent(in) :: b
    mulrv = vector3d(a*b%x, a*b%y, a*b%z)
  end function mulrv

  elemental type(vector3d) function mulvr(a, b)
    type(vector3d), intent(in) :: a
    double precision, intent(in) :: b
    mulvr = vector3d(b*a%x, b*a%y, b*a%z)
  end function mulvr

!  elemental function mulvr(a, b)
!    type(vector3d), intent(in) :: a
!    double precision, intent(in) :: b
!    mulvr = vector3d(b*a%x, b*a%y, b*a%z)
!  end function mulvr

! Division of vector by real
  elemental type(vector3d) function divvr(v, r)
    double precision, intent(in) :: r
    type(vector3d), intent(in) :: v
    divvr = vector3d(v%x/r, v%y/r, v%z/r)
  end function divvr

! Norm of a vector
  pure double precision function norm(a)
    type(vector3d), intent(in) :: a
    norm = sqrt(a%x**2 + a%y**2 + a%z**2)
  end function norm

! Vector from points
  pure type(vector3d) function vecpp(a,b)
    type(point3d), intent(in)  :: a, b
    vecpp = vector3d(b%x-a%x, b%y-a%y, b%z-a%z)
  end function vecpp

! Dot product of two vectors
  pure double precision function dotvv(a,b)
    type(vector3d), intent(in)  :: a, b
    dotvv = a%x*b%x +a%y*b%y + a%z*b%z
  end function dotvv
    
! Distance calculation between two points,
!i.e. norm of the vector that joins the two points
  pure double precision function distance(a, b)
    type(point3d), intent(in) :: a, b
    distance = norm(vecpp(a,b))
  end function distance

! Angle between two vectors (in radians)
  pure double precision function angle(a,b)
    type(vector3d), intent(in)  :: a, b
    angle = acos(a.dot.b / (norm(a) * norm(b)))
  end function angle

! normalization of vector
  pure type(vector3d) function normalize(a)
    type(vector3d), intent(in)  :: a
    normalize = a/norm(a)
  end function normalize

! Cross product of two vectors
  pure type(vector3d) function cross_product(a, b)
    type(vector3d), intent(in)  :: a,b
    cross_product = vector3d(a%y*b%z-a%z*b%y, -a%x*b%z+a%z*b%x, a%x*b%y-a%y*b%x)
  end function cross_product

! Returns an orthogonal normalized vector
  pure type(vector3d) function orthv(a,b)
    type(vector3d), intent(in)  :: a,b
    orthv = normalize(cross_product(a,b))
  end function orthv
  
end module geometry
