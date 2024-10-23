module geometry
  use iso_fortran_env
  implicit none

  type :: vector3d !3d vector
     real(real64) :: x, y, z
  end type vector3d

  type :: point3d !3d point
     real(real64) :: x, y, z
  end type point3d

  interface operator(+) 
     module procedure sumvp, sumpv, sumvv
  end interface operator(+)

  interface operator(-)
     module procedure subvp, subpv, subvv, subpp
  end interface operator(-)

  interface operator(*)
     module procedure mulrv, mulvr
  end interface operator(*)

  interface operator(/)
     module procedure divvr
  end interface operator(/)

contains
  pure type(point3d) function sumvp(v, p) !adds a 3d vector and a 3d point
    type(vector3d), intent(in) :: v
    type(point3d), intent(in) :: p
    sumvp = point3d(v%x + p%x, v%y + p%y, v%z + p%z)
  end function sumvp

  pure type(point3d) function sumpv(p, v) !adds a 3d point and a 3d vector
    type(point3d), intent(in) :: p
    type(vector3d), intent(in) :: v
    sumpv = point3d(p%x + v%x, p%y + v%y, p%z + v%z)
  end function sumpv

  pure type(vector3d) function sumvv(v1, v2) !adds two 3d vectors
    type(vector3d), intent(in) :: v1, v2
    sumvv = vector3d(v1%x + v2%x, v1%y + v2%y, v1%z + v2%z)
  end function sumvv

  pure type(point3d) function subvp(v, p) !substracts a 3d vector and a 3d point
    type(vector3d), intent(in) :: v
    type(point3d), intent(in) :: p
    subvp = point3d(v%x - p%x, v%y - p%y, v%z - p%z)
  end function subvp

  pure type(point3d) function subpv(p, v) !substracts a 3d point and a 3d vector
    type(point3d), intent(in) :: p
    type(vector3d), intent(in) :: v
    subpv = point3d(p%x - v%x, p%y - v%y, p%z - v%z)
  end function subpv

  pure type(vector3d) function subvv(v1, v2) !substracts two 3d vectors
    type(vector3d), intent(in) :: v1, v2
    subvv = vector3d(v1%x - v2%x, v1%y - v2%y, v1%z - v2%z)
  end function subvv

  pure type(vector3d) function subpp(p1, p2) !substracts two 3d points
    type(point3d), intent(in) :: p1, p2
    subpp = vector3d(p1%x - p2%x, p1%y - p2%y, p1%z - p2%z)
  end function subpp

  pure type(vector3d) function mulrv(r, v) !multiplies a real number by a 3d vector
    real(real64), intent(in) :: r
    type(vector3d), intent(in) :: v
    mulrv = vector3d(r*v%x, r*v%y, r*v%z)
  end function mulrv

  pure type(vector3d) function mulvr(v, r) !multiplies a 3d vector by a real number
    type(vector3d), intent(in) :: v
    real(real64), intent(in) :: r
    mulvr = vector3d(v%x*r, v%y*r, v%z*r)
  end function mulvr

  pure type(vector3d) function divvr(v, r) !divides a 3d vector by a real number
    type(vector3d), intent(in) :: v
    real(real64), intent(in) :: r
    if (r == 0.0) then
       divvr = vector3d(0.0, 0.0, 0.0) !Returns a null 3d vector in case of dividing by zero
    else
       divvr = vector3d(v%x/r, v%y/r, v%z/r)
    end if   
  end function divvr

  pure real(real64) function distance(p1, p2) !distance between two 3d points
    type(point3d), intent(in) :: p1, p2
    distance = sqrt((p2%x - p1%x)**2 + (p2%y - p1%y)**2 + (p2%z - p1%z))
  end function distance
    
  pure real(real64) function norm(v) !norm of a vector
    type(vector3d), intent(in) :: v
    norm = sqrt(v%x**2 + v%y**2 + v%z**2)
  end function norm

  pure type(vector3d) function normalize(v) !normalize a 3d vector
    type(vector3d), intent(in) :: v
    if (norm(v) == 0.0) then
       normalize = vector3d(0.0, 0.0, 0.0) !Returns a null 3d vector in case of normalizing a zero norm vector
    else
       normalize = divvr(v,  norm(v))
    end if
  end function normalize

  pure real(real64) function dot_product(v1, v2) !multiplies two 3d vectors
    type(vector3d), intent(in) :: v1, v2
    dot_product  = v1%x * v2%x + v1%y * v2%y + v1%z * v2%z
  end function dot_product

  pure real(real64) function angle(v1, v2) !angle (in radians) between two 3d vectors 
    type(vector3d), intent(in) :: v1, v2
    angle = acos(dot_product(normalize(v1), normalize(v2)))
  end function angle
  
  pure type(vector3d) function cross_product(v1, v2) !cross product of two 3d vectors
    type(vector3d), intent(in) :: v1, v2
    cross_product = vector3d(v1%y * v2%z - v2%y * v1%z, v1%z * v2%x - v1%x * v2%z, v1%x * v2%y - v2%x * v1%y)
  end function cross_product

  pure type(vector3d) function orthv(v1, v2) !orthogonal vector of two 3d vectors
    type(vector3d), intent(in) :: v1, v2
    orthv = cross_product(v1, v2)
  end function orthv

end module geometry

