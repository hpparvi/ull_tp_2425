module geometry
  use iso_fortran_env
  implicit none

  !Define a type for a 3D vector with components x, y and z
  type :: vector3d
     real(real64) :: x, y, z
  end type vector3d

  !Define a type for a 3D point with coordinates x, y and z
  type :: point3d
     real(real64) :: x, y, z
  end type point3d

  !Overload operator (+) for adding vectors and points
  interface operator(+) 
     module procedure sumvp, sumpv, sumvv
  end interface operator(+)

  !Overload operator (-) for subtracting vectors and points
  interface operator(-)
     module procedure subvp, subpv, subvv, subpp
  end interface operator(-)

  !Overload operator (*) for multiplying a vector by a scalar
  interface operator(*)
     module procedure mulrv, mulvr
  end interface operator(*)

  !Overload operator (/) for dividing a vector by a scalar
  interface operator(/)
     module procedure divvr
  end interface operator(/)

contains

  !Function to add a 3D vector to a 3D point, resulting in a new point
  pure type(point3d) function sumvp(v, p)
    type(vector3d), intent(in) :: v
    type(point3d), intent(in) :: p
    sumvp = point3d(v%x + p%x, v%y + p%y, v%z + p%z)
  end function sumvp

  !Function to add a 3D point to a 3D vector, resulting in a new point
  pure type(point3d) function sumpv(p, v)
    type(point3d), intent(in) :: p
    type(vector3d), intent(in) :: v
    sumpv = point3d(p%x + v%x, p%y + v%y, p%z + v%z)
  end function sumpv

  !Function to add two 3D vectors, resulting in a new vector
  pure type(vector3d) function sumvv(v1, v2)
    type(vector3d), intent(in) :: v1, v2
    sumvv = vector3d(v1%x + v2%x, v1%y + v2%y, v1%z + v2%z)
  end function sumvv

  !Function to subtract a 3D point from a 3D vector, resulting in a new point
  pure type(point3d) function subvp(v, p)
    type(vector3d), intent(in) :: v
    type(point3d), intent(in) :: p
    subvp = point3d(v%x - p%x, v%y - p%y, v%z - p%z)
  end function subvp

  !Function to subtract a 3D vector from a 3D point, resulting in a new point
  pure type(point3d) function subpv(p, v)
    type(point3d), intent(in) :: p
    type(vector3d), intent(in) :: v
    subpv = point3d(p%x - v%x, p%y - v%y, p%z - v%z)
  end function subpv

  !Function to subtract two 3D vectors, resulting in a new vector
  pure type(vector3d) function subvv(v1, v2)
    type(vector3d), intent(in) :: v1, v2
    subvv = vector3d(v1%x - v2%x, v1%y - v2%y, v1%z - v2%z)
  end function subvv

  !Function to subtract two 3D points, resulting in a new vector
  pure type(vector3d) function subpp(p1, p2)
    type(point3d), intent(in) :: p1, p2
    subpp = vector3d(p1%x - p2%x, p1%y - p2%y, p1%z - p2%z)
  end function subpp

  !Function to multiply a scalar by a 3D vector, resulting in a new vector
  pure type(vector3d) function mulrv(r, v)
    real(real64), intent(in) :: r
    type(vector3d), intent(in) :: v
    mulrv = vector3d(r*v%x, r*v%y, r*v%z)
  end function mulrv

  !Function to multiply a 3D vector by a scalar, resulting in a new vector
  pure type(vector3d) function mulvr(v, r)
    type(vector3d), intent(in) :: v
    real(real64), intent(in) :: r
    mulvr = vector3d(v%x*r, v%y*r, v%z*r)
  end function mulvr

  !Function to divide a 3D vector by a scalar, resulting in a new vector
  pure type(vector3d) function divvr(v, r)
    type(vector3d), intent(in) :: v
    real(real64), intent(in) :: r
    if (r == 0.0) then
       divvr = vector3d(0.0, 0.0, 0.0) !Return a null vector if division by zero
    else
       divvr = vector3d(v%x/r, v%y/r, v%z/r)
    end if   
  end function divvr

  !Function to calculate the distance between two 3D points
  pure real(real64) function distance(p1, p2)
    type(point3d), intent(in) :: p1, p2
    distance = sqrt((p2%x - p1%x)**2 + (p2%y - p1%y)**2 + (p2%z - p1%z))
  end function distance
    
  !Function to calculate the norm of a 3D vector
  pure real(real64) function norm(v)
    type(vector3d), intent(in) :: v
    norm = sqrt(v%x**2 + v%y**2 + v%z**2)
  end function norm

  !Function to normalize a 3D vector, returning a unit vector
  pure type(vector3d) function normalize(v)
    type(vector3d), intent(in) :: v
    if (norm(v) == 0.0) then
       normalize = vector3d(0.0, 0.0, 0.0) !Return a null vector if norm is zero
    else
       normalize = divvr(v,  norm(v))
    end if
  end function normalize

  !Function to calculate the dot product of two 3D vectors
  pure real(real64) function dot_product(v1, v2)
    type(vector3d), intent(in) :: v1, v2
    dot_product  = v1%x * v2%x + v1%y * v2%y + v1%z * v2%z
  end function dot_product

  !Function to calculate the angle (in radians) between two 3D vectors
  pure real(real64) function angle(v1, v2)
    type(vector3d), intent(in) :: v1, v2
    angle = acos(dot_product(normalize(v1), normalize(v2)))
  end function angle
  
  !Function to calculate the cross product of two 3D vectors
  pure type(vector3d) function cross_product(v1, v2)
    type(vector3d), intent(in) :: v1, v2
    cross_product = vector3d(v1%y * v2%z - v2%y * v1%z, v1%z * v2%x - v1%x * v2%z, v1%x * v2%y - v2%x * v1%y)
  end function cross_product

  !Function to calculate an orthogonal vector to two 3D vectors (using cross product)
  pure type(vector3d) function orthv(v1, v2)
    type(vector3d), intent(in) :: v1, v2
    orthv = cross_product(v1, v2)
  end function orthv

end module geometry
