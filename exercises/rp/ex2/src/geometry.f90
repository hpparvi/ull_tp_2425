module geometry
   implicit none

   type :: vector3d
      real :: x, y, z
   end type vector3d

   type :: point3d
      real :: x, y, z
   end type point3d

   interface operator (+)
      module procedure sumvp  ! vector + point
   end interface
   interface operator (+)
      module procedure sumpv  ! point + vector
   end interface
   interface operator (+)
      module procedure sumpp  ! point + point
   end interface
   interface operator (+)
      module procedure sumvv  ! point + vector
   end interface
   interface operator (-)
      module procedure subvp  ! vector - point
   end interface
   interface operator (-)
      module procedure subpv  ! point - vector
   end interface
   interface operator (-)
      module procedure subpp  ! point - point
   end interface
   interface operator (-)
      module procedure subvv  ! point + vector
   end interface
   interface operator (*)
      module procedure mulrv  ! constant * vector
   end interface
   interface operator (*)
      module procedure mulvr  ! vector * constant
   end interface
   interface operator (*)
      module procedure mulrp  ! constant * point
   end interface
   interface operator (*)
      module procedure mulpr  ! point * constant
   end interface
   interface operator (/)
      module procedure divvr  ! vector / constant
   end interface

contains

   pure function sumvp(v, p) result(r)
      ! vector + point
      type(vector3d), intent(in) :: v
      type(point3d), intent(in) :: p
      type(point3d) :: r

      r%x = p%x + v%x
      r%y = p%y + v%y
      r%z = p%z + v%z
   end function sumvp

   pure function sumpv(p, v) result(r)
      ! point + vector
      type(vector3d), intent(in) :: v
      type(point3d), intent(in) :: p
      type(point3d) :: r

      r%x = p%x + v%x
      r%y = p%y + v%y
      r%z = p%z + v%z
   end function sumpv

   pure function sumvv(v1, v2) result(r)
      ! vector + vector
      type(vector3d), intent(in) :: v1
      type(vector3d), intent(in) :: v2
      type(vector3d) :: r

      r%x = v1%x + v2%x
      r%y = v1%y + v2%y
      r%z = v1%z + v2%z
   end function sumvv

   pure function sumpp(p1, p2) result(r)
      ! point + point
      type(point3d), intent(in) :: p1
      type(point3d), intent(in) :: p2
      type(point3d) :: r

      r%x = p1%x + p2%x
      r%y = p1%y + p2%y
      r%z = p1%z + p2%z
   end function sumpp

   pure function subvp(v, p) result(r)
      ! vector - point
      type(vector3d), intent(in) :: v
      type(point3d), intent(in) :: p
      type(point3d) :: r

      r%x = v%x - p%x
      r%y = v%y - p%y
      r%z = v%z - p%z
   end function subvp

   pure function subpv(p, v) result(r)
      ! point - vector
      type(vector3d), intent(in) :: v
      type(point3d), intent(in) :: p
      type(point3d) :: r

      r%x = p%x - v%x
      r%y = p%y - v%y
      r%z = p%z - v%z
   end function subpv


   pure function subpp(p1, p2) result(r)
      ! point - point
      type(point3d), intent(in) :: p1
      type(point3d), intent(in) :: p2
      type(vector3d) :: r

      r%x = p1%x - p2%x
      r%y = p1%y - p2%y
      r%z = p1%z - p2%z
   end function subpp


   pure function subvv(v1, v2) result(r)
      ! point - point
      type(vector3d), intent(in) :: v1
      type(vector3d), intent(in) :: v2
      type(vector3d) :: r

      r%x = v1%x - v2%x
      r%y = v1%y - v2%y
      r%z = v1%z - v2%z
   end function subvv

   pure function mulrv(s, v) result(r)
      ! constant * vector
      type(vector3d), intent(in) :: v
      real, intent(in) :: s
      type(vector3d) :: r

      r%x = s * v%x
      r%y = s * v%y
      r%z = s * v%z
   end function mulrv

   pure function mulvr(v, s) result(r)
      ! vector * constant
      type(vector3d), intent(in) :: v
      real, intent(in) :: s
      type(vector3d) :: r

      r%x = s * v%x
      r%y = s * v%y
      r%z = s * v%z
   end function mulvr


   pure function mulrp(s, v) result(r)
      ! constant * vector
      type(point3d), intent(in) :: v
      real, intent(in) :: s
      type(point3d) :: r

      r%x = s * v%x
      r%y = s * v%y
      r%z = s * v%z
   end function mulrp

   pure function mulpr(v, s) result(r)
      ! vector * constant
      type(point3d), intent(in) :: v
      real, intent(in) :: s
      type(point3d) :: r

      r%x = s * v%x
      r%y = s * v%y
      r%z = s * v%z
   end function mulpr


   pure function divvr(v, s) result(r)
      ! vector / constant
      type(vector3d), intent(in) :: v
      real, intent(in) :: s
      type(vector3d) :: r

      r%x = v%x / s
      r%y = v%y / s
      r%z = v%z / s
   end function divvr

   pure function distance(p1, p2) result(r)
      ! Distance between points
      type(point3d), intent(in) :: p1
      type(point3d), intent(in) :: p2
      real :: r

      r = sqrt((p1%x - p2%x)**2+(p1%y - p2%y)**2+(p1%z - p2%z)**2)

   end function distance

   pure function dot(v1, v2) result(r)
      ! scalar product of two vectors
      type(vector3d), intent(in) :: v1
      type(vector3d), intent(in) :: v2
      real :: r
      r = (v1%x)*(v2%x) + (v1%y)*(v2%y) + (v1%z)*(v2%z)
   end function dot

   pure function norm(v) result(r)
      ! norm of a vector
      type(vector3d), intent(in) :: v
      real :: r
      r = sqrt(dot(v, v))
   end function norm

   pure function angle(v1, v2) result(theta)
      ! angle between two vectors
      type(vector3d), intent(in) :: v1
      type(vector3d), intent(in) :: v2
      real :: theta
      theta = acos(dot(v1, v2) / (norm(v1) * norm(v2)))
   end function angle

   pure function normalize(v) result(v_norm)
      ! normalize a vector
      type(vector3d), intent(in) :: v
      type(vector3d) :: v_norm
      v_norm = v / norm(v)
   end function normalize

   pure function cross_product(v1, v2) result(v_cross)
      ! cross product between two vectors
      type(vector3d), intent(in) :: v1
      type(vector3d), intent(in) :: v2
      type(vector3d) :: v_cross
      v_cross%x = (v1%y)*(v2%z) - (v1%z)*(v2%y)
      v_cross%y = (v1%z)*(v2%x) - (v1%x)*(v2%z)
      v_cross%z = (v1%x)*(v2%y) - (v1%y)*(v2%x)
   end function cross_product

   pure function orthv(v1, v2) result(v_orth)
      ! orthogonal vector
      type(vector3d), intent(in) :: v1
      type(vector3d), intent(in) :: v2
      type(vector3d) :: v_orth
      v_orth = cross_product(v1, v2)
      v_orth = v_orth / norm(v_orth)
   end function orthv

end module geometry
