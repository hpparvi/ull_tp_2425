module particle
  use geometry
  implicit none
  
  type ::  particle3d
     type(point3d) :: p
     type(vector3d) :: v  !vector3d instead of 2d
     type(vector3d) :: a  !needed for next part
     real :: m
  end type particle3d

end module particle
