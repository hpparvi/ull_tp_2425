module particl
  use geometry
  implicit none
  
  type ::  particle
     type(point3d) :: p
     type(vector3d) :: v  !vector3d instead of 2d
     type(vector3d) :: a  !needed for next part
     real :: m
  end type particle

end module particl
