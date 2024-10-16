module particle
  use geometry
  implicit none
  type particle3d
     type(point3d) :: p !position
     type(vector3d) :: v !velocity
     real :: m !mass
  end type particle3d
end module particle
