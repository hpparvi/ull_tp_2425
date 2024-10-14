module particle
  use geometry
  implicit none
  type particle
     type(point3d) :: position
     type(vector3d) :: velocity
     real :: mass
  end type particle
end module particle
