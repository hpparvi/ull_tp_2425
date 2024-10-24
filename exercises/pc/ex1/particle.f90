module particle
  use geometry
  implicit none
  
  type :: particle3d
     type(point3d) :: p !particle's position
     type(vector3d) :: v !particle's velocity
     type(vector3d) :: a !particle's acceleration
     real(real64) :: m !particles's mass
  end type particle3d
  
end module particle
