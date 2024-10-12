module particle
  use geometry
  implicit none
  type particle
     type(point3d) :: p !osition
     type(vector3d) :: v !elocity
     real :: m !ass
  end type particle
end module particle
