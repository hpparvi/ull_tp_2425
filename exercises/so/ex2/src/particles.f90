module particles
  ! module that contains
  !   - data type: particle, with position, velocity and mass

  use geometry

  implicit none

  type :: particle
     type(point3d)     :: p  ! position
     type(vector3d)    :: v  ! velocity
     type(vector3d)    :: a  ! acceleration
     real  :: m  ! mass
  end type particle

end module particles
