module particle        !This module defines the particle structure with position, velocity and mass.
  use iso_fortran_env  !This module ensures all variables are defined as 64-bit.
  use geometry         !This module defines 3D vector and point operations for vector3d and point3d types.
  implicit none
  
  type ::  particle3d     !Definition of the particle.
     type(point3d) :: p   !Postion.
     type(vector3d) :: v  !Velocity.
     real(real64) :: m    !Mass.
  end type particle3d

end module particle
