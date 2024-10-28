module particle
  use geometry !Importing the geometry module
  use iso_fortran_env !Importing iso_fortran_env to specify the number of bits of variables
  implicit none

  !Define a type for a 3D particle with position, velocity, acceleration and mass
  type :: particle3d
     type(point3d) :: p  !Position of particle
     type(vector3d) :: v !Velocity of particle
     type(vector3d) :: a !Acceleration of particle
     real(real64) :: m   !Mass of particle
  end type particle3d

end module particle
