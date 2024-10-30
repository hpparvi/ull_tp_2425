module particle
  use iso_fortran_env
  use geometry
  implicit none
  
  type particle3d
    type(point3d) :: p ! particle position
    type(vector3d) :: v ! velocity of the particle
    real(real64) :: m  ! mass of the particle
  end type particle3d
  
end module particle
