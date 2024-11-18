MODULE particle
  USE, INTRINSIC :: iso_fortran_env
  USE geometry
  IMPLICIT NONE

  TYPE particle3d
     TYPE(point3d)  :: p ! Position
     TYPE(vector3d) :: v ! Velocity
     REAL(real64) :: m ! Mass
  END TYPE particle3d


END MODULE particle
