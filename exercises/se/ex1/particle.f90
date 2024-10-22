MODULE particle
  USE geometry
  IMPLICIT NONE

  TYPE particle3d
     TYPE(point3d)  :: p ! Position
     TYPE(vector3d) :: v ! Velocity
     REAL :: m ! Mass
  END TYPE particle3d


END MODULE particle
