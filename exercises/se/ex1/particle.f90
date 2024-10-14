MODULE particles
  USE geometry
  IMPLICIT NONE

  TYPE particle
     TYPE(point3d)  :: p ! Position
     TYPE(vector3d) :: v ! Velocity
     REAL           :: m ! Mass
  END TYPE particle


END MODULE particles
