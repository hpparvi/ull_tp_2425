MODULE particlepvm
  use geometry
  TYPE particle
    TYPE(point3d)  :: p ! position
	  TYPE(vector3d) :: v ! velocity
    REAL :: m           ! mass
  END TYPE particle
END MODULE particlepvm