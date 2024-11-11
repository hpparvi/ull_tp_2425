MODULE particle
  USE geometry
  USE ISO_FORTRAN_ENV
  TYPE particle3d
    TYPE(point3d)  :: p ! position
	  TYPE(vector3d) :: v ! velocity
    REAL(REAL64) :: m           ! mass
  END TYPE particle3d
END MODULE particle