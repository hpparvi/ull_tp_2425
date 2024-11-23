MODULE EX2_deriv_types
  USE geometry
  USE ISO_FORTRAN_ENV
  
  TYPE particle3d
    TYPE(point3d)  :: p ! position
	  TYPE(vector3d) :: v ! velocity
    REAL(REAL64) :: m           ! mass
  END TYPE particle3d
  
  ! Define a interval between two point
  TYPE RANGE
    TYPE(point3d)  :: min,max
  END TYPE RANGE
  
  TYPE CPtr
    TYPE(CELL), POINTER :: ptr
  END TYPE CPtr
  
  TYPE CELL
    TYPE (RANGE) :: range
    INTEGER :: pos
    INTEGER :: type !! 0 = no particle; 1 = particle; 2 = conglomerado
    TYPE(particle3d) :: pt 
    TYPE(vector3d)  :: c_o_m ! center of mase. Is a vector!
    TYPE (CPtr), DIMENSION(2,2,2) :: subcell 
  END TYPE CELL
  
END MODULE EX2_deriv_types