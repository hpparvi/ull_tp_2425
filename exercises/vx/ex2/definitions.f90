module definitions
  
  use iso_fortran_env
  use geometry
  implicit none
  
  TYPE RANGE
     REAL(real64), DIMENSION(3) :: min,max
  END TYPE RANGE
  
  TYPE CPtr
     TYPE(CELL), POINTER :: ptr
  END TYPE CPtr
  
  TYPE CELL
     TYPE(RANGE) :: range
     REAL(real64), DIMENSION(3) :: part
     INTEGER(int64) :: pos
     INTEGER(int64) :: type !! 0 = no particle; 1 = particle; 2 = conglomerado
     REAL(real64) :: mass
     REAL(real64), DIMENSION(3) :: c_o_m
     TYPE (CPtr), DIMENSION(2,2,2) :: subcell
  END TYPE CELL

  type ::  particle3d     !Definition of the particle.
     type(point3d) :: p   !Postion.
     type(vector3d) :: v  !Velocity.
     real(real64) :: m    !Mass.
  end type particle3d

end module definitions
