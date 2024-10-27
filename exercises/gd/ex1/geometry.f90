MODULE geometry
  
  USE ISO_FORTRAN_ENV
  ! Point and vector have different behavior in maths operations
  TYPE vector3d
    REAL(REAL64) :: x,y,z
  END TYPE vector3d

  TYPE point3d
    REAL(REAL64) :: x,y,z
  END TYPE point3d

  INTERFACE operator(+)
    MODULE procedure sumvv, sumvp, sumpv
  END INTERFACE

  INTERFACE operator(-)
    MODULE procedure subvv, subvp, subpv, subpp
  END INTERFACE
  
  INTERFACE operator(*)
    MODULE procedure mulrv, mulvr, dotmul
  END INTERFACE

  INTERFACE operator(/)
    MODULE procedure divpr, divvr
  END INTERFACE
  
CONTAINS
  
  !sum of vectors is vector
  TYPE(vector3d) ELEMENTAL FUNCTION sumvv(v1,v2)
    IMPLICIT NONE
    TYPE(vector3d), INTENT(IN) :: v1, v2
    sumvv%x = v1%x + v2%x
    sumvv%y = v1%y + v2%y
    sumvv%z = v1%z + v2%z
  END FUNCTION sumvv
  
  !Sum of vector ad a point is a translation (of the point)
  TYPE(point3d) ELEMENTAL FUNCTION sumvp(v,p)
    IMPLICIT NONE
    TYPE(point3d), INTENT(IN) :: p
    TYPE(vector3d), INTENT(IN) :: v
    sumvp%x = p%x + v%x
    sumvp%y = p%y + v%y
    sumvp%z = p%z + v%z
  END FUNCTION sumvp
  
  TYPE(point3d) ELEMENTAL FUNCTION sumpv(p,v)
    IMPLICIT NONE
    TYPE(point3d), INTENT(IN) :: p
    TYPE(vector3d), INTENT(IN) :: v
    sumpv%x = p%x + v%x
    sumpv%y = p%y + v%y
    sumpv%z = p%z + v%z
  END FUNCTION sumpv
  !two point cannot be added
  
  !Subs of two vector is a vector
  TYPE(vector3d) ELEMENTAL FUNCTION subvv(v1,v2)
    IMPLICIT NONE
    TYPE(vector3d), INTENT(IN) :: v1, v2
    subvv%x = v1%x - v2%x
    subvv%y = v1%y - v2%y
    subvv%z = v1%z - v2%z
  END FUNCTION subvv
  
  !Sub of vector ad a point is a translation (of the point)
  TYPE(point3d) ELEMENTAL FUNCTION subvp(v,p)
    IMPLICIT NONE
    TYPE(point3d), INTENT(IN) :: p
    TYPE(vector3d), INTENT(IN) :: v
    subvp%x = v%x - p%x
    subvp%y = v%y - p%y
    subvp%z = v%z - p%z
  END FUNCTION subvp
  
  TYPE(point3d) ELEMENTAL FUNCTION subpv(p,v)
    IMPLICIT NONE
    TYPE(point3d), INTENT(IN) :: p
    TYPE(vector3d), INTENT(IN) :: v
    subpv%x = p%x - v%x
    subpv%y = p%y - v%y
    subpv%z = p%z - v%z
  END FUNCTION subpv
  
  !Subs of two points is the vector that joins
  TYPE(vector3d) ELEMENTAL FUNCTION subpp(p1,p2)
    IMPLICIT NONE
    TYPE(point3d), INTENT(IN) :: p1, p2
    subpp%x = p2%x - p1%x
    subpp%y = p2%y - p1%y
    subpp%z = p2%z - p1%z
  END FUNCTION subpp
  
  !Mults and division are done coordinate to coordinate
  TYPE(vector3d) ELEMENTAL FUNCTION mulrv(r,v)
    REAL(REAL64), INTENT(IN) :: r
    TYPE(vector3d), INTENT(IN) :: v
    mulrv%x = r * v%x
    mulrv%y = r * v%y
    mulrv%z = r * v%z
  END FUNCTION mulrv
  
  TYPE(vector3d) ELEMENTAL FUNCTION mulvr(v,r)
    REAL(REAL64), INTENT(IN) :: r
    TYPE(vector3d), INTENT(IN) :: v
    mulvr%x = r * v%x
    mulvr%y = r * v%y
    mulvr%z = r * v%z
  END FUNCTION mulvr
  
  TYPE(vector3d) ELEMENTAL FUNCTION divpr(p,r)
    REAL(REAL64), INTENT(IN) :: r
    TYPE(point3d), INTENT(IN) :: p
    divpr%x = p%x / r
    divpr%y = p%y / r
    divpr%z = p%z / r
  END FUNCTION divpr
  
  TYPE(vector3d) ELEMENTAL FUNCTION divvr(v,r)
    REAL(REAL64), INTENT(IN) :: r
    TYPE(vector3d), INTENT(IN) :: v
    divvr%x = v%x / r
    divvr%y = v%y / r
    divvr%z = v%z / r
  END FUNCTION divvr
  
  ! definition of the norm of a vector
  REAL(REAL64) ELEMENTAL FUNCTION norm(v)
    TYPE(vector3d), INTENT(IN) :: v
    norm = SQRT(v%x**2 + v%y**2 + v%z**2)
  END FUNCTION norm
  
  !Distance between two point is the norm of the vector that separates them
  REAL(REAL64) ELEMENTAL FUNCTION distance(p1,p2)
    TYPE(point3d), INTENT(IN) :: p1, p2
    distance = norm(p1-p2) !SQRT((p1%x-p2%x)**2 + (p1%y-p2%y)**2 + (p1%z-p2%z)**2)
  END FUNCTION distance
  
  ! Definition of dot product
  REAL(REAL64) ELEMENTAL FUNCTION dotmul(v1,v2)
    TYPE(vector3d), INTENT(IN) :: v1, v2
    dotmul = v1%x * v2%x + v1%y * v2%y + v1%z * v2%z
  END FUNCTION dotmul
  
  ! Sefinition of angle between two vector
  REAL(REAL64) ELEMENTAL FUNCTION angle(v1, v2)
    TYPE(vector3d), INTENT(IN) :: v1, v2
    angle = ACOS((v1*v2)/(norm(v1)*norm(v2)))
  END FUNCTION angle
  
  ! Normalice is div the vector by its norm
  TYPE(vector3d) ELEMENTAL FUNCTION normalize(v)
    TYPE(vector3d), INTENT(IN) :: v
    normalize = v/norm(v)
  END FUNCTION normalize
  
  ! Definition of cross product
  TYPE(vector3d) ELEMENTAL FUNCTION cross_product(v1,v2)
    TYPE(vector3d), INTENT(IN) :: v1, v2
    cross_product%x = v1%y*v2%z - v1%z*v2%y
    cross_product%y = v1%z*v2%x - v1%x*v2%z
    cross_product%z = v1%x*v2%y - v1%y*v2%x
  END FUNCTION cross_product
  
  ! normalize the cross product
  TYPE(vector3d) ELEMENTAL FUNCTION orthv(v1,v2)
    TYPE(vector3d), INTENT(IN) :: v1, v2
    orthv = normalize(cross_product(v1,v2))
  END FUNCTION orthv
  
END module