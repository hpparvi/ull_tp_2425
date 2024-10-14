! This module contains the definitions of several
! useful operators

MODULE geometry
  IMPLICIT NONE

  ! Define vectors
  TYPE vector3d
     DOUBLE PRECISION :: x, y, z
  END TYPE vector3d

  ! Define points
  TYPE point3d
     DOUBLE PRECISION :: x, y, z
  END TYPE point3d

  ! Operators for code clarity
  INTERFACE OPERATOR(+)
     MODULE PROCEDURE sumvp, sumpv
  END INTERFACE

  INTERFACE OPERATOR(-)
     MODULE PROCEDURE subvp, subpv
  END INTERFACE

  INTERFACE OPERATOR(*)
     MODULE PROCEDURE mulrv, mulvr
  END INTERFACE

  INTERFACE OPERATOR(/)
     MODULE PROCEDURE divvr
  END INTERFACE


  
! Functions
CONTAINS
  
  ! Adding a vector and a point
  FUNCTION sumvp(vec, poi)
    TYPE(vector3d) :: sumvp
    TYPE(vector3d) :: vec
    TYPE(point3d)  :: poi

    sumvp%x = vec%x + poi%x
    sumvp%y = vec%y + poi%y
    sumvp%z = vec%z + poi%z
    
  END FUNCTION sumvp


  ! Adding a point and a vector
  ! (Fortran does not assume commutative property)
  FUNCTION sumpv(poi, vec)
    TYPE(vector3d) :: sumpv
    TYPE(vector3d) :: vec
    TYPE(point3d)  :: poi
    
    sumpv%x = vec%x + poi%x
    sumpv%y = vec%y + poi%y
    sumpv%z = vec%z + poi%z
    
  END FUNCTION sumpv

  
  ! Subtracting a vector minus a point
  FUNCTION subvp(vec, poi)
    TYPE(vector3d) :: subvp
    TYPE(vector3d) :: vec
    TYPE(point3d)  :: poi
    
    subvp%x = vec%x - poi%x
    subvp%y = vec%y - poi%y
    subvp%z = vec%z - poi%z
    
  END FUNCTION subvp

  
  ! Subtracting a point minus a vector
  FUNCTION subpv(poi, vec)
    TYPE(vector3d) :: subpv
    TYPE(vector3d) :: vec
    TYPE(point3d)  :: poi
    
    subpv%x = poi%x - vec%x
    subpv%y = poi%y - vec%y
    subpv%z = poi%z - vec%z
    
  END FUNCTION subpv

  
  ! Multiplying a real and a vector
  FUNCTION mulrv(re, vec)
    TYPE(vector3d) :: mulrv
    TYPE(vector3d) :: vec
    DOUBLE PRECISION :: re
    
    mulrv%x = re * vec%x
    mulrv%y = re * vec%y
    mulrv%z = re * vec%z
    
  END FUNCTION mulrv

  
  ! Multiplying a vector and a real
  FUNCTION mulvr(vec, re)
    TYPE(vector3d) :: mulvr
    TYPE(vector3d) :: vec
    DOUBLE PRECISION :: re
    
    mulvr%x = re * vec%x
    mulvr%y = re * vec%y
    mulvr%z = re * vec%z
    
  END FUNCTION mulvr

  
  ! Dividing a vector by a real
  FUNCTION divvr
    TYPE(vector3d) :: divvr
    TYPE(vector3d) :: vec
    DOUBLE PRECISION :: re
    
    divvr%x = vec%x / re
    divvr%y = vec%y / re
    divvr%z = vec%z / re
    
  END FUNCTION divvr


  ! Getting the magnitude of a vector
  FUNCTION magnitude(vec)
    DOUBLE PRECISION :: magnitude
    TYPE(vector3d) :: vec

    magnitude = SQRT(vec%x**2 + vec%y**2 + vec%z**2)
    
  END FUNCTION magnitude


  ! Calculating the dot product of a vector
  ! Note: DOT_PRODUCT is an intrinsic function, but
  ! it works on arrays, so I redefine it
  FUNCTION dot_prod(vec1, vec2)
    DOUBLE PRECISION :: dot_prod
    TYPE(vector3d) :: vec1, vec2

    dot_prod = vec1%x * vec2%x + &
               vec1%y * vec2%y + &
               vec1%z * vec2%z
    
  END FUNCTION dot_prod
  

  ! Distance between points
  FUNCTION distance(p1, p2)
    DOUBLE PRECISION :: distance
    TYPE(point3d) :: p1, p2

    distance = SQRT((p2%x - p1%x)**2 + &
                    (p2%y - p1%y)**2 + &
                    (p2%z - p1%z)**2)
    
  END FUNCTION distance

  
  ! Angle between vectors (in radians)
  FUNCTION angle(vec1, vec2)
    DOUBLE PRECISION :: angle
    TYPE(vector3d) :: vec1, vec2
    DOUBLE PRECISION :: mag1, mag2 ! For the vector magnitudes
    DOUBLE PRECISION :: dot        ! For the dot product
    
    mag1 = magnitude(vec1)
    mag2 = magnitude(vec2)

    dot = dot_prod(vec1, vec2)
    
    angle = ACOS(dot/(mag1*mag2))
    
  END FUNCTION angle

  
  ! Normalizing one vector
  FUNCTION normalize(vec)
    TYPE(vector3d) :: normalize
    TYPE(vector3d) :: vec
    DOUBLE PRECISION :: mag

    mag = magnitude(vec)

    normalize = vec / mag ! Should work because of divvr
    
  END FUNCTION normalize

  
  ! Cross product of two vectors
  FUNCTION cross_product(vec1, vec2)
    TYPE(vector3d) :: cross_product
    TYPE(vector3d) :: vec1, vec2

    cross_product%x = vec1%y * vec2%z - vec2%y * vec1%z
    cross_product%y = vec1%z * vec2%x - vec2%x * vec1%z
    cross_product%z = vec1%x * vec2%y - vec2%y * vec1%x
  END FUNCTION cross_product

  
  ! Returns an orthogonal vector to the input ones
  FUNCTION otrhv(vec1, vec2)
    TYPE(vector3d) :: orthv
    TYPE(vector3d) :: vec1, vec2

    ! We normalize it, otherwise it's the same as cross_product
    orthv = normalize(cross_product(vec1,vec2))
    
  END FUNCTION orthv

      
END MODULE geometry
