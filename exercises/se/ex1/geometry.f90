! This module contains the definitions of several
! useful operators

MODULE geometry
  IMPLICIT NONE

  ! Define vectors
  TYPE vector3d
     REAL :: x, y, z
  END TYPE vector3d

  ! Define points
  TYPE point3d
     REAL :: x, y, z
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
    REAL           :: re
    
    mulrv%x = re * vec%x
    mulrv%y = re * vec%y
    mulrv%z = re * vec%z
  END FUNCTION mulrv

  
  ! Multiplying a vector and a real
  FUNCTION mulvr(vec, re)
    TYPE(vector3d) :: mulvr
    TYPE(vector3d) :: vec
    REAL           :: re
    
    mulvr%x = re * vec%x
    mulvr%y = re * vec%y
    mulvr%z = re * vec%z
  END FUNCTION mulvr

  
  ! Dividing a vector by a real
  FUNCTION divvr
    TYPE(vector3d) :: divvr
    TYPE(vector3d) :: vec
    REAL           :: re
    
    divvr%x = vec%x / re
    divvr%y = vec%y / re
    divvr%z = vec%z / re
  END FUNCTION divvr



END MODULE geometry
