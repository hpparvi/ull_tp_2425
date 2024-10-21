module geometry
    implicit none
    type :: vector3d
       real :: x, y, z
    end type vector3d
    type :: point3d
       real :: x, y, z
    end type point3d
 
    interface operator (+)
       module procedure sumvp
    end interface
    interface operator (+)
       module procedure sumpv
    end interface
    interface operator (-)
       module procedure subvp
    end interface
    interface operator (-)
       module procedure subpv
    end interface
 
    interface operator (*)
       module procedure mulrv
    end interface
    interface operator (*)
       module procedure mulvr
    end interface
 
    interface operator (/)
       module procedure divvr
    end interface
 
 contains
 
    pure function sumvp(v, p) result(r)
       type(vector3d), intent(in) :: v
       type(point3d), intent(in) :: p
       type(point3d) :: r
 
       r%x = p%x + v%x
       r%y = p%y + v%y
       r%z = p%z + v%z
    end function sumvp
 
    pure function sumpv(p, v) result(r)
       type(vector3d), intent(in) :: v
       type(point3d), intent(in) :: p
       type(point3d) :: r
 
       r%x = p%x + v%x
       r%y = p%y + v%y
       r%z = p%z + v%z
    end function sumpv
 
    pure function subvp(v, p) result(r)
       type(vector3d), intent(in) :: v
       type(point3d), intent(in) :: p
       type(point3d) :: r
 
       r%x = v%x - p%x
       r%y = v%y - p%y
       r%z = v%z - p%z
    end function subvp
 
    pure function subpv(p, v) result(r)
       type(vector3d), intent(in) :: v
       type(point3d), intent(in) :: p
       type(point3d) :: r
 
       r%x = p%x - v%x
       r%y = p%y - v%y
       r%z = p%z - v%z
    end function subpv
 
    pure function mulrv(s, v) result(r)
       type(vector3d), intent(in) :: v
       real, intent(in) :: s
       type(vector3d) :: r
 
       r%x = s * v%x
       r%y = s * v%y
       r%z = s * v%z
    end function mulrv
 
    pure function mulvr(v, s) result(r)
       type(vector3d), intent(in) :: v
       real, intent(in) :: s
       type(vector3d) :: r
 
       r%x = s * v%x
       r%y = s * v%y
       r%z = s * v%z
    end function mulvr
 
    pure function divvr(v, s) result(r)
       type(vector3d), intent(in) :: v
       real, intent(in) :: s
       type(vector3d) :: r
 
       r%x = v%x / s
       r%y = v%y / s
       r%z = v%z / s
 
    end function divvr
 
    pure function distance(p1, p2) result(r)
       type(point3d), intent(in) :: p1
       type(point3d), intent(in) :: p2
       real :: r
 
       r = sqrt((p1%x - p2%x)**2+(p1%y - p2%y)**2+(p1%z - p2%z)**2)
 
    end function distance
 
    pure function dot(v1, v2) result(r)
       type(vector3d), intent(in) :: v1
       type(vector3d), intent(in) :: v2
       real :: r
       r = (v1%x)*(v2%x) + (v1%y)*(v2%y) + (v1%z)*(v2%z) 
    end function dot
 
    pure function norm(v) result(r)
    type(vector3d), intent(in) :: v
    real :: r
    r = sqrt(dot(v, v))
    end function norm
 
    pure function angle(v1, v2) result(theta)
       type(vector3d), intent(in) :: v1
       type(vector3d), intent(in) :: v2
       real :: theta
       theta = acos(dot(v1, v2) / (norm(v1) * norm(v2)))
    end function angle
 
    pure function normalize(v) result(v_norm)
    type(vector3d), intent(in) :: v
    type(vector3d) :: v_norm
    v_norm = v / norm(v)
    end function normalize
 
    pure function cross_product(v1, v2) result(v_cross)
       type(vector3d), intent(in) :: v1
       type(vector3d), intent(in) :: v2
       type(vector3d) :: v_cross
       v_cross%x = (v1%y)*(v2%z) - (v1%z)*(v2%y)
       v_cross%y = (v1%z)*(v2%x) - (v1%x)*(v2%z)
       v_cross%z = (v1%x)*(v2%y) - (v1%y)*(v2%x)
    end function cross_product
 
    pure function orthv(v1, v2) result(v_orth)
       type(vector3d), intent(in) :: v1
       type(vector3d), intent(in) :: v2
       type(vector3d) :: v_orth
       v_orth = cross_product(v1, v2)
       v_orth = v_orth / norm(v_orth)
    end function orthv
 
 end module geometry