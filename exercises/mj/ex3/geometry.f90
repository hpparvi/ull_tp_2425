module geometry
    implicit none
    
    ! types vector3d and point3d
    type :: vector3d
        real(kind=8) :: x, y, z
    end type vector3d

    type :: point3d
        real(kind=8) :: x, y, z
    end type point3d

    ! function operators
    interface operator(+)
        module procedure sumvp, sumpvv, sumpv, sumpp
    end interface operator(+)

    interface operator(-)
        module procedure subvp, subvv, subpv, subpp
    end interface operator(-)

    interface operator(*)
        module procedure mulvr, mulrv, mulrp, mulpr, dot
    end interface operator(*)

    interface operator(/)
        module procedure divvr, divpr
    end interface operator(/)

contains
  
  FUNCTION point_to_vector(p) RESULT(v)
    TYPE(point3d), INTENT(IN) :: p
    TYPE(vector3d) :: v
    v%x = p%x
    v%y = p%y
    v%z = p%z
  END FUNCTION point_to_vector
  
  FUNCTION vector_to_point(v) RESULT(p)
    TYPE(vector3d), INTENT(IN) :: v
    TYPE(point3d) :: p
    p%x = v%x
    p%y = v%y
    p%z = v%z
  END FUNCTION vector_to_point


    ! operators matching
    function sumpv(p, v) result(r)
        type(point3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        type(point3d) :: r
        r%x = p%x + v%x
        r%y = p%y + v%y
        r%z = p%z + v%z
    end function sumpv

    function sumvp(v, p) result(r)
        type(point3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        type(point3d) :: r
        r%x = v%x + p%x
        r%y = v%y + p%y
        r%z = v%z + p%z
    end function sumvp

    function sumpvv(v1, v2) result(r)
        type(vector3d), intent(in) :: v1, v2
        type(vector3d) :: r
        r%x = v1%x + v2%x
        r%y = v1%y + v2%y
        r%z = v1%z + v2%z
    end function sumpvv

    function sumpp(p1, p2) result(r)
        type(point3d), intent(in) :: p1, p2
        type(vector3d) :: r
        r%x = p2%x + p1%x
        r%y = p2%y + p1%y
        r%z = p2%z + p1%z
    end function sumpp

    function subpv(p, v) result(r)
        type(point3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        type(point3d) :: r
        r%x = p%x - v%x
        r%y = p%y - v%y
        r%z = p%z - v%z
    end function subpv

function subvp(v, p) result(r)
        type(point3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        type(point3d) :: r
        r%x = v%x - p%x
        r%y = v%y - p%y
        r%z = v%z - p%z
    end function subvp

    function subvv(v1, v2) result(r)
        type(vector3d), intent(in) :: v1, v2
        type(vector3d) :: r
        r%x = v1%x - v2%x
        r%y = v1%y - v2%y
        r%z = v1%z - v2%z
    end function subvv

 function subpp(p1, p2) result(r)
        type(point3d), intent(in) :: p1, p2
        type(vector3d) :: r
        r%x = p2%x - p1%x
        r%y = p2%y - p1%y
        r%z = p2%z - p1%z
    end function subpp

    function mulvr(v, s) result(r)
        type(vector3d), intent(in) :: v
        real(kind=8), intent(in) :: s
        type(vector3d) :: r
        r%x = v%x * s
        r%y = v%y * s
        r%z = v%z * s
    end function mulvr

    function mulrv(s, v) result(r)
        type(vector3d), intent(in) :: v
        real(kind=8), intent(in) :: s
        type(vector3d) :: r
        r%x = s * v%x
        r%y = s * v%y
        r%z = s * v%z
    end function mulrv

    function mulrp(s, p) result(r)
        type(point3d), intent(in) :: p
        real(kind=8), intent(in) :: s
        type(point3d) :: r
        r%x = s * p%x
        r%y = s * p%y
        r%z = s * p%z
    end function mulrp

    function mulpr(p, s) result(r)
        type(point3d), intent(in) :: p
        real(kind=8), intent(in) :: s
        type(point3d) :: r
        r%x = p%x * s
        r%y = p%y * s
        r%z = p%z * s
    end function mulpr

    function dot(v1, v2) result(r)
        type(vector3d), intent(in) :: v1, v2
        real(kind=8) :: r
        r = v1%x * v2%x + v1%y * v2%y + v1%z * v2%z
    end function dot
    
    function divvr(v, s) result(r)
        type(vector3d), intent(in) :: v
        real(kind=8), intent(in) :: s
        type(vector3d) :: r
        r%x = v%x / s
        r%y = v%y / s
        r%z = v%z / s
    end function divvr

    function divpr(p, s) result(r)
        type(point3d), intent(in) :: p
        real(kind=8), intent(in) :: s
        type(vector3d) :: r
        r%x = p%x / s
        r%y = p%y / s
        r%z = p%z / s
    end function divpr

    ! function distance
    function distance(p1, p2) result(d)
        type(point3d), intent(in) :: p1, p2
        real(kind=8) :: d
        d = sqrt((p2%x - p1%x)**2 + (p2%y - p1%y)**2 + (p2%z - p1%z)**2)
    end function distance

    ! function angle
    function angle(a, b) result(theta)
        type(vector3d), intent(in) :: a, b
        real(kind=8) :: theta
        real(kind=8) :: dot_product
        real(kind=8) :: norm_a, norm_b

        dot_product = a%x * b%x + a%y * b%y + a%z * b%z
        norm_a = sqrt(a%x**2 + a%y**2 + a%z**2)
        norm_b = sqrt(b%x**2 + b%y**2 + b%z**2)
        theta = acos(dot_product / (norm_a * norm_b))
    end function angle

    ! function normalize
    function normalize(v) result(n)
        type(vector3d), intent(in) :: v
        type(vector3d) :: n
        real(kind=8) :: length

        length = sqrt(v%x**2 + v%y**2 + v%z**2)
        n = divvr(v, length)
    end function normalize

    ! function cross_product
    function cross_product(a, b) result(c)
        type(vector3d), intent(in) :: a, b
        type(vector3d) :: c
        c%x = a%y * b%z - a%z * b%y
        c%y = a%z * b%x - a%x * b%z
        c%z = a%x * b%y - a%y * b%x
    end function cross_product

    ! function orthv
    function orthv(a, b) result(o)
        type(vector3d), intent(in) :: a, b
        type(vector3d) :: o  
        o = cross_product(a, b)
    end function orthv

end module geometry
