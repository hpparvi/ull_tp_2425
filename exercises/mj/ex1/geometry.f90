module geometry
    implicit none

    ! types vector3d and point3d
    type :: vector3d
        real(kind=8) :: x, y, z
    end type vector3d

    type :: point3d
        real(kind=8) :: x, y, z
    end type point3d

    ! operator functions
    interface operator(+)
        module procedure sumvp, sumpv
    end interface operator(+)

    interface operator(-)
        module procedure subvp
    end interface operator(-)

    interface operator(*)
        module procedure mulvr
    end interface operator(*)

    interface operator(/)
        module procedure divvr
    end interface operator(/)

contains

    ! operators matching
    function sumvp(p, v) result(r)
        type(point3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        type(point3d) :: r
        r%x = p%x + v%x
        r%y = p%y + v%y
        r%z = p%z + v%z
    end function sumvp

    function sumpv(p, v) result(r)
        type(vector3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        type(vector3d) :: r
        r%x = p%x + v%x
        r%y = p%y + v%y
        r%z = p%z + v%z
    end function sumpv

    function subvp(p, v) result(r)
        type(point3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        type(point3d) :: r
        r%x = p%x - v%x
        r%y = p%y - v%y
        r%z = p%z - v%z
    end function subvp

    function mulvr(v, s) result(r)
        type(vector3d), intent(in) :: v
        real(kind=8), intent(in) :: s
        type(vector3d) :: r
        r%x = v%x * s
        r%y = v%y * s
        r%z = v%z * s
    end function mulvr

    function divvr(v, s) result(r)
        type(vector3d), intent(in) :: v
        real(kind=8), intent(in) :: s
        type(vector3d) :: r
        r%x = v%x / s
        r%y = v%y / s
        r%z = v%z / s
    end function divvr

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

    ! function cross_profuct
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
