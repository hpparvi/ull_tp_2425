module geometry
    use iso_fortran_env

    implicit none
    
    !Vector type
    type :: vector3d
        real(real64) :: x, y, z
    end type vector3d

    !Point type
    type :: point3d
        real(real64) :: x, y, z
    end type point3d

    !Opeartions between defined types
    interface operator(+)
        module procedure sumvp, sumpv, sumvv
    end interface operator(+)

    interface operator(-)
        module procedure subvp, subpv, subpp, subvv
    end interface operator(-)

    interface operator(*)
        module procedure mulrv, mulvr, mulrp, mulpr
    end interface operator(*)

    interface operator(/)
        module procedure divvr
    end interface operator(/)

    contains
        !vector3d + point3d = point3d
        function sumvp(a, b)
            type(vector3d), intent(in) :: a
            type(point3d), intent(in) :: b
            type(point3d) :: sumvp
            sumvp = point3d(a%x + b%x, a%y + b%y, a%z + b%z)
        end function sumvp
        
        !point3d + vector3d = point3d
        function sumpv(a, b)
            type(point3d), intent(in) :: a
            type(vector3d), intent(in) :: b
            type(point3d) :: sumpv
            sumpv = point3d(a%x + b%x, a%y + b%y, a%z + b%z)
        end function sumpv

        !vector3d + vector3d = vector3d
        function sumvv(a, b)
            type(vector3d), intent(in) :: a, b
            type(vector3d) :: sumvv
            sumvv = vector3d(a%x + b%x, a%y + b%y, a%z + b%z)
        end function sumvv

        !vector3d - point3d = point3d
        function subvp(a, b)
            type(vector3d), intent(in) :: a
            type(point3d), intent(in) :: b
            type(point3d) :: subvp
            subvp = point3d(a%x - b%x, a%y - b%y, a%z - b%z)
        end function subvp

        !point3d - vector3d = point3d
        function subpv(a, b)
            type(point3d), intent(in) :: a
            type(vector3d), intent(in) :: b
            type(point3d) :: subpv
            subpv = point3d(a%x - b%x, a%y - b%y, a%z - b%z)
        end function subpv

        !vector3d - vector3d = vector3d
        function subvv(a, b)
            type(vector3d), intent(in) :: a, b
            type(vector3d) :: subvv
            subvv = vector3d(a%x - b%x, a%y - b%y, a%z - b%z)
        end function subvv

        !point3d - point3d = vector3d
        function subpp(a, b)
            type(point3d), intent(in) :: a, b
            type(vector3d) :: subpp
            subpp = vector3d(a%x - b%x, a%y - b%y, a%z - b%z)
        end function subpp

        !r * vector3d = vector3d
        function mulrv(r, v)
            real(real64), intent(in) :: r
            type(vector3d), intent(in) :: v
            type(vector3d) :: mulrv
            mulrv = vector3d(r * v%x, r * v%y, r * v%z)
        end function mulrv

        !vector3d * r = vector3d
        function mulvr(v, r)
            real(real64), intent(in) :: r
            type(vector3d), intent(in) :: v
            type(vector3d) :: mulvr
            mulvr = vector3d(v%x * r, v%y * r, v%z * r)
        end function mulvr

        !r * point3d = point3d
        function mulrp(r, p)
            real(real64), intent(in) :: r
            type(point3d), intent(in) :: p
            type(point3d) :: mulrp
            mulrp = point3d(r * p%x, r * p%y, r * p%z)
        end function mulrp

        !point3d * r = point3d
        function mulpr(p, r)
            real(real64), intent(in) :: r
            type(point3d), intent(in) :: p
            type(point3d) :: mulpr
            mulpr = point3d(p%x * r, p%y * r, p%z * r)
        end function mulpr

        !vector3d / r = vector3d
        function divvr(v, r)
            real(real64), intent(in) :: r
            type(vector3d), intent(in) :: v
            type(vector3d) :: divvr
            divvr = vector3d(v%x / r, v%y / r, v%z / r)
        end function divvr
        
        !Distance between two points
        function distance(p1, p2)
            type(point3d), intent(in) :: p1, p2
            real(real64) :: distance
            distance = sqrt((p1%x - p2%x)**2 + (p1%y - p2%y)**2 + (p1%z - p2%z)**2)
        end function distance

        !Square of the norm of a vector
        function normsquare(v)
            type(vector3d), intent(in) :: v
            real(real64) :: normsquare
            normsquare = v%x**2 + v%y**2 + v%z**2
        end function normsquare

        !Angle between two vectors
        function angle(a, b)
            type(vector3d), intent(in) :: a, b
            real(real64) :: angle, ab
            ab = a%x * b%x + a%y * b%y + a%z * b%z
            angle = acos(ab / (normsquare(a)**2 * normsquare(b)**2))
        end function angle
        
        !Normalize a vector
        function normalize(v)
            type(vector3d), intent(in) :: v
            type(vector3d) :: normalize
            normalize = divvr(v, sqrt(v%x**2 + v%y**2 + v%z**2))
        end function normalize

        !Cross product between two vectors
        function cross_product(a, b)
            type(vector3d), intent(in) :: a, b
            type(vector3d) :: cross_product
            cross_product = vector3d(a%y * b%z - a%z * b%y, &
                                     a%z * b%x - a%x * b%z, &
                                     a%x * b%y - a%y * b%x)
        end function cross_product

        !Orthogonal vector to two vectors
        function orthv(a, b)
            type(vector3d), intent(in) :: a, b
            type(vector3d) :: orthv
            orthv = cross_product(a, b)
        end function orthv
end module geometry