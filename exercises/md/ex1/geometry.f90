module geometry
    use iso_fortran_env
    
    implicit none
    

    type :: vector3d
        real(real64) :: x, y, z
    end type vector3d

    type :: point3d
        real(real64) :: x, y, z
    end type point3d

    interface operator(+)
        module procedure sumvp, sumpv, sumvv
    end interface operator(+)

    interface operator(-)
        module procedure subvp, subpv, subpp, subvv
    end interface operator(-)

    interface operator(*)
        module procedure mulrv, mulvr
    end interface operator(*)

    interface operator(/)
        module procedure divvr
    end interface operator(/)

    contains
        function sumvp(a, b)
            type(vector3d), intent(in) :: a
            type(point3d), intent(in) :: b
            type(point3d) :: sumvp
            sumvp = point3d(a%x + b%x, a%y + b%y, a%z + b%z)
        end function sumvp

        function sumpv(a, b)
            type(point3d), intent(in) :: a
            type(vector3d), intent(in) :: b
            type(point3d) :: sumpv
            sumpv = point3d(a%x + b%x, a%y + b%y, a%z + b%z)
        end function sumpv

        function sumvv(a, b)
            type(vector3d), intent(in) :: a, b
            type(vector3d) :: sumvv
            sumvv = vector3d(a%x + b%x, a%y + b%y, a%z + b%z)
        end function sumvv

        function subvp(a, b)
            type(vector3d), intent(in) :: a
            type(point3d), intent(in) :: b
            type(point3d) :: subvp
            subvp = point3d(a%x - b%x, a%y - b%y, a%z - b%z)
        end function subvp

        function subpv(a, b)
            type(point3d), intent(in) :: a
            type(vector3d), intent(in) :: b
            type(point3d) :: subpv
            subpv = point3d(a%x - b%x, a%y - b%y, a%z - b%z)
        end function subpv

        function subvv(a, b)
            type(vector3d), intent(in) :: a, b
            type(vector3d) :: subvv
            subvv = vector3d(a%x - b%x, a%y - b%y, a%z - b%z)
        end function subvv

        function subpp(a, b)
            type(point3d), intent(in) :: a, b
            type(vector3d) :: subpp
            subpp = vector3d(a%x - b%x, a%y - b%y, a%z - b%z)
        end function subpp

        function mulrv(r, v)
            real(real64), intent(in) :: r
            type(vector3d), intent(in) :: v
            type(vector3d) :: mulrv
            mulrv = vector3d(r * v%x, r * v%y, r * v%z)
        end function mulrv

        function mulvr(v, r)
            real(real64), intent(in) :: r
            type(vector3d), intent(in) :: v
            type(vector3d) :: mulvr
            mulvr = vector3d(v%x * r, v%y * r, v%z * r)
        end function mulvr

        function divvr(v, r)
            real(real64), intent(in) :: r
            type(vector3d), intent(in) :: v
            type(vector3d) :: divvr
            divvr = vector3d(v%x / r, v%y / r, v%z / r)
        end function divvr
        
        function distance(p1, p2)
            type(point3d), intent(in) :: p1, p2
            real(real64) :: distance
            distance = sqrt((p1%x - p2%x)**2 + (p1%y - p2%y)**2 + (p1%z - p2%z)**2)
        end function distance

        function normsquare(v)
            type(vector3d), intent(in) :: v
            real(real64) :: normsquare
            normsquare = v%x**2 + v%y**2 + v%z**2
        end function normsquare

        function angle(a, b)
            type(vector3d), intent(in) :: a, b
            real(real64) :: angle, ab
            ab = a%x * b%x + a%y * b%y + a%z * b%z
            angle = acos(ab / (normsquare(a)**2 * normsquare(b)**2))
        end function angle
        
        function normalize(v)
            type(vector3d), intent(in) :: v
            type(vector3d) :: normalize
            normalize = divvr(v, sqrt(v%x**2 + v%y**2 + v%z**2))
        end function normalize

        function cross_product(a, b)
            type(vector3d), intent(in) :: a, b
            type(vector3d) :: cross_product
            cross_product = vector3d(a%y * b%z - a%z * b%y, &
                                     a%z * b%x - a%x * b%z, &
                                     a%x * b%y - a%y * b%x)
        end function cross_product

        function orthv(a, b)
            type(vector3d), intent(in) :: a, b
            type(vector3d) :: orthv
            orthv = cross_product(a, b)
        end function orthv
end module geometry