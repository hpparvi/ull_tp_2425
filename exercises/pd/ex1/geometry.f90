module geometry
    implicit none

    ! Define a 3D vector type with 64-bit real components
    type :: vector3d
        real(8) :: x
        real(8) :: y
        real(8) :: z
    end type vector3d

    ! Define a 3D point type with 64-bit real components
    type :: point3d
        real(8) :: x
        real(8) :: y
        real(8) :: z
    end type point3d

    ! Operator overloading for vector + point
    interface operator(+)
        module procedure sumvp
    end interface operator(+)

    ! Operator overloading for point + vector
    interface operator(+)
        module procedure sumpv
    end interface operator(+)

    ! Operator overloading for vector - point
    interface operator(-)
        module procedure subvp
    end interface operator(-)

    ! Operator overloading for point - vector
    interface operator(-)
        module procedure subpv
    end interface operator(-)

    ! Operator overloading for vector * real
    interface operator(*)
        module procedure mulvr
    end interface operator(*)

    ! Operator overloading for real * vector
    interface operator(*)
        module procedure mulrv
    end interface operator(*)

    ! Operator overloading for vector / real
    interface operator(/)
        module procedure divvr
    end interface operator(/)

contains

    ! Function to normalize a vector
    function normalize(v) result(v_normalized)
        type(vector3d), intent(in) :: v
        type(vector3d) :: v_normalized
        real(8) :: magnitude

        ! Calculate the magnitude of the vector
        magnitude = sqrt(v%x**2 + v%y**2 + v%z**2)

        ! Check for zero magnitude to avoid division by zero
        if (magnitude > 0.0) then
            v_normalized%x = v%x / magnitude
            v_normalized%y = v%y / magnitude
            v_normalized%z = v%z / magnitude
        else
            print *, "Error: Cannot normalize a zero vector."
            v_normalized%x = 0.0
            v_normalized%y = 0.0
            v_normalized%z = 0.0
        end if
    end function normalize
    !Function to calculate de vector between two bodies
    function vector_between(p1, p2) result(rel_vec)
        type(point3d), intent(in) :: p1, p2
        type(vector3d) :: rel_vec

        rel_vec%x = p2%x - p1%x
        rel_vec%y = p2%y - p1%y
        rel_vec%z = p2%z - p1%z
    end function vector_between

    !Function to sum two vectors
    function vsum(v1, v2) result(sumv)
        type(vector3d), intent(in) :: v1, v2
        type(vector3d) :: sumv

        sumv%x = v1%x + v2%x
        sumv%y = v1%y + v2%y
        sumv%z = v1%z + v2%z
    end function vsum

    !Function to substract two vectors
    function vsub(v1, v2) result(subv)
        type(vector3d), intent(in) :: v1, v2
        type(vector3d) :: subv

        subv%x = v2%x - v1%x
        subv%y = v2%y - v1%y
        subv%z = v2%z - v1%z
    end function vsub

    ! Function to calculate the distance between two points
    function distance(p1, p2) result(dist)
        type(point3d), intent(in) :: p1, p2
        real(8) :: dist

        dist = sqrt((p2%x - p1%x)**2 + (p2%y - p1%y)**2 + (p2%z - p1%z)**2)
    end function distance

    ! Function to calculate the angle between two vectors in radians
    function angle(a, b) result(theta)
        type(vector3d), intent(in) :: a, b
        real(8) :: theta
        real(8) :: dot_product, mag_a, mag_b

        ! Calculate the dot product
        dot_product = a%x * b%x + a%y * b%y + a%z * b%z

        ! Calculate the magnitudes of the vectors
        mag_a = sqrt(a%x**2 + a%y**2 + a%z**2)
        mag_b = sqrt(b%x**2 + b%y**2 + b%z**2)

        ! Calculate the angle in radians using the arccos function
        if (mag_a > 0.0 .and. mag_b > 0.0) then
            theta = acos(dot_product / (mag_a * mag_b))
        else
            print *, "Error: One or both vectors have zero magnitude."
            theta = 0.0
        end if
    end function angle

    ! Function to calculate the cross product of two vectors
    function cross_product(a, b) result(c)
        type(vector3d), intent(in) :: a, b
        type(vector3d) :: c

        c%x = a%y * b%z - a%z * b%y  ! i-component
        c%y = a%z * b%x - a%x * b%z  ! j-component
        c%z = a%x * b%y - a%y * b%x  ! k-component
    end function cross_product

    ! Function to calculate a vector orthogonal to two given vectors
    function orthv(a, b) result(c)
        type(vector3d), intent(in) :: a, b
        type(vector3d) :: c

        c = cross_product(a, b)  ! Calculate the orthogonal vector using cross product
    end function orthv

    

    ! Function to add a vector to a point
    function sumvp(p, v) result(p_out)
        type(point3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        type(point3d) :: p_out

        p_out%x = p%x + v%x
        p_out%y = p%y + v%y
        p_out%z = p%z + v%z
    end function sumvp

    ! Function to add a point to a vector
    function sumpv(v, p) result(v_out)
        type(vector3d), intent(in) :: v
        type(point3d), intent(in) :: p
        type(vector3d) :: v_out

        v_out%x = v%x + p%x
        v_out%y = v%y + p%y
        v_out%z = v%z + p%z
    end function sumpv

    ! Function to subtract a vector from a point
    function subvp(p, v) result(p_out)
        type(point3d), intent(in) :: p
        type(vector3d), intent(in) :: v
        type(point3d) :: p_out

        p_out%x = p%x - v%x
        p_out%y = p%y - v%y
        p_out%z = p%z - v%z
    end function subvp

    ! Function to subtract a point from a vector
    function subpv(v, p) result(v_out)
        type(vector3d), intent(in) :: v
        type(point3d), intent(in) :: p
        type(vector3d) :: v_out

        v_out%x = v%x - p%x
        v_out%y = v%y - p%y
        v_out%z = v%z - p%z
    end function subpv

    ! Function to multiply a vector by a real number
    function mulrv(v, r) result(v_out)
        type(vector3d), intent(in) :: v
        real(8), intent(in) :: r
        type(vector3d) :: v_out

        v_out%x = v%x * r
        v_out%y = v%y * r
        v_out%z = v%z * r
    end function mulrv

    ! Function to multiply a real number by a vector
    function mulvr(r, v) result(v_out)
        real(8), intent(in) :: r
        type(vector3d), intent(in) :: v
        type(vector3d) :: v_out

        v_out%x = r * v%x
        v_out%y = r * v%y
        v_out%z = r * v%z
    end function mulvr

    ! Function to divide a vector by a real number
    function divvr(v, r) result(v_out)
        type(vector3d), intent(in) :: v
        real(8), intent(in) :: r
        type(vector3d) :: v_out

        if (r /= 0.0) then
            v_out%x = v%x / r
            v_out%y = v%y / r
            v_out%z = v%z / r
        else
            print *, "Error: Division by zero"
            v_out%x = 0.0
            v_out%y = 0.0
            v_out%z = 0.0
        end if
    end function divvr

end module geometry


