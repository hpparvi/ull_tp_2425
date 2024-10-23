module particle
    use geometry
    implicit none
    type :: particle3d
        type(point3d) :: p
        type(vector3d) :: v
        real :: m
    end type particle3d
end module particle
