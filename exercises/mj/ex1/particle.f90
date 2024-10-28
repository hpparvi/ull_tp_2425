module particle
    use geometry 
    implicit none

    ! particle3d
    type :: particle3d
        type(point3d) :: p   ! position
        type(vector3d) :: v   ! velocity
        real(kind=8) :: m     ! mass
    end type particle3d

end module particle
