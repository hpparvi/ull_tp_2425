module particle
    use geometry  ! Use the geometry module defined earlier
    implicit none

    ! Define a type for a 3D particle
    type :: particle3d
        type(point3d) :: p      ! Position of the particle
        type(vector3d) :: v     ! Velocity of the particle
        real(8) :: m            ! Mass of the particle
    end type particle3d

end module particle
