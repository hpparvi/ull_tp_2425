module particle
    use geometry
    implicit none

    type :: particle3d
        type(point3d) :: p
        type(vector3d) :: v
        type(vector3d) :: a
        real :: m
    end type particle3d

    contains
        subroutine reset_a(particles)
            type(particle3d), intent(inout) :: particles(:)
            integer :: i
            do i = 1, size(particles)
                particles(i)%a = vector3d(0.0, 0.0, 0.0)
            end do
        end subroutine reset_a

        subroutine update_pos(particles, dt)
            type(particle3d), intent(inout) :: particles(:)
            real, intent(in) :: dt
            integer :: i
            do i = 1, size(particles)
                particles(i)%p = particles(i)%p + particles(i)%v * dt
            end do
        end subroutine update_pos

        subroutine update_vel(particles, dt)
            type(particle3d), intent(inout) :: particles(:)
            real, intent(in) :: dt
            integer :: i
            do i = 1, size(particles)
                particles(i)%v = particles(i)%v + particles(i)%a * dt * 0.5
            end do
        end subroutine update_vel
end module particle