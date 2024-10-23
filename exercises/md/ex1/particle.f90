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

        subroutine update_a(particles)
            type(particle3d), intent(inout) :: particles(:)
            integer :: i, j
            real :: r2, r3
            type(vector3d):: rji
            
            do i = 1, size(particles)
                do j = i+1, size(particles)
                    rji = particles(j)%p - particles(i)%p
                    r2 = normsquare(rji)
                    r3 = r2 * sqrt(r2)
                    particles(i)%a = particles(i)%a + particles(j)%m * rji / r3
                    particles(j)%a = particles(j)%a - particles(i)%m * rji / r3
                end do
            end do
        end subroutine update_a
end module particle