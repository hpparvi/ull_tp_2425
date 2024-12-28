module particle
    !$use omp_lib
    use iso_fortran_env
    use geometry

    implicit none

    !Particle type: it has position, velocity, aceleration and mass
    type :: particle3d
        type(point3d) :: p
        type(vector3d) :: v
        type(vector3d) :: a
        real(real64) :: m
    end type particle3d
    

    contains
        !Reset acceleration of the particles to zero
        subroutine reset_a(particles)
            type(particle3d), intent(inout) :: particles(:)
            integer :: i
            !$omp do private(i)
            do i = 1, size(particles)
                particles(i)%a = vector3d(0.0, 0.0, 0.0)
            end do
            !$omp end do
        end subroutine reset_a

        !Update the position of the particles
        subroutine update_pos(particles, dt)
            type(particle3d), intent(inout) :: particles(:)
            real(real64), intent(in) :: dt
            integer :: i
            !$omp do private(i)
            do i = 1, size(particles)
                particles(i)%p = particles(i)%p + particles(i)%v * dt
            end do
            !$omp end do
        end subroutine update_pos

        !Update the velocity of the particles in a middle step
        subroutine update_vel(particles, dt)
            type(particle3d), intent(inout) :: particles(:)
            real(real64), intent(in) :: dt
            integer :: i
            !$omp do private(i)
            do i = 1, size(particles)
                particles(i)%v = particles(i)%v + particles(i)%a * dt * real(0.5, kind=real64)
            end do
            !$omp end do
        end subroutine update_vel

        !Update the acceleration of the particles
        subroutine update_a(particles)
            type(particle3d), intent(inout) :: particles(:)
            integer :: i, j
            real(real64) :: r2, r3
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