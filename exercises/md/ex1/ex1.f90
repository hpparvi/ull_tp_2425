program leapfrog
    use geometry
    use particle

    implicit none

    integer :: i, j
    integer :: n
    real :: dt, t_end, t, dt_out, t_out
    real :: r2, r3

    type(particle3d), allocatable :: particles(:)
    type(vector3d) :: rji

    print*, "Enter the time step:"
    read*, dt
    print*, "Enter the output time step:"
    read*, dt_out
    print*, "Enter the final time:"
    read*, t_end
    print*, "Enter the number of particles:"
    read*, n

    allocate(particles(n))

    do i = 1, n
        print*, "Particle ", i, " position: x, y, z"
        read*, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
        print*, "Particle ", i, " velocity: x, y, z"
        read*, particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
        print*, "Particle ", i, " mass"
        read*, particles(i)%m
        print*, "____________________________________"
    end do

    call reset_a(particles)

    do i = 1, n
        do j = i+1, n 
            rji = particles(j)%p - particles(i)%p
            r2 = normsquare(rji)
            r3 = r2 * sqrt(r2)
            particles(i)%a = particles(i)%a + particles(j)%m * rji / r3
            particles(j)%a = particles(j)%a - particles(i)%m * rji / r3
        end do
    end do

    t_out = 0.
    open(unit=1, file='output/positions.txt', action='write')
        do while (t < t_end)
            call update_pos(particles, dt)
            call update_vel(particles, dt)
            call reset_a(particles)

            do i = 1, n
                do j = i+1, n 
                    rji = particles(j)%p - particles(i)%p
                    r2 = normsquare(rji)
                    r3 = r2 * sqrt(r2)
                    particles(i)%a = particles(i)%a + particles(j)%m * rji / r3
                    particles(j)%a = particles(j)%a - particles(i)%m * rji / r3
                end do
            end do

            call update_vel(particles, dt)
            t_out = t_out + dt
            if  (t_out >= dt_out) then
                print*, "Positions at time ", t
                    do i = 1, n
                        print*, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
                        write(1, "(I5, 1X, ES12.5, 1X, ES12.5, 1X, ES12.5)") i, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
                    end do
                print*, "____________________________________"
                t_out = 0.
            end if

            t = t + dt
        end do
    close(1)
end program leapfrog