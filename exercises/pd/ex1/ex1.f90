program ex1
    use particle  ! Import the particle module
    use geometry  ! Import the geometry module
    implicit none

    integer :: i, j
    integer :: n
    real(8) :: dt, t_end, t, dt_out, t_out
    real(8) ::  r3
    type(particle3d), dimension(:), allocatable :: particles  ! Array of particles
    type(vector3d) :: rji  ! Relative position vector
    type(vector3d) , dimension(:), allocatable :: a  ! Array of accelerations
    integer :: output_unit  ! Declare the output unit for file operations


    ! Read inputs
    print*, "Timestep value:"
    read*, dt
    print*, "Output time value:"
    read*, dt_out
    print*, "Value of the final time:"
    read*, t_end
    print*, "Number of particles:"
    read*, n

    ! Allocate the particles and accelerations arrays
    allocate(particles(n))
    allocate(a(n))

    ! Read mass, position (x,y,z), and velocity (vx,vy,vz) for each particle
    print*,"Write mass, position (x,y,z), and velocity (vx,vy,vz) for each particle in this order: m,x,y,z,vx,vy,vz:"
    do i = 1, n
        read(*,*) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
        particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
    end do


    ! Initialize acceleration to zero
    a = vector3d(0.0,0.0,0.0)

    do i = 1,n
        do j = i+1, n 
            rji = vector_between(particles(i)%p, particles(j)%p)  ! Calculate the vector connecting every pair of particles: rj - ri
            r3 = distance(particles(i)%p, particles(j)%p)**3  ! Calculate distance between two particles and raise it to the third power
            a(i) = vsum(a(i) , mulvr(particles(j)%m , divvr(rji , r3)))  !Calculate acceleration of particle i
            a(j) = vsub(mulvr(particles(i)%m , divvr(rji , r3)),a(j)) !Calculate acceleration of particle j
        end do

    end do


    ! Open output file for writing positions
    open(unit=output_unit, file='orbits.txt', status='replace')

    t_out = 0.0
    t = 0.0
    do while(t <= t_end)
        do i = 1, n
            particles(i)%v = vsum(particles(i)%v, mulvr(dt / 2, a(i)))
            particles(i)%p = sumvp(particles(i)%p , mulvr(dt, particles(i)%v))
        end do

        a = vector3d(0.0,0.0,0.0)
        do i = 1,n 
            do j=i+1,n 
                rji = vector_between(particles(i)%p, particles(j)%p)
                r3 = distance(particles(i)%p, particles(j)%p)**3
                a(i) = vsum(a(i) , mulvr(particles(j)%m , divvr(rji , r3)))
                a(j) = vsub(mulvr(particles(i)%m , divvr(rji , r3)),a(j))
            end do
        end do
        do i = 1, n
            particles(i)%v = vsum(particles(i)%v, mulvr(dt / 2, a(i)))
        end do

        t_out = t_out + dt

        if (t_out >= dt_out) then
            print*, "t =", t
            do i = 1,n 
                ! Print positions to file in column format: time, particle id, x, y, z
                WRITE(output_unit,*) t, i ,particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
            end do
            t_out = 0.0
        end if
        t = t + dt
    end do
    

end program ex1


