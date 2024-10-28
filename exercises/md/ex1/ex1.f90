program n_body
    use iso_fortran_env
    use particle
    use geometry

    implicit none

    integer :: i, n
    real(real64) :: dt, t_end, t, dt_out, t_out
    character(len=100) :: sim_name, file_out

    type(particle3d), allocatable :: particles(:)
    character(len=100) :: file_in


    if (command_argument_count() == 1) then
        call get_command_argument(1, file_in)
        print*, 'Using input file: ', trim(file_in)

        open(unit=1, file=file_in, action='read')
            read(1, *) 
            read(1, '(A)') sim_name
            read(1, *)
            read(1, *)
            read(1, *) dt
            read(1, *) 
            read(1, *)
            read(1, *) dt_out
            read(1, *) 
            read(1, *)
            read(1, *) t_end
            read(1, *) 
            read(1, *)
            read(1, *) n
            read(1, *)

            allocate(particles(n))

            read(1, *) 
            do i = 1, n
                read(1, *) particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
            end do
            read(1, *)
            
            read(1, *) 
            do i = 1, n
                read(1, *) particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
            end do
            read(1, *)

            read(1, *) 
            do i = 1, n
                read(1, *) particles(i)%m
            end do
        close(1)

    else if (command_argument_count() < 1) then
        print*, "No input file provided, using terminal input"
        print*, "Enter the simulation name:"
        read*, sim_name
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
    else
        print*, "Error: input file or terminal input not provided"
        print*, "Use for input file: ex1 <input_file>"
        print*, "Use for terminal input: ex1"
        stop
    end if

    file_out = 'output/' // trim(adjustl(sim_name)) // '.txt'
    
    print*, "Simulation name: ", sim_name
    print*, "Time step: ", dt
    print*, "Output time step: ", dt_out
    print*, "Final time: ", t_end
    print*, "Number of particles: ", n
    print*, "Output file: ", file_out
    print*, "____________________________________"

    do i = 1, n
        print*, "____________________________________"
        print*, "Particle ", i, " initial position: x, y, z"
        print*, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
        print*, "Particle ", i, " initial velocity: x, y, z"
        print*, particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
        print*, "Particle ", i, " mass"
        print*, particles(i)%m
    end do

    call reset_a(particles)
    call update_a(particles)

    t_out = 0.
    open(unit=2, file=file_out, action='write')
        write(2, "(A6, 1X, A12, 1X, A12, 1X, A12, 1X, A12)") "#Index", "X", "Y", "Z", "Mass"
        do while (t < t_end)
            call update_vel(particles, dt)
            call update_pos(particles, dt)
            call reset_a(particles)
            call update_a(particles)
            call update_vel(particles, dt)


            t_out = t_out + dt
            if  (t_out >= dt_out) then
                print*, "Positions at time ", t
                do i = 1, n
                    print*, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, particles(i)%m
                    write(2, "(I6, 1X, ES12.5, 1X, ES12.5, 1X, ES12.5, 1X, ES12.5)") i, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, particles(i)%m
                end do
                print*, "____________________________________"
                t_out = 0.
            end if

            t = t + dt
            
        end do
    close(2)
end program n_body