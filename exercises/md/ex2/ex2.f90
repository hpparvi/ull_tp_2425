program tree
    use iso_fortran_env
    use particle
    use geometry
    use calculations
    use utils

    implicit none

    integer :: i, n
    integer, parameter :: in = 1, out = 2
    real(real64) :: dt, t_end, t, dt_out, t_out
    type(particle3d), allocatable :: particles(:)
    character(len=100) :: sim_name, file_out
    !Input file of initial conditions
    character(len=100) :: file_in
    type(cell), pointer:: head, temp_cell

    !Check if an input file was provided or use terminal input
    if (command_argument_count() == 1) then
        !Get the input file name
        call get_command_argument(1, file_in)
        print*, 'Using input file: ', trim(file_in)

        !Read the input file in the format as the template in ics/
        call read_data(file_in, particles, sim_name, dt, dt_out, t_end, n, in)
    
    !If no input file was provided, use terminal input
    else if (command_argument_count() < 1) then
        print*, 'Using terminal input, printing positions, velocities and masses'
        call read_data_terminal(particles, sim_name, dt, dt_out, t_end, n)
    
        do i = 1, n
            print*, "Particle ", i, " position: x, y, z"
            read*, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
            print*, "Particle ", i, " velocity: x, y, z"
            read*, particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
            print*, "Particle ", i, " mass"
            read*, particles(i)%m
            print*, "____________________________________"
        end do
    !If more than one input file was provided, stop the program
    else
        print*, "Error: input file or terminal input not provided"
        print*, "Use for input file: ex1 <input_file>"
        print*, "Use for terminal input: ex1"
        stop
    end if

    !File to write the output, adjustl is used to put initial spaces at final, 
    !trim is used to remove the spaces at the end of the string
    file_out = 'output/' // trim(adjustl(sim_name)) // '.txt'
    
    print*, "Simulation name: ", sim_name
    print*, "Time step: ", dt
    print*, "Output time step: ", dt_out
    print*, "Final time: ", t_end
    print*, "Number of particles: ", n
    print*, "Output file: ", file_out
    print*, "____________________________________"

    print*, "Starting simulation..."

    call reset_a(particles)

    allocate(head)

    call calculate_ranges(particles, head)
    head%type = 0
    call nullify_pointers(head)

    do i = 1, n
        call find_cell(head, temp_cell, particles(i)%p)
        call place_cell(temp_cell, particles(i)%p, i)
    end do
    
    call borrar_empty_leaves(head)
    call calculate_masses(particles, head)
    call calculate_forces(particles, head)

    t_out = 0.0
    open(unit=out, file=file_out, action='write')
        write(out, "(A6, 1X, A12, 1X, A12, 1X, A12, 1X, A12, 1X, A12)") "#Index", "X", "Y", "Z", "Mass", "Time"
        do while (t < t_end)
            call update_vel(particles, dt)
            call update_pos(particles, dt)

            call borrar_tree(head)

            call calculate_ranges(particles, head)
            head%type = 0
            call nullify_pointers(head)

            do i = 1, n
                call find_cell(head, temp_cell, particles(i)%p)
                call place_cell(temp_cell, particles(i)%p, i)
            end do

            call reset_a(particles)

            call borrar_empty_leaves(head)
            call calculate_masses(particles, head)
            call calculate_forces(particles, head)
            call update_vel(particles, dt)

            t_out = t_out + dt
            if (t_out >= dt_out) then
                call write_data(particles, t, n, out)
                t_out = 0.0
            end if

            t = t + dt
        end do

    print*, "Done!"
end program tree
