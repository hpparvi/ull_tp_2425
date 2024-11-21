program tree
    use iso_fortran_env
    use particle
    use geometry
    use calculations

    implicit none

    integer :: i, n
    real(real64) :: dt, t_end, t, dt_out, t_out
    type(particle3d), allocatable :: particles(:)

    type(cell), pointer:: head, temp_cell

    read*, dt
    read*, dt_out
    read*, t_end
    read*, n

    allocate(particles(n))

    do i = 1, n
        read*, particles(i)%m, particles(i)%p, particles(i)%v
    end do

    call reset_a(particles)

    allocate(head)

    call calculate_range(head)
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

        call borrar_empty_leaves(head)
        call calculate_masses(particles, head)
        call calculate_forces(particles, head)
        call update_vel(particles, dt)

        t_out = t_out + dt
        if (t_out >= dt_out) then
            do i = 1, 10
                print*, particles(i)%p
            end do
            print*, "--------------------------------"
            print*, ""
            t_out = 0.0
        end if
    end do
end program tree
