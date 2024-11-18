module calculations
    use iso_fortran_env
    use geometry
    use particle

    implicit none
    real(real64), parameter :: theta = 1.0

    type range 
    real(real64), dimension(3) :: min, max
    end type range

    type CPtr
        type(cell), pointer :: ptr
    end type CPtr

    type cell
        type(range) :: range
        type(point3d) :: part
        integer :: pos
        integer :: type ! 0 = no particle; 1 = particle; 2 = conglomerado
        real(real64) :: mass
        type(point3d) :: c_o_m
        type(CPtr), dimension(2,2,2) :: subcell
    end type cell


    contains
        subroutine calculate_ranges(particles, goal)
            type(cell), pointer :: goal
            type(particle3d), intent(in) :: particles(:)
            
            
            real(real64), dimension(3) :: mins, maxs, medios
            real(real64) :: span

            mins(1) = minval(particles%p%x)
            maxs(1) = maxval(particles%p%x)
            mins(2) = minval(particles%p%y)
            maxs(2) = maxval(particles%p%y)
            mins(3) = minval(particles%p%z)
            maxs(3) = maxval(particles%p%z)

            span = maxval(maxs - mins) * 1.1
            medios = (mins + maxs) * 0.5
            goal%range%min = medios - span * 0.5
            goal%range%max = medios + span * 0.5
        end subroutine calculate_ranges

        recursive subroutine find_cell(root, goal, part)
            type(point3d) :: part
            type(cell), pointer :: root, goal, temp
            integer :: i, j, k

            select case (root%type)
            case (2)
                out: do i = 1,2
                    do j = 1,2
                        do k = 1,2
                            if (belongs(part, root%subcell(i,j,k)%ptr)) then
                                call find_cell(root%subcell(i,j,k)%ptr, temp, part)
                                goal => temp
                                exit out
                            end if
                        end do
                    end do
                end do out
            case default
                goal => root
            end select
        end subroutine find_cell

        recursive subroutine place_cell(goal, part, n)
            type(cell), pointer :: goal, temp
            type(point3d) :: part
            integer :: n

            select case (goal%type)
            case (0)
                goal%type = 1
                goal%part = part
                goal%pos = n
            case (1)
                call crear_subcells(goal)
                call find_cell(goal, temp, part)
                call place_cell(temp, part, n)
            case default
                print*, "SHOULD NOT BE HERE, ERROR"
            end select
        end subroutine place_cell

        subroutine crear_subcells(goal)
            type(cell), pointer :: goal
            type(point3d) :: part
            integer :: i, j, k
            integer, dimension(3) :: octant

            part = goal%part
            goal%type = 2

            do i = 1,2
                do j = 1,2
                    do k = 1,2 
                        octant = (/i,j,k/)
                        allocate(goal%subcell(i,j,k)%ptr)
                        goal%subcell(i,j,k)%ptr%range%min = calcular_range(0, goal, octant)
                        goal%subcell(i,j,k)%ptr%range%max = calcular_range(1, goal, octant)

                        if (belongs(part, goal%subcell(i,j,k)%ptr)) then
                            goal%subcell(i,j,k)%ptr%part = part
                            goal%subcell(i,j,k)%ptr%type = 1
                            goal%subcell(i,j,k)%ptr%pos = goal%pos
                        else
                            goal%subcell(i,j,k)%ptr%type = 0
                        end if
                        call nullify_pointers(goal%subcell(i,j,k)%ptr)
                    end do
                end do
            end do
        end subroutine crear_subcells

        subroutine nullify_pointers(goal)
            type(cell), pointer :: goal
            integer :: i,j,k

            do i = 1,2
                do j = 1,2
                    do k = 1,2
                        nullify(goal%subcell(i,j,k)%ptr)
                    end do
                end do
            end do
        end subroutine nullify_pointers

        function belongs(part, goal)
            type(point3d) :: part
            type(cell), pointer :: goal
            logical :: belongs

            if (part%x >= goal%range%min(1) .and. part%x <= goal%range%max(1) &
                .and. part%y >= goal%range%min(2) .and. part%y <= goal%range%max(2) &
                .and. part%z >= goal%range%min(3) .and. part%z <= goal%range%max(3)) then
                belongs = .true.
            else
                belongs = .false.
            end if
        end function belongs

        function calcular_range(what, goal, octant)
            integer :: what
            type(cell), pointer :: goal
            integer, dimension(3) :: octant
            real(real64), dimension(3) :: calcular_range, valor_medio
            
            valor_medio = (goal%range%min + goal%range%max) * 0.5
            
            select case (what)
                case (0)
                    where (octant == 1)
                        calcular_range = goal%range%min
                    elsewhere
                        calcular_range = valor_medio
                    end where
                case (1)
                    where (octant == 1)
                        calcular_range = valor_medio
                    elsewhere
                        calcular_range = goal%range%max
                    end where
            end select 
        end function calcular_range

        recursive subroutine borrar_empty_leaves(goal)
            type(cell), pointer :: goal
            integer :: i, j, k

            if (associated(goal%subcell(1,1,1)%ptr)) then
                do i = 1,2
                    do j = 1,2
                        do k = 1,2
                            call borrar_empty_leaves(goal%subcell(i,j,k)%ptr)
                            if (goal%subcell(i,j,k)%ptr%type == 0) then
                                deallocate(goal%subcell(i,j,k)%ptr)
                            end if
                        end do
                    end do
                end do
            end if
        end subroutine borrar_empty_leaves

        recursive subroutine borrar_tree(goal)
            type(cell), pointer :: goal
            integer :: i, j, k

            do i = 1,2
                do j = 1,2
                    do k = 1,2
                        if (associated(goal%subcell(i,j,k)%ptr)) then
                            call borrar_tree(goal%subcell(i,j,k)%ptr)
                            deallocate(goal%subcell(i,j,k)%ptr)
                        end if
                    end do
                end do
            end do
        end subroutine borrar_tree

        recursive subroutine calculate_masses(particles, goal)
            type(particle3d), intent(in) :: particles(:)
            type(cell), pointer :: goal
            integer :: i, j, k
            real(real64) :: mass

            goal%mass = 0.0
            goal%c_o_m%x = 0.0
            goal%c_o_m%y = 0.0
            goal%c_o_m%z = 0.0

            select case (goal%type)
                case(1)
                    goal%mass = particles(goal%pos)%m
                    goal%c_o_m = particles(goal%pos)%p
                case(2)
                    do i = 1,2
                        do j = 1,2
                            do k = 1,2
                                if (associated(goal%subcell(i,j,k)%ptr)) then
                                    call calculate_masses(particles, goal%subcell(i,j,k)%ptr)
                                    mass = goal%mass
                                    goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
                                    goal%c_o_m%x = (mass * goal%c_o_m%x + goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m%x) / goal%mass
                                    goal%c_o_m%y = (mass * goal%c_o_m%y + goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m%y) / goal%mass
                                    goal%c_o_m%z = (mass * goal%c_o_m%z + goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m%z) / goal%mass
                                end if
                            end do
                        end do
                    end do
            end select
        end subroutine calculate_masses

        subroutine calculate_forces(particles, head)
            type(particle3d), intent(inout) :: particles(:)
            type(cell), pointer :: head
            integer :: i, n
            
            n = size(particles)
            do i = 1, n
                call calculate_forces_aux(particles, i, head)
            end do
        end subroutine calculate_forces

        recursive subroutine calculate_forces_aux(particles, goal, tree)
            type(particle3d), intent(inout) :: particles(:)
            type(cell), pointer :: tree
            integer :: i, j, k, goal
            real(real64) :: l, D, r2, r3
            type(vector3d) :: rji

            select case (tree%type)
            case(1)
                if (goal .ne. tree%pos) then
                    rji = tree%c_o_m - particles(goal)%p
                    r2 = normsquare(rji)
                    r3 = r2*sqrt(r2)
                    particles(goal)%a = particles(goal)%a + particles(tree%pos)%m * rji / r3
                end if
            case(2)
                l = tree%range%max(1) - tree%range%min(1)
                rji = tree%c_o_m - particles(goal)%p
                r2 = normsquare(rji)
                D = sqrt(r2)
                if (1/D < theta) then
                    r3 = r2 * D
                    particles(goal)%a = particles(goal)%a + tree%mass * rji / r3
                else
                    do i = 1,2
                        do j = 1,2
                            do k = 1,2
                                if (associated(tree%subcell(i,j,k)%ptr)) then
                                    call calculate_forces_aux(particles, goal, tree%subcell(i,j,k)%ptr)
                                end if
                            end do
                        end do
                    end do
                end if
            end select
        end subroutine calculate_forces_aux
end module calculations