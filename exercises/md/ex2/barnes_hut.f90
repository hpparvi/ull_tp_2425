module barnes_hut
    !$use omp_lib
    use iso_fortran_env
    use particle
    use geometry
    use calculations

    implicit none
    type(cell), pointer:: temp_cell

    contains 
    subroutine barnes_hut_tree(particles, head, n)
        type(cell), pointer, intent(inout):: head
        integer :: i
        integer, intent(in) :: n
        type(particle3d), intent(inout) :: particles(:)

        call calculate_ranges(particles, head)
        head%type = 0
        call nullify_pointers(head)
        do i = 1, n
            call find_cell(head, temp_cell, particles(i)%p)
            call place_cell(temp_cell, particles(i)%p, i)
        end do

        call borrar_empty_leaves(head)
    end subroutine barnes_hut_tree
end module barnes_hut

