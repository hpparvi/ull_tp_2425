module barnes_hut
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

        !First, we need to calculate the ranges of the particles
        call calculate_ranges(particles, head)
        head%type = 0
        call nullify_pointers(head)
        !For each particle, we need to place it in the tree
        do i = 1, n
            call find_cell(head, temp_cell, particles(i)%p)
            call place_cell(temp_cell, particles(i)%p, i)
        end do
        !Delete empty leaves
        call borrar_empty_leaves(head)
    end subroutine barnes_hut_tree
end module barnes_hut

