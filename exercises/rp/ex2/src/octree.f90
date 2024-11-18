module octree
   use geometry
   use particle
   implicit none

   type range
      real, dimension(3) :: r_min, r_max
   end type range

   type octant_pointer
      type(octant), pointer :: ptr
   end type octant_pointer

   type octant
      type(range) :: range
      integer :: category !! 0 = no bodies; 1 = one body; 2 = many bodies
      real :: mass  ! Total mass of octant
      real, dimension(3) :: com  ! Center of mass of octant

!! these two still dont know what are they used for
!   real, dimension(3) :: part
! integer :: pos

      type(octant_pointer), dimension(2,2,2) :: subcell

   end type octant

contains

   ! Given ALL particles, defines the coordinates of the box that contains all the particles
   subroutine get_cell_range(cell, n_particles, particles)
      type(octant),pointer :: cell
      integer :: n_particles
      type(particle3d) :: particles(n_particles)

      real, dimension(3) :: mins, maxs, span, mids

      mins(1) = minval(particles%p%x)
      mins(2) = minval(particles%p%y)
      mins(3) = minval(particles%p%z)

      maxs(1) = maxval(particles%p%x)
      maxs(2) = maxval(particles%p%y)
      maxs(3) = maxval(particles%p%z)

      ! shape (3) all elements equal to the max of maxs-mins
      span = maxval(maxs-mins)*1.01  
      mids = (maxs + mins) / 2.0

      cell%range%r_min = mids - span/2.0
      cell%range%r_max = mids + span/2.0

   end subroutine get_cell_range

   subroutine nullify_pointers(cell)
      type(octant), pointer :: cell
      integer :: i,j,k
      do i = 1,2
         do j = 1,2
            do k = 1,2
               nullify(cell%subcell(i,j,k)%ptr)
            end do
         end do
      end do
   end subroutine nullify_pointers
   
end module octree
