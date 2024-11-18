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
      type(point3d) :: com  ! Center of mass of octant

      type(particle3d) :: particle  ! particle in octant
      integer :: idx  ! index of the particle in octant

      type(octant_pointer), dimension(2,2,2) :: suboctant

   end type octant

contains

   ! To be called ONCE for the root_octant and ALL particles
   ! assigns attributes %r_min, %r_max of the root_octant
   ! As we are dealing with Axis-Aligned Bounding Boxes (AABBs)
   ! this info is enough to get the cooridnates of all the edges

   subroutine give_root_octant_a_range(tronco, n_particles, particles)
      type(octant),pointer :: tronco
      integer :: n_particles
      type(particle3d) :: particles(n_particles)

      real, dimension(3) :: mins, maxs, span, mids

      ! get global minimum x,y,z within all particles
      mins(1) = minval(particles%p%x)
      mins(2) = minval(particles%p%y)
      mins(3) = minval(particles%p%z)

      ! get global maxmimu x,y,z within all particles
      maxs(1) = maxval(particles%p%x)
      maxs(2) = maxval(particles%p%y)
      maxs(3) = maxval(particles%p%z)

      ! Get span and add 1%
      span = maxval(maxs-mins)*1.01

      ! Get half span
      mids = (maxs + mins) / 2.0

      ! Assign info of edges
      tronco%range%r_min = mids - span/2.0
      tronco%range%r_max = mids + span/2.0

   end subroutine give_root_octant_a_range


   ! Nullification of suboctant pointers
   ! Important because:
   ! 1. Uninitialized pointers in suboctantscould lead to undefined behavior if accessed.
   ! 2. It makes it safe to check pointer association later using associated.
   subroutine nullify_suboctant_pointers(cell)
      type(octant), pointer :: cell
      integer :: i,j,k
      do i = 1,2
         do j = 1,2
            do k = 1,2
               nullify(cell%suboctant(i,j,k)%ptr)
            end do
         end do
      end do
   end subroutine nullify_suboctant_pointers


   ! Finds octant where to put a SINGLE particle (particl)
   ! We provide:
   ! :base_octant: to chech if particle belongs there
   ! :new_octant: may or may not be a new octant

   ! IF base_octant has 0 particles (category=0)
   ! THEN we assign the particle to that octant
   !      and we assign new_octant to base_octant

   ! ELIF base_octant has 1 particle (category=1)
   ! THEN we ALSO assign the particle to that octant!
   !
   !      Why? Because we allow maximum of 2 particles per octant
   !      And then when calling place_cell() later on
   !      then the place_cell() function
   !      will, if the octant has indeed explicitely 2 particles
   !      subdivide it in suboctants and reassign its 2 particles
   !
   !
   ! ELIF base_octant has 2 particles or more (category=2)
   ! THEN we loop through each sub-octant and call find_octant() recursively

   recursive subroutine find_octant(base_octant,new_octant,particl)

      ! Initial octant.
      type(octant),pointer :: base_octant

      ! To be the octant of the particle
      ! If 0 or 1 particle in it
      ! new_octant is base_octant
      type(octant),pointer :: new_octant

      ! comes into play if 2 or more particles in octant
      ! It is one of the possible eight sub-octants of base_octant
      ! then new_octant is asisgned to the sub-octant that contains the particle
      type(octant),pointer :: temp_octant

      type(particle3d) :: particl

      integer :: i,j,k

      select case (base_octant%category)
       case (2)  ! if it has more than 1 particle
         out: do i = 1,2
            do j = 1,2
               do k = 1,2
                  if (belongs(particl,base_octant%suboctant(i,j,k)%ptr)) then
                     call find_octant(base_octant%suboctant(i,j,k)%ptr,temp_octant,particl)
                     new_octant => temp_octant
                     exit out
                  end if
               end do
            end do
         end do out
       case default
         new_octant => base_octant
      end select
   end subroutine find_octant


   ! Pretty straightforward, outputs bool if particle is in octant or not
   function belongs(particl,goal)
      type(particle3d) :: particl
      type(octant), pointer :: goal
      logical :: belongs
      if (particl%p%x >= goal%range%r_min(1) .and. &
         particl%p%x <= goal%range%r_max(1) .and. &
         particl%p%y >= goal%range%r_min(2) .and. &
         particl%p%y <= goal%range%r_max(2) .and. &
         particl%p%z >= goal%range%r_min(3) .and. &
         particl%p%z <= goal%range%r_max(3)) then
         belongs = .true.
      else
         belongs = .false.
      end if
   end function belongs




   recursive subroutine place_particle_in_octant(octant_of_particle, particl, i)

      ! The current octant of the particle
      ! It can only have either:
      ! 1 particle: the one being called (particl)
      ! 2 particles (particl) and anothere one placed before
      type(octant),pointer :: octant_of_particle

      ! To be used for the case of 2 particles in octant_of_particle
      type(octant),pointer :: temp_octant

      type(particle3d) :: particl  ! the particle
      integer :: i  ! 1 the idx of particle in the array of particles_3d

      select case (octant_of_particle%category)

         ! Simplest case, pctant was empty, we fill it with particl
       case (0)
         octant_of_particle%category = 1
         octant_of_particle%particle = particl
         octant_of_particle%idx = i

         ! Octant already has a particle in it
       case (1)
         ! We need to subdivide it in 8 sub-octants
         ! This will place the particle that was already in octant_od_particle
         ! in one of the suboctants
         call create_suboctants(octant_of_particle)
         ! Find to which suboctant does the particle belong to
         call find_octant(octant_of_particle,temp_octant,particl)
         ! Place it there
         call place_particle_in_octant(temp_octant,particl,i)

         ! find_octant() prevents that we end up in this situation!
       case default
         print*,"should not be here. error!"
      end select

   end subroutine place_particle_in_octant


   ! This is called EXCLUSIVELY in place_particle_in_octant()
   ! ONLY when we were going to place a particle to an octant that already had
   ! exclusively one particle in it.

   ! So we are in a situation where we have an octant with an already assigned particle:
   ! the previous one
   ! not the one passed as argument in place_particle_in_octant()

   subroutine create_suboctants(base_octant)

      ! octant to subdivide
      type(octant), pointer :: base_octant

      ! particle that already was in that octant
      type(particle3d) :: particl

      ! integer coordinates of sub-octants (i,j,k) with i,j,k \in {1,2}
      integer, dimension(3) :: sub_octant_pos

      ! dummy indices
      integer :: i,j,k

      ! We get the particle that was already in that octant
      particl = base_octant%particle

      ! We set the category to 2 (two or more particles) because if we are calling subdivide is because we were going to place another particle to base_octant that already had one
      base_octant%category=2

      ! loop though each sub-octant
      do i = 1,2
         do j = 1,2
            do k = 1,2

               sub_octant_pos = (/i,j,k/)

               ! Important to avoid errors
               allocate(base_octant%suboctant(i,j,k)%ptr)

               ! Assigning r_min of suboctant i,j,k
               base_octant%suboctant(i,j,k)%ptr%range%r_min = &
                  compute_suboctant_ranges(0,base_octant,sub_octant_pos)

               ! Assigning r_max of suboctant i,j,k
               base_octant%suboctant(i,j,k)%ptr%range%r_max =&
                  compute_suboctant_ranges(1,base_octant,sub_octant_pos)

               ! Checking to what octant
               ! the particle that was already there
               ! belongs

               ! if belongs
               if (belongs(particl,base_octant%suboctant(i,j,k)%ptr)) then
                  ! assign particle to that sub-octant
                  base_octant%suboctant(i,j,k)%ptr%particle = particl
                  ! update category
                  base_octant%suboctant(i,j,k)%ptr%category = 1
                  ! provide the idx of particle
                  base_octant%suboctant(i,j,k)%ptr%idx = base_octant%idx

               else
                  ! category set to no particles
                  base_octant%suboctant(i,j,k)%ptr%category = 0
               end if

               ! important to avoid issues when allocating later on
               call nullify_suboctant_pointers(base_octant%suboctant(i,j,k)%ptr)

            end do
         end do
      end do
   end subroutine create_suboctants


   function compute_suboctant_ranges(min_or_max, parent_octant, octant_pos)&
      result(the_range)

      ! Values are either 0 or 1
      ! Are we computing the r_min edge (0) or r_max edge (1)
      integer :: min_or_max

      ! The parent octant
      type(octant), pointer :: parent_octant

      ! Binary position of sub-octant in parent_octant
      ! for which we are getting the range
      ! i.e. and array (i,j,k) with i,j \in {1,2}
      integer, dimension(3) :: octant_pos

      real, dimension(3)  :: half_side_of_parent

      ! output of the function
      real, dimension(3)  :: the_range

      ! dummy idx
      integer :: n

      half_side_of_parent = (parent_octant%range%r_min + parent_octant%range%r_max) / 2

      ! we initialize it like this and overwrite it next
      the_range = parent_octant%range%r_min

      select case (min_or_max)
         ! getting the r_min
       case (0)
         where (octant_pos == 1)
            the_range = parent_octant%range%r_min
         elsewhere
            the_range = half_side_of_parent
         endwhere

         ! getting the r_max
       case (1)
         where (octant_pos == 1)
            the_range = half_side_of_parent
         elsewhere
            the_range = parent_octant%range%r_max
         endwhere
      end select

   end function compute_suboctant_ranges

   ! pretty self explanatory
   recursive subroutine delete_empty_leaves(an_octant)
      type(octant),pointer :: an_octant
      integer :: i,j,k

      ! if any of suboctants is associated, all others are
      if (associated(an_octant%suboctant(1,1,1)%ptr)) then
         do i = 1,2
            do j = 1,2
               do k = 1,2
                  call delete_empty_leaves(an_octant%suboctant(i,j,k)%ptr)
                  if (an_octant%suboctant(i,j,k)%ptr%category == 0) then
                     deallocate (an_octant%suboctant(i,j,k)%ptr)
                  end if
               end do
            end do
         end do
      end if
   end subroutine delete_empty_leaves

   ! pretty self explanatory
   recursive subroutine delete_tree(an_octant)
      type(octant),pointer :: an_octant
      integer :: i,j,k
      do i = 1,2
         do j = 1,2
            do k = 1,2
               if (associated(an_octant%suboctant(i,j,k)%ptr)) then
                  call delete_tree(an_octant%suboctant(i,j,k)%ptr)
                  deallocate (an_octant%suboctant(i,j,k)%ptr)
               end if
            end do
         end do
      end do
   end subroutine delete_tree


!    recursive subroutine compute_masses(goal, n_particles, particles)
!       type(octant),pointer :: goal
!       integer :: n_particles
!       type(particle3d) :: particles(n_particles)
!       integer :: i,j,k
!       real :: mass

!       goal%mass = 0.

!       goal%com%x = 0.
!       goal%com%y = 0
!       goal%com%z = 0



!       select case (goal%category)
!        case (1)
!          goal%mass = particles(goal%idx)%m
!          goal%com%x = particles(goal%idx)%p%x
!          goal%com%y = particles(goal%idx)%p%y
!          goal%com%z = particles(goal%idx)%p%z
!        case (2)
!          do i = 1,2
!             do j = 1,2
!                do k = 1,2
!                   if (associated(goal%subcell(i,j,k)%ptr)) then
!                      call compute_masses(goal%subcell(i,j,k)%ptr, n_particles, particles)
!                      mass = goal%mass
!                      goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
!                      goal%com = (goal%com * mass + goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%com)
!                   end if
!                end do
!             end do
!          end do
!       end select
!    end subroutine compute_masses



!    subroutine compute_forces(head, n_particles,  particles, e, theta, a)
!       type(octant),pointer :: head
!       type(particle3d) :: particles(n_particles)
!       integer :: n_particles, i,j,k,start,end
!       real :: theta, e
!       type(vector3d), dimension(n_particles) :: a

!       do i = 1,n_particles
!          call compute_forces_aux(i,head, n_particles,  particles, e, theta, a)
!       end do
!    end subroutine compute_forces


!    recursive subroutine compute_forces_aux(goal, tree, n_particles, particles, e, theta, a)
!       integer :: goal, n_particles, i,j,k
!       type(octant), pointer :: tree
!       type(particle3d), dimension(n_particles) :: particles
!       type(vector3d), dimension(n_particles) :: a
!       type(vector3d) :: d_vec
!       real :: l,d, theta, e, dist

!       select case (tree%category)
!        case (1)
!          if (goal .ne. tree%idx) then
!             d_vec = tree%com - particles(goal)%p
!             dist = norm(d_vec)
!             dist = (dist**2 + e**2)**(1.5)
!             a(goal) = a(goal) + particles(tree%idx)%m * d_vec / dist
!          end if
!        case (2)

!          l = tree%range%r_max(1) - tree%range%r_min(1)
!          d_vec = tree%com - particles(goal)%p
!          dist = norm(d_vec)
!          if (l/d < theta) then
!             dist = (dist**2 + e**2)**(1.5)
!             a(goal) = a(goal) + tree%mass * d_vec / dist
!          else
!             do i = 1,2
!                do j = 1,2
!                   do k = 1,2
!                      if (associated(tree%subcell(i,j,k)%ptr)) then
!                         call compute_forces_aux(goal,tree%subcell(i,j,k)%ptr, n_particles, particles, e, theta,a)
!                      end if
!                   end do
!                end do
!             end do
!          end if
!       end select
!    end subroutine compute_forces_aux

end module octree
