program tree
   use geometry
   use particle
   use i_o_utils
   use octree
   implicit none

   character(len=500) :: ics_file
   character(len=500) :: savefolder

   real :: G ! Gravitational constant
   integer :: N ! number of particles
   real :: T ! Total integration time
   real :: dt ! Time-step
   real :: e ! Smoothing length
   real :: theta ! octree threshold
   type(particle3d), allocatable :: bodies(:)

   integer :: N_time_steps ! Number of time-steps
   integer :: N_snapshots  ! Number of snapshots to save
   integer :: snap_step  ! factor between two above

   type(octant), pointer :: root_octant, temp_octant

   type(vector3d), allocatable :: a(:) ! accelerations

   integer :: i, j, k, ii  ! dummy loop varables

   ! to get the time
   real :: start_time, end_time, elapsed_time

   ! progress marker stuff
   integer :: progress_marker, progress_marker_aux
   progress_marker = 1  ! every 1%
   progress_marker_aux = progress_marker


   call read_config(ics_file, savefolder, N_snapshots, G, T, dt, e, theta)
   N_time_steps = int(T/dt)

   call get_ics_from_file(ics_file, N, bodies)
   call create_snapshot_folder(savefolder)
   call save_data(0, savefolder, bodies, N)  ! Save positions at time=0

   call get_N_snapshots(N_snapshots, N_time_steps)
   snap_step = N_time_steps/N_snapshots

   print*, ' '
   print*, "Number of particles:", N
   print*, 'Total integration time:', T
   print*, 'Time step:', dt
   print*, 'Number of steps:', N_time_steps
   print*, 'Number of snapshots:', N_snapshots


   allocate(root_octant)
   call give_root_octant_a_range(root_octant, N, bodies)
   root_octant%category = 0  ! No bodies in `root_octant`, yet
   call nullify_suboctant_pointers(root_octant)  ! setting subcell pointers to null


   do i = 1,N
      call find_octant(root_octant,temp_octant,bodies(i))
      call place_particle_in_octant(temp_octant, bodies(i),i)
   end do

   call delete_empty_leaves(root_octant)

   call compute_masses(root_octant) !, N, bodies)

   allocate(a(N))
   ! Initialize accelerations to zero
   do i = 1, N
      a(i)%x = 0.0
      a(i)%y = 0.0
      a(i)%z = 0.0
   end do

   do i = 1, N
      call compute_forces(i, root_octant)
   end do

   ! save ICS
   call save_data(0, savefolder, bodies, N)


   call cpu_time(start_time)

   j = 1
   k = 1
   do i = 1, N_time_steps

      ! First part of leapfrog step
      do ii=1, N
         bodies(ii)%v = bodies(ii)%v + (a(ii) * dt) * 0.5
         bodies(ii)%p = bodies(ii)%p + bodies(ii)%v * dt
      end do

      ! RECOMPUTE ALL ACCELERATIONS
      call delete_tree(root_octant)
      call give_root_octant_a_range(root_octant, N, bodies)
      root_octant%category = 0
      call nullify_suboctant_pointers(root_octant)
      do ii = 1,N
         call find_octant(root_octant,temp_octant,bodies(ii))
         call place_particle_in_octant(temp_octant, bodies(ii),ii)
      end do
      call delete_empty_leaves(root_octant)
      call compute_masses(root_octant)
      do ii = 1, N
         a(ii)%x = 0.0
         a(ii)%y = 0.0
         a(ii)%z = 0.0
      end do
      do ii = 1, N
         call compute_forces(ii, root_octant)
      end do

      ! Second part of leapfrog step
      do ii=1, N
         bodies(ii)%v = bodies(ii)%v + (a(ii) * dt) * 0.5
      end do

      ! Save data
      if (MOD(i-1, snap_step) == 0) then
         call save_data(j, savefolder, bodies, N)
         j=j+1
      end if

      if (int(real(i) / N_time_steps * 100) >= progress_marker) then
         write(*, '(I3, A)', ADVANCE="NO") progress_marker, '%'  ! Print progress without advancing
         write(*, '(A)', ADVANCE="NO") CHAR(13)                 ! Carriage return to overwrite
         progress_marker = progress_marker + progress_marker_aux
      end if

   end do

   call cpu_time(end_time)

   elapsed_time = end_time - start_time
   print *, "Done! Elapsed CPU time:", elapsed_time, 's'


contains
   recursive subroutine compute_forces(idx, base_octant)

      ! idx of particle
      integer, intent(in) :: idx

      type(octant), pointer, intent(in) :: base_octant

      type(vector3d) :: d_vec ! Relatve position between particles
      real :: d, denominator ! Relatve distance between particles
      real :: octant_span

      integer :: iii, jjj, kkk

      select case (base_octant%category)

         ! Only one particle in octant
       case (1)
         ! Avoid self-interaction
         if (idx .ne. base_octant%idx) then
            d_vec = bodies(idx)%p - bodies(base_octant%idx)%p
            d = norm(d_vec)
            denominator = (d**2 + e**2)**(1.5)
            a(idx) = a(idx) - G * bodies(base_octant%idx)%m * d_vec / denominator
         end if

         ! More than one particle in octant
       case (2)
         ! Span is the same in all three dimensions
         octant_span = base_octant%range%r_max(1) - base_octant%range%r_min(1)

         ! Distance to COM of octant
         d_vec = bodies(idx)%p - base_octant%com
         d = norm(d_vec)

         ! Apply opening angle criterion
         if (octant_span / d < theta) then
            denominator = (d**2 + e**2)**(1.5)
            a(idx) = a(idx) - G * base_octant%mass * d_vec / denominator
         else
            ! Recurse into suboctants
            do iii = 1, 2
               do jjj = 1, 2
                  do kkk = 1, 2
                     if (associated(base_octant%suboctant(iii, jjj, kkk)%ptr)) then
                        call compute_forces(idx, base_octant%suboctant(iii, jjj, kkk)%ptr)
                     end if
                  end do
               end do
            end do
         end if
      end select
   end subroutine compute_forces

end program tree
