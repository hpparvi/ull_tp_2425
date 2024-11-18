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
   integer :: snap_step  ! number of time-steps to ignore between saves

   type(octant), pointer :: root_octant, temp_octant

   integer :: i, j, k  ! dummy loop varables



   call read_config(ics_file, N_snapshots, G, T, dt, e, theta)
   N_time_steps = int(T/dt)

   call get_ics_from_file(ics_file, N, bodies)
   savefolder = create_snapshot_folder()
   call save_data(0, savefolder, bodies, N)  ! Save positions at time=0

   call get_N_snapshots(N_snapshots, N_time_steps)

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

   ! call compute_masses(tronco, N, bodies)

   ! call nullify_pointers(temp_cell)
   ! ! allocate(temp_cell)

   ! allocate(a(N))

   ! call compute_forces(head, N, bodies, e, theta, a)

   ! N_snapshots = 100
   ! k = 1
   ! snap_step = int(M/N_snapshots)

   ! print*, 'snap_step', snap_step

   ! do kk = 1, M

   !    print*, 'STEP', k

   !    do jj=1, N
   !       bodies(jj)%v = bodies(jj)%v + (a(jj) * dt) * 0.5
   !       bodies(jj)%p = bodies(jj)%p + bodies(jj)%v * dt
   !    end do

   !    CALL delete_tree(head)
   !    CALL get_cell_range(head, N, bodies)
   !    head%category = 0
   !    call nullify_pointers(head)
   !    call nullify_pointers(temp_cell)

   !    do i = 1,N
   !       call find_cell(head,temp_cell,bodies(i))
   !       call place_cell(temp_cell,bodies(i),i)
   !    end do

   !    CALL delete_empty_leaves(head)
   !    CALL compute_masses(head,  N, bodies)

   !    do ii=1, N
   !       a(ii)%x = 0.
   !       a(ii)%y = 0.
   !       a(ii)%z = 0.
   !    end do
   !    CALL compute_forces(head, N, bodies, e, theta, a)
   !    ! do jj=1, N
   !    !    bodies(jj)%v = bodies(jj)%v + (a(jj) * dt) * 0.5
   !    ! end do

   !    if (MOD(kk-1, snap_step) == 0) then
   !       call save_data(k, savefolder, bodies, N)
   !       k=k+1
   !    end if


   ! end do

   ! print*, 'Done!'

contains


end program tree


 !a = 0.0



 !belong_bool = belongs(bodies(1), head)

 ! print*, theta
 ! print*, head%range%r_min
 ! print*, head%range%r_max


 ! print*, G
 ! print*, T
 ! print*, dt
 ! print*, e
 ! print*, ics_file
 ! print*, N
 !print*, 'shape of bodies%p', shape(bodies)
 ! do i=1, 10
 ! print*, bodies(i)%m
 ! end do
 ! print*, bodies(0)%p
 ! print*, bodies(0)%v
