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
   type(particle3d), allocatable :: bodies(:)

   type(octant), pointer :: head, temp_cell

   integer :: i, j, k, ii, jj, kk

   call read_config(ics_file, G, T, dt, e)

   call get_ics_from_file(ics_file, N, bodies)

   ! Create folder to save data
   savefolder = create_snapshot_folder()
   ! Save positions at time=0
   call save_data(0, savefolder, bodies, N)

   allocate(head)
   call get_cell_range(head, N, bodies)
   head%category = 0
   call nullify_pointers(head)


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

end program tree
