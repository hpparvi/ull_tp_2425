program ex1
   use geometry
   use particle
   use ics_module
   implicit none

   character(len=3) :: answer
   character(len=200) :: filename
   logical :: read_from_file
   integer :: N ! number of particles
   real :: T ! Total integration time
   real :: dt ! Time-step
   real :: e ! Smoothing length
   type(particle3d), allocatable :: particles(:)
   integer :: M ! Number of steps

   integer :: i, j, k
   type(vector3d) :: d ! Relatve position between particles
   real, dimension(:,:), allocatable :: a ! accelerations


   call ask_ics(answer, read_from_file, filename)

   if (read_from_file .eqv. .false.) then
      call get_manual_ics(N, T, dt, e, particles)
   else
      call get_ics_from_file(filename, N, T, dt, e, particles)
   end if

   print*, ' '
   print*, 'Integrating trajectories for', N, ' particles.'
   print*, 'Total integration time ', T
   print*, 'Time step:', dt 
   M = int(T/dt)
   print*, 'Number of steps:', M 
   print*, 'Smoothing length:', e 
   print*, ' '

   allocate(a(N, 3))
   a = 0






end program ex1
