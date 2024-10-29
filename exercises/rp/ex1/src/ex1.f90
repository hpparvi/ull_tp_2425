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
 ! integer :: i, j
 ! integer :: n
 ! real :: dt, t_end, t, dt_out, t_out
 ! type(particle3d), allocatable :: particles(:)
 ! type(vector3d) :: rji
 ! real :: r2, r3

 ! ! Read input values for time step, output interval, end time, and number of particles
 ! read*, dt
 ! read*, dt_out
 ! read*, t_end
 ! read*, n

 !! Allocate array of particles
 ! allocate(particles(n))

 ! ! Read initial data for each particle (mass, position, velocity)
 ! do i = 1, n
 !     read*, particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
 !           particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
 ! end do

 ! ! Initialize accelerations to zero
 ! do i = 1, n
 !     particles(i)%v = vector3d(0.0, 0.0, 0.0)
 ! end do

 ! ! Initial force calculation
 ! do i = 1, n
 !     do j = i + 1, n
 !         rji = particles(j)%p - particles(i)%p
 !         r2 = dot(rji, rji)
 !         r3 = r2 * sqrt(r2)

 !         particles(i)%v = particles(i)%v + (particles(j)%m / r3) * rji
 !         particles(j)%v = particles(j)%v - (particles(i)%m / r3) * rji
 !     end do
 ! end do

 ! ! Integration loop using leapfrog
 ! t_out = 0.0
 ! do t = 0.0, t_end, dt

 !     ! Update velocities by half a step
 !     do i = 1, n
 !         particles(i)%v = particles(i)%v + (0.5 * dt) * particles(i)%v
 !     end do

 !     ! Update positions by a full step
 !     do i = 1, n
 !         particles(i)%p = particles(i)%p + dt * particles(i)%v
 !     end do

 !     ! Reset accelerations to zero
 !     do i = 1, n
 !         particles(i)%v = vector3d(0.0, 0.0, 0.0)
 !     end do

 !     ! Recalculate forces
 !     do i = 1, n
 !         do j = i + 1, n
 !             rji = particles(j)%p - particles(i)%p
 !             r2 = dot(rji, rji)
 !             r3 = r2 * sqrt(r2)

 !             particles(i)%v = particles(i)%v + (particles(j)%m / r3) * rji
 !             particles(j)%v = particles(j)%v - (particles(i)%m / r3) * rji
 !         end do
 !     end do

 !     ! Update velocities by another half a step
 !     do i = 1, n
 !         particles(i)%v = particles(i)%v + (0.5 * dt) * particles(i)%v
 !     end do

 !     ! Output positions at the specified interval
 !     t_out = t_out + dt
 !     if (t_out >= dt_out) then
 !         do i = 1, n
 !             print*, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
 !         end do
 !         t_out = 0.0
 !     end if

 ! end do




! program nbody_program
!    use geometry
!    use particle
!    implicit none

!    type ( particle3d ) :: b1

!    b1%v = vector3d(1, 2, 3)
!    b1%p = point3d(1, 2, 3)
!    b1%m = 10

!    print *, b1

! end program nbody_program

