module ics_module
   use geometry
   use particle
   implicit none
contains

   recursive subroutine ask_ics(answer, read_from_file, filename)
      implicit none
      character(len=3), intent(inout) :: answer
      character(len=200), intent(out) :: filename
      logical, intent(out) :: read_from_file

      print*, 'Do you have a .dat file with initial conditions? (y/n)'
      read*, answer

      if (trim(adjustl(answer)) == 'y') then

         read_from_file = .true.
         print*, 'Good! Provide the path:'
         read '(A)', filename
      else if (trim(adjustl(answer)) == 'n') then

         read_from_file = .false.
      else
         print*, "Please, answer 'y' or 'n'."
         call ask_ics(answer, read_from_file, filename)
      end if
   end subroutine ask_ics

   subroutine get_manual_ics(n_particles, integration_time, time_step, smoothing_length, particles)
      implicit none
      integer, intent(out) :: n_particles
      real, intent(out) :: integration_time
      real, intent(out) :: time_step
      real, intent(out) :: smoothing_length

      type(particle3d), allocatable, intent(out) :: particles(:)
      integer :: i


      print*, 'Enter the total number of particles (maximum is 5):'
      read*, n_particles

      if (n_particles>5) then
         print*, 'Those are too many to enter manually. We stick to 5.'
         n_particles = 5
      else if (n_particles<0) then
         print*, 'Nothing to do here! Exiting...'
         stop
      end if

      print*, 'Enter the total integration time:'
      read*, integration_time
      if (integration_time < 0) then
         print*, 'I was not programmed to integrate back in time. Exiting...'
         stop
      end if

      print*, 'Enter the time-step:'
      read*, time_step

      print*, 'Enter the smoothing length:'
      read*, smoothing_length


      allocate(particles(n_particles))

      do i = 1, n_particles
         print*, 'Enter the initial position (x, y, z) for particle', i, ':'
         read*, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z

         print*, 'Enter the initial velocity (v_x, v_y, v_z) for particle', i, ':'
         read*, particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
      end do
   end subroutine get_manual_ics


   subroutine get_ics_from_file(filename, n_particles, integration_time, time_step, smoothing_length, particles)
      implicit none
      character(len=200), intent(out) :: filename
      integer, intent(out) :: n_particles
      real, intent(out) :: integration_time
      real, intent(out) :: time_step
      real, intent(out) :: smoothing_length
      type(particle3d), allocatable, intent(out) :: particles(:)
      integer :: unit_num, i

      unit_num = 10
      open(unit=unit_num, file=filename, status='old', action='read')
      read(unit_num, *) n_particles
      read(unit_num, *) integration_time
      read(unit_num, *) time_step
      read(unit_num, *) smoothing_length


      allocate(particles(n_particles))

      do i = 1, n_particles
         read(unit_num, *) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
            particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
      end do

      close(unit_num)

   end subroutine get_ics_from_file

end module ics_module
