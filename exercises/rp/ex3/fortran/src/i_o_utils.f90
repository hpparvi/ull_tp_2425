module i_o_utils
   use geometry
   use particle
   use mpi
   use mpi_utils

   implicit none

contains
   subroutine read_config(config_path, ics_path, savefolder,  N_snapshots, G, integration_time, time_step, smoothing_length, theta)
      implicit none
      character(len=500), intent(out) :: ics_path
      character(len=500), intent(out) :: savefolder
      real, intent(out) :: G
      real, intent(out) :: integration_time
      real, intent(out) :: time_step
      real, intent(out) :: smoothing_length
      real, intent(out) :: theta
      integer, intent(out) :: N_snapshots

      character(len=500) :: config_path, line, key, value
      integer :: ios, idx

      ! Open the configuration file
      open(unit=10, file=config_path, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         print *, "Error opening file: ", config_path
         stop
      end if

      ! Loop through the lines of the configuration file
      do
         read(10, '(A)', iostat=ios) line
         if (ios < 0) exit ! Exit the loop at EOF
         if (trim(adjustl(line)) == '' .or. line(1:1) == '#') cycle ! Skip blank or comment lines

         idx = index(line, "=")
         if (idx > 0) then
            key = trim(adjustl(line(1:idx-1)))
            value = trim(adjustl(line(idx+1:)))

            select case (key)
             case ("G")
               read(value, *, iostat=ios) G
             case ("T")
               read(value, *, iostat=ios) integration_time
             case ("DT")
               read(value, *, iostat=ios) time_step
             case ("EPSILON")
               read(value, *, iostat=ios) smoothing_length
             case ("THETA")
               read(value, *, iostat=ios) theta
            case ("N_SNAPSHOTS")
               read(value, *, iostat=ios) N_snapshots
             case ("ICS_FILE")
               ics_path = value
            case ("SAVE_FOLDER")
               savefolder = value
             case default
               print*, "Unknown key: ", key
            end select
         end if
      end do

      close(10)
   end subroutine read_config


   subroutine get_ics_from_file(filename, n_particles, particles)
      ! Provide file of ICS
      implicit none
      character(len=500), intent(in) :: filename
      integer, intent(out) :: n_particles
      type(particle3d), allocatable, intent(out) :: particles(:)
      integer :: unit_num, i

      unit_num = 10
      open(unit=unit_num, file=filename, status='old', action='read')
      read(unit_num, *) n_particles

      allocate(particles(n_particles))

      do i = 1, n_particles
         read(unit_num, *) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
            particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
      end do

      close(unit_num)

   end subroutine get_ics_from_file


   subroutine create_snapshot_folder(savefolder)
      implicit none
      character(len=300), intent(out) :: savefolder

      call execute_command_line('rm -rf ' // TRIM(savefolder))
      call execute_command_line('mkdir -p ' // TRIM(savefolder))

      print*, "Saving data in..."
      print *, TRIM(savefolder)

   end subroutine create_snapshot_folder

   subroutine save_data(s, savefolder, bodies, N)
      implicit none
      integer :: s
      integer :: N
      integer :: i
      character(len=500) :: savefolder
      character(len=500) :: filename
      type(particle3d), allocatable :: bodies(:)

      write(filename, "('snapshot_', I4.4, '.dat')") s
      filename = TRIM(savefolder) // filename
      open(unit=1, file=filename, status="unknown")

      do i=1, N
         write(1, *) bodies(i)%p
      end do
      close(1)

   end subroutine save_data

   subroutine get_N_snapshots(N_snapshots, N_timesteps)

      integer, intent(out):: N_snapshots
      integer, intent(in) :: N_timesteps

      real :: fact
      
      fact = real(N_timesteps) / N_snapshots

      if (fact <= 1) then 
         N_snapshots = N_timesteps
      else 
         N_snapshots = int(N_timesteps/fact)
      end if

   end subroutine get_N_snapshots

end module i_o_utils
