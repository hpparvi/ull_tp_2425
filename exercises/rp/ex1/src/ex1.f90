program ex1
   use geometry
   use particle
   use ics_module
   implicit none

   character(len=3) :: answer
   character(len=200) :: ics_file
   character(len=300) :: savefolder
   logical :: read_from_file
   integer :: N ! number of particles
   real :: T ! Total integration time
   real :: dt ! Time-step
   real :: e ! Smoothing length
   type(particle3d), allocatable :: bodies(:)
   integer :: M ! Number of steps
   integer :: N_snapshots ! Number of snapshots
   integer :: snap_step

   integer :: i, j, k, ii, jj, kk
   type(vector3d) :: d_vec ! Relatve position between particles
   real :: d ! Relatve distance between particles
   type(vector3d), allocatable :: a(:) ! accelerations


   call ask_ics(answer, read_from_file, ics_file)

   if (read_from_file .eqv. .false.) then
      call get_manual_ics(N, T, dt, e, bodies)
   else
      call get_ics_from_file(ics_file, N, T, dt, e, bodies)
   end if


   print*, ' '
   print*, "Let's integrate trajectories for', N, ' particles."
   print*, 'Total integration time ', T
   print*, 'Time step:', dt
   M = int(T/dt)
   print*, 'Number of steps:', M
   print*, 'Smoothing length:', e
   print*, ' '

   print*, 'How many snapshots do you want to save?'
   read*, N_snapshots

   allocate(a(N))

   ! Initialize accelerations before loop
   call compute_accelerations()
   savefolder = create_results_folder()

   snap_step = int(M/N_snapshots)
   k = 0
   do i = 1, M
      call leapfrog_step()
      if (MOD(i, snap_step) == 0) then
         call save_data(k)
         k=k+1
      end if
   end do







contains

   subroutine compute_accelerations()
      implicit none

      ! Initialize accelerations to zero
      do ii = 1, N
         a(ii)%x = 0.0
         a(ii)%y = 0.0
         a(ii)%z = 0.0
      end do

      ! Compute accelerations
      do ii = 1, N
         do jj = ii + 1, N
            d_vec = bodies(jj)%p - bodies(ii)%p
            d = norm(d_vec)
            d = d**3

            ! Update acceleration for particle i
            a(ii) = a(ii) + (bodies(jj)%m * d_vec) / d

            ! Update acceleration for particle j (Newton's third law)
            a(jj) = a(jj) - (bodies(ii)%m * d_vec) / d
         end do
      end do

   end subroutine compute_accelerations

   subroutine leapfrog_step()
      implicit none

      do jj=1, N
         bodies(jj)%v = bodies(jj)%v + (a(jj) * dt) * 0.5
         bodies(jj)%p = bodies(jj)%p + bodies(jj)%v * dt
      end do
      call compute_accelerations
      do jj=1, N
         bodies(jj)%v = bodies(jj)%v + (a(jj) * dt) * 0.5
      end do

   end subroutine leapfrog_step

   function create_results_folder() result(savepath)
      implicit none
      character(len=300) :: savepath
      character(len=300) :: temp_path
      integer :: unit_num

      call execute_command_line('pwd > tmpfile.txt')

      unit_num = 10
      open(unit=unit_num, file='tmpfile.txt', status='old', action='read')
      read(unit_num, '(A)') temp_path
      close(unit_num)

      call execute_command_line('rm tmpfile.txt')

      savepath = TRIM(temp_path) // '/data/results/'

      call execute_command_line('rm -rf ' // TRIM(savepath))
      call execute_command_line('mkdir -p ' // TRIM(savepath))

      print*, "Saving data in..."
      print *, TRIM(savepath)

   end function create_results_folder

   subroutine save_data(s)
      implicit none
      integer :: s
      character(len=300) :: filename

      write(filename, "('snapshot_', I4.4, '.dat')") s
      filename = TRIM(savefolder) // filename
      open(unit=1, file=filename, status="unknown")

      do kk=1, N
      write(1, *) bodies(kk)%p 
      end do
      close(1)

   end subroutine save_data

end program ex1
