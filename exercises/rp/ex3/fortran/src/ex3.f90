program tree
   use geometry
   use particle
   use i_o_utils
   use octree
   use mpi_utils

   use mpi

   implicit none

   ! MPI variables
   integer :: ierr, rank, size

   ! dummy loop varables
   integer :: i, j, k, ii, jj, global_i
   integer :: step, snap_idx
   integer :: start_idx, end_idx


   character(len=500) :: config_path, ics_file, savefolder

   integer :: N, N_snapshots, N_time_steps, snap_step, local_N
   real :: G, T, dt, e, theta

   type(particle3d), allocatable :: bodies(:)

   integer :: remainder
   integer, allocatable :: counts(:), displs(:)

   type(octant), pointer :: root_octant, temp_octant
   type(vector3d), allocatable :: a(:) ! accelerations

   real :: start_time, end_time, elapsed_time

   ! progress marker stuff
   integer :: progress_marker, progress_marker_aux
   progress_marker = 1  ! every 1%
   progress_marker_aux = progress_marker

   ! Initialize the MPI environment
   call MPI_INIT(ierr)
   ! Get the rank (ID) of this process
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   ! Get the total number of processes
   call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

   ! Create custom MPI type for particle3d
   call create_particle_mpi_type()
   call create_vector_mpi_type()

   ! Only first process reads config file
   if (rank == 0) then
      call get_command_argument(1, config_path)
      call read_config(config_path, ics_file, savefolder, N_snapshots, G, T, dt, e, theta)
   end if

   ! Broadcast the simulation parameters to all processes
   call MPI_Bcast(G, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(T, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(e, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(theta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(N_snapshots, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(savefolder, len(savefolder), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

   ! Get N_snapshots
   N_time_steps = int(T/dt)
   call get_N_snapshots(N_snapshots, N_time_steps)
   snap_step = N_time_steps/N_snapshots


   ! Only first process creates snapshot folder
   ! Only first process reads ics from .dat file
   ! Only first process saves data
   ! (get_ics_... func allocates variable bodies)
   if (rank == 0) then
      call create_snapshot_folder(savefolder)
      call get_ics_from_file(ics_file, N, bodies)
      call save_data(0, savefolder, bodies, N)
   end if

   call MPI_Bcast(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   if (rank == 0) then
      print*, ' '
      print*, "Number of particles:", N
      print*, 'Total integration time:', T
      print*, 'Time step:', dt
      print*, 'Number of steps:', N_time_steps
      print*, 'Number of snapshots:', N_snapshots
   end if


   ! Allocate bodies in other processes and broadcast
   if (.not. allocated(bodies)) allocate(bodies(N))
   call MPI_Bcast(bodies, N, mpi_particle_type, 0, MPI_COMM_WORLD, ierr)

   ! Accounting for non equal number of particles per process
   ! Allocate counts and displacements with 1-based indexing
   allocate(counts(size), displs(size))
   remainder = mod(N, size)

   do i=1, size
      counts(i) = N / size
   end do
   do i=1, remainder
      counts(i) = counts(i) + 1
   end do

   displs(1) = 0
   do i=2, size
      displs(i) = displs(i-1) + counts(i-1)
   end do

   local_N = counts(rank+1)  ! portion that rank "owns"

   start_idx = displs(rank+1) + 1
   end_idx   = displs(rank+1) + counts(rank+1)

   allocate(a(N))


   ! check number of particles per process
   print *, 'Rank', rank, 'has', local_N, 'particles'


   ! Every process builds the tree for the *full* set:
   allocate(root_octant)

   call give_root_octant_a_range(root_octant, N, bodies)
   root_octant%category = 0  ! No bodies in `root_octant`, yet
   call nullify_suboctant_pointers(root_octant)  ! setting subcell pointers to null
   do i = 1,N
      call find_octant(root_octant,temp_octant,bodies(i))
      call place_particle_in_octant(temp_octant, bodies(i),i)
   end do
   call delete_empty_leaves(root_octant)
   call compute_masses(root_octant) !

   ! Initialize accelerations to zero
   do i = 1, N
      a(i)%x = 0.0
      a(i)%y = 0.0
      a(i)%z = 0.0
   end do

   ! Each rank calls compute_forces FOR ITS CHUNK ONLY
   do i = start_idx, end_idx
      call compute_forces(i, root_octant)
   end do



   ! Start the timer
   start_time = mpi_wtime()
   snap_idx = 1
   do step = 1, N_time_steps

      ! First half-step: each rank updates only its local subset
      do i = start_idx, end_idx
         bodies(i)%v = bodies(i)%v + a(i)*dt*0.5
         bodies(i)%p = bodies(i)%p + bodies(i)%v*dt
      end do

      ! comunicate updates
      call MPI_Allgatherv( &
         bodies(start_idx), counts(rank+1), mpi_particle_type, &
         bodies,            counts,         displs, mpi_particle_type, &
         MPI_COMM_WORLD, ierr )

      ! Rebuild the tree (whole tree for each rank)
      call delete_tree(root_octant)
      call give_root_octant_a_range(root_octant, N, bodies)
      root_octant%category = 0
      call nullify_suboctant_pointers(root_octant)
      do jj = 1, N
         call find_octant(root_octant, temp_octant, bodies(jj))
         call place_particle_in_octant(temp_octant, bodies(jj), jj)
      end do
      call delete_empty_leaves(root_octant)
      call compute_masses(root_octant)

      ! Zero local accelerations
      do ii = 1, N
         a(ii)%x = 0.0
         a(ii)%y = 0.0
         a(ii)%z = 0.0
      end do

      ! Compute new forces in each particle subset
      do i = start_idx, end_idx
         call compute_forces(i, root_octant)
      end do

      ! update velocities in each particle subset
      do i = start_idx, end_idx
         bodies(i)%v = bodies(i)%v + a(i)*dt*0.5
      end do

      ! Allgather velocities again
      call MPI_Allgatherv( &
         bodies(start_idx), counts(rank+1), mpi_particle_type, &
         bodies,            counts,         displs, mpi_particle_type, &
         MPI_COMM_WORLD, ierr )

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! Save data only in one process
      if (rank == 0) then
         if (MOD(step-1, snap_step) == 0) then
            call save_data(snap_idx, savefolder, bodies, N)
            snap_idx = snap_idx + 1
         end if
      end if


   end do

   ! End the timer after the loop
   end_time = mpi_wtime()
   elapsed_time = end_time - start_time

   ! Print elapsed time (only done by rank 0)
   if (rank == 0) then
      print *, "Done! Elapsed CPU time:", elapsed_time, "s"
   end if

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)

   ! Free the custom MPI type
   call free_particle_mpi_type()
   ! Finalize the MPI environment
   call MPI_FINALIZE(ierr)
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
