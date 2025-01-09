program tree
    use mpi_f08
    use iso_fortran_env
    use particle
    use geometry
    use calculations
    use utils
    use barnes_hut

    implicit none
    
    integer :: rank, comsize, ierr, local_n, extra_n, index_start, index_end
    integer, dimension(:), allocatable :: recvcounts, desplz
    integer :: i, n
    integer, parameter :: in = 1, out = 2
    real(real64) :: dt, t_end, dt_out, t_out, exec_time, t=0.0
    real(real64) :: epsilon, theta_local
    integer(int64) :: start_time,  end_time, time_rate
    type(particle3d), allocatable :: particles(:)
    character(len=100) :: sim_name, file_out
    ! Input file of initial conditions
    character(len=100) :: file_in
    type(cell), pointer:: head
    

    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, comsize, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    
    ! Create the MPI type for the particle3d type to share the data between the processes
    call create_mpi_particle_type(mpi_particle3d, ierr)
    
    ! Only the root process reads the initial conditions
    if (rank == 0) then
        ! Check if an input file was provided or use terminal input
        if (command_argument_count() == 1) then
            !Get the input file name
            call get_command_argument(1, file_in)
            print*, 'Using input file: ', trim(file_in)

            ! Read the input file in the format as the template in ics/
            call read_data(file_in, particles, sim_name, dt, dt_out, t_end, n, in, epsilon, theta_local)
        
        ! If no input file was provided, use terminal input
        else if (command_argument_count() < 1) then
            print*, 'Using terminal input, printing positions, velocities and masses'
            call read_data_terminal(particles, sim_name, dt, dt_out, t_end, n, epsilon, theta_local)
        
            do i = 1, n
                print*, "Particle ", i, " position: x, y, z"
                read*, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
                print*, "Particle ", i, " velocity: x, y, z"
                read*, particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
                print*, "Particle ", i, " mass"
                read*, particles(i)%m
                print*, "____________________________________"
            end do
            
            ! Write the initial conditions to a file
            open(unit=in, file='ics/ic_' // trim(adjustl(sim_name)) // '.txt', action='write')
                write(in, '(A)') '! Simulation name:'
                write(in, '(A)') sim_name
                write(in, *)
                write(in, '(A)') '! Time step:'
                write(in, *) dt
                write(in, *)
                write(in, '(A)') '! Output time step:'
                write(in, *) dt_out
                write(in, *)
                write(in, '(A)') '! Final time:'
                write(in, *) t_end
                write(in, *)
                write(in, '(A)') '! Number of particles:'
                write(in, *) n
                write(in, *)
                write(in, '(A)') '! Theta:'
                write(in, *) theta_local
                write(in, *)
                write(in, '(A)') '! Epsilon (softening length):'
                write(in, *) epsilon
                write(in, *)

                ! Write particle positions
                write(in, '(A)') '! Positions (x,y,z) (each line, a particle):'
                do i = 1, n
                write(in, *) particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
                end do
                write(in, *)

                ! Write particle velocities
                write(in, '(A)') '! Velocity (vx, vy, vz) (each line, a particle):'
                do i = 1, n
                write(in, *) particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
                end do
                write(in, *)

                ! Write particle masses
                write(in, '(A)') '! Mass (each line, a particle):'
                do i = 1, n
                write(in, *) particles(i)%m
                end do
            close(in)
        ! If more than one input file was provided, stop the program
        else
            print*, "Error: input file or terminal input not provided"
            print*, "Use for input file: ex1 <input_file>"
            print*, "Use for terminal input: ex1"
            stop
        end if

        ! File to write the output, adjustl is used to put initial spaces at final, 
        ! trim is used to remove the spaces at the end of the string
        file_out = 'output/' // trim(adjustl(sim_name)) // '.dat'

        print*, "____________________________________"
        print*, "Simulation name: ", sim_name
        print*, "Time step: ", dt
        print*, "Output time step: ", dt_out
        print*, "Final time: ", t_end
        print*, "Number of particles: ", n
        print*, "Output file: ", file_out
        print*, "Theta: ", theta_local
        print*, "____________________________________"
        print*, "Starting simulation..."

        call system_clock(start_time, time_rate)
    end if

    
    
    ! Broadcast the simulation parameters to all processes
    call mpi_bcast(dt, 1, mpi_real8, 0, mpi_comm_world, ierr)
    call mpi_bcast(dt_out, 1, mpi_real8, 0, mpi_comm_world, ierr)
    call mpi_bcast(t_end, 1, mpi_real8, 0, mpi_comm_world, ierr)
    call mpi_bcast(n, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_bcast(epsilon, 1, mpi_real8, 0, mpi_comm_world, ierr)
    call mpi_bcast(theta_local, 1, mpi_real8, 0, mpi_comm_world, ierr)
    call mpi_bcast(file_out, 100, mpi_character, 0, mpi_comm_world, ierr)

    ! Stop the program if something went wrong with the broadcast
    if (ierr /= MPI_SUCCESS) then
        print*, 'Error in MPI_Bcast:', ierr
        call MPI_Abort(MPI_COMM_WORLD, ierr)
    end if

    ! Allocate the particles array in all processes except the root because it was already allocated
    if (rank /= 0) then
        allocate(particles(n))
    end if

    allocate(head)
    
    ! Broadcast the particles array to all processes from root
    call mpi_bcast(particles, n, mpi_particle3d, 0, mpi_comm_world, ierr)

    ! Number of particles per process 
    local_n = n / comsize
    ! Number of extra particles for the first processes
    extra_n = mod(n, comsize)
    
    ! Calculate the start and end index of the particles for each process
    ! First processes will have 1 more particle if there are extra particles
    if (rank < extra_n) then
        index_start = rank * (local_n + 1) + 1
        index_end = index_start + local_n
    else
        index_start = rank * local_n + extra_n + 1
        index_end = index_start + local_n - 1
    end if

    ! Array to store the number of particles each process will receive, each element represents a rank
    allocate(recvcounts(comsize))
    do i = 1, comsize
        recvcounts(i) = local_n
    end do
    
    ! Adding extra particles to the previous array
    do i = 1, extra_n
        recvcounts(i) = recvcounts(i) + 1
    end do
    
    ! Array to store the displacement of the particles for each process, each element represents a rank
    allocate(desplz(comsize))
    desplz(1) = 0
    do i = 2, comsize
        desplz(i) = desplz(i-1) + recvcounts(i-1)
    end do

    ! Create the octree
    call barnes_hut_tree(particles, head, n)
    ! Calculate the masses of the cells
    call calculate_masses(particles, head)
    
    ! Calculate the forces for the particles, each process will calculate the forces for its particles
    call reset_a(particles(index_start:index_end))
    call calculate_forces(particles, index_start, index_end, head, epsilon, theta_local, rank)

    ! Start the simulation
    t_out = 0.0
    if (rank==0) then
        write(*,'(A)', advance='no') char(13) // "["
    end if
    open(unit=out, file=file_out, action='write')
        if (rank == 0) then
            write(out, "(A6, 1X, A12, 1X, A12, 1X, A12, 1X, A12, 1X, A12, 1X)") &
                "#Index", "X", "Y", "Z", "Mass", "Time"
        end if
        do while (t < t_end)

            call update_vel(particles(index_start:index_end), dt)
            call update_pos(particles(index_start:index_end), dt)
            
            ! Synchronize the data between the processes, so after that all processes have the updated positions and &
            ! velocities of all particles
            call syncronize_data(particles, index_start, index_end, recvcounts, desplz, comsize, n, rank, ierr)

            call borrar_tree(head)
            call barnes_hut_tree(particles, head, n)
            call calculate_masses(particles, head)

            call reset_a(particles(index_start:index_end))
            call calculate_forces(particles, index_start, index_end, head, epsilon, theta_local, rank)
            ! it is not necessary to synchronize the data here because the update of the positions and velocities is of &
            ! a particle does not depend on the other particles (just the forces depends on the other particles)
            call update_vel(particles(index_start:index_end), dt)

            
            t_out = t_out + dt
            if (rank == 0) then 
                if (t_out >= dt_out) then
                    call write_data(particles, t, n, out)
                    t_out = 0.0
                end if
            end if

            t = t + dt

            if (rank == 0) then
                !Progress bar
                call show_progress(t, t_end)
            end if
        end do

        if (rank == 0) then
            call system_clock(end_time)
            exec_time = real(end_time - start_time, kind=real64) / real(time_rate, kind=real64)
            write(out, "(A, F12.6)") "#Simulation execution time: ", exec_time
        end if
    close(out)
    if (rank == 0) then
        print*, " "
        print*, "Done!"
        if (exec_time > 60.0) then
            print "(A, F6.3, A)", "Simulation execution time: ", exec_time / 60.0, " minutes"
        else
            print "(A, F6.3, A)", "Simulation execution time: ", exec_time, " seconds"
        end if
        print*, "____________________________________"
    end if
    call mpi_finalize(ierr)
end program tree
