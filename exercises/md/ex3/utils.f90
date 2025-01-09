module utils
    use mpi_f08
    use iso_fortran_env
    use particle
    implicit none
    
    contains
        subroutine read_data(file_in, particles, sim_name, dt, dt_out, t_end, n, in, epsilon, theta_local)
            character(len=100), intent(in) :: file_in
            character(len=100), intent(inout) :: sim_name
            real(real64), intent(out) :: dt, t_end, dt_out, theta_local, epsilon
            integer, intent(inout) :: n
            integer, intent(in) :: in
            integer :: i
            type(particle3d), allocatable, intent(inout) :: particles(:)
            
            !Read the input file in the format as the template in ics/
            open(unit=in, file=file_in, action='read')
                read(in, *) 
                read(in, '(A)') sim_name
                read(in, *)
                read(in, *)
                read(in, *) dt
                read(in, *) 
                read(in, *)
                read(in, *) dt_out
                read(in, *) 
                read(in, *)
                read(in, *) t_end
                read(in, *) 
                read(in, *)
                read(in, *) n
                read(in, *)
                read(in, *) 
                read(in, *) theta_local
                read(in, *)
                read(in, *)
                read(in, *) epsilon
                read(in, *)

                !Now I know the number of particles, so I can allocate the array
                allocate(particles(n))

                !Read particle positions
                read(in, *)
                do i = 1, n
                    read(in, *) particles(i)%p%x, particles(i)%p%y, particles(i)%p%z
                end do
                read(in, *)

                !Read particle velocities
                read(in, *)
                do i = 1, n
                    read(in, *) particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
                end do
                read(in, *)

                !Read particle masses
                read(in, *)
                do i = 1, n
                    read(in, *) particles(i)%m
                end do
            close(in)
        end subroutine read_data

        subroutine read_data_terminal(particles, sim_name, dt, dt_out, t_end, n, epsilon, theta_local)
            character(len=100), intent(inout) :: sim_name
            real(real64), intent(inout) :: dt, t_end, dt_out, theta_local, epsilon
            integer, intent(inout) :: n
            type(particle3d), allocatable, intent(inout) :: particles(:)
            
            print*, "Enter the simulation name:"
            read*, sim_name
            print*, "Enter the time step:"
            read*, dt
            print*, "Enter the output time step:"
            read*, dt_out
            print*, "Enter the final time:"
            read*, t_end
            print*, "Enter the number of particles:"
            read*, n
            print*, "Enter theta:"
            read*, theta_local
            print*, "Enter epsilon:"
            read*, epsilon

            allocate(particles(n))
        end subroutine read_data_terminal

        subroutine write_data(particles, t, n, out)
            type(particle3d), intent(in) :: particles(:)
            real(real64), intent(in) :: t
            integer, intent(in) :: n, out
            integer :: i
            do i = 1, n
                write(out, "(I6, 1X, ES12.5, 1X, ES12.5, 1X, ES12.5, 1X, ES12.5, 1X, ES12.5)") &
                i, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, particles(i)%m, t
            end do
        end subroutine write_data

        subroutine show_progress(current_time, total_time)
            real(real64), intent(in) :: current_time, total_time
            integer :: bar_width, position, i
            real(real64) :: progress

            bar_width = 50
            progress = current_time / total_time
            position = int(bar_width * progress)


            write(*,'(A)', advance='no') char(13) // "["
            do i = 1, bar_width
                if (i <= position) then
                    write(*,'(A)', advance='no') "#"
                else
                    write(*,'(A)', advance='no') " "
                end if
            end do
            write(*,'(A, F5.1, A, A, F6.1, A)', advance='no') "] ", progress*100.0, "% "
    
        end subroutine show_progress

        subroutine syncronize_data(particles, index_start, index_end, recvcounts, desplz, comsize, n, rank, ierr)
            integer, intent(in) :: comsize, n, rank
            type(particle3d), intent(inout) :: particles(n)
            integer, intent(in) :: index_start, index_end
            integer, intent(in) :: recvcounts(comsize), desplz(comsize)
            integer, intent(inout) :: ierr

            ! Gather the particles to the root process
            call mpi_gatherv(particles(index_start:index_end), recvcounts(rank+1), mpi_particle3d, particles, &
                recvcounts, desplz, mpi_particle3d, 0, mpi_comm_world, ierr)
            ! Broadcast the particles to all processes from root
            call mpi_bcast(particles, n, mpi_particle3d, 0, mpi_comm_world, ierr)

            ! Stop the program if something went wrong with the broadcast
            if (ierr /= MPI_SUCCESS) then
                print*, 'Error in gather/broadcast of data:', ierr
                call MPI_Abort(MPI_COMM_WORLD, ierr)
            end if
        end subroutine syncronize_data

end module utils