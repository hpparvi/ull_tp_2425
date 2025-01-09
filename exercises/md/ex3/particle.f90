module particle
    use mpi_f08
    use iso_fortran_env
    use geometry

    implicit none
    

    !Particle type: it has position, velocity, aceleration and mass
    type :: particle3d
        type(point3d) :: p
        type(vector3d) :: v
        type(vector3d) :: a
        real(real64) :: m
    end type particle3d

    type(MPI_Datatype) :: mpi_particle3d
    

    contains
    !Create the MPI type for the particle3d type
    subroutine create_mpi_particle_type(mpi_particle3d, ierr)
        type(MPI_Datatype) :: mpi_particle3d
        type(MPI_Datatype) :: mpi_point3d, mpi_vector3d
        type(vector3d) :: vect3d
        type(point3d) :: pont3d
        integer, intent(inout) :: ierr
        integer, parameter :: block_count = 4
        integer :: block_lengths(block_count)
        integer(kind=MPI_ADDRESS_KIND) :: displacements(block_count)
        type(MPI_Datatype) :: block_types(block_count)

        ! Create the MPI types for the point3d and vector3d types
        call MPI_Type_contiguous(3, mpi_real8, mpi_point3d, ierr)
        call MPI_Type_commit(mpi_point3d, ierr)

        ! Create the MPI types for the point3d and vector3d types
        call MPI_Type_contiguous(3, mpi_real8, mpi_vector3d, ierr)
        call MPI_Type_commit(mpi_vector3d, ierr)


        block_lengths = [1, 1, 1, 1]

        block_types = [mpi_point3d, mpi_vector3d, mpi_vector3d, mpi_real8]

        displacements = [0_mpi_address_kind, &
                        sizeof(pont3d), &
                        sizeof(pont3d) + sizeof(vect3d), &
                        sizeof(pont3d) + 2*sizeof(vect3d)]


        call mpi_type_create_struct(block_count, block_lengths, displacements, block_types, mpi_particle3d, ierr)
        call mpi_type_commit(mpi_particle3d, ierr)
    end subroutine create_mpi_particle_type


        !Reset acceleration of the particles to zero
        subroutine reset_a(particles)
            type(particle3d), intent(inout) :: particles(:)
            integer :: i
            do i = 1, size(particles)
                particles(i)%a = vector3d(0.0, 0.0, 0.0)
            end do
        end subroutine reset_a

        !Update the position of the particles
        subroutine update_pos(particles, dt)
            type(particle3d), intent(inout) :: particles(:)
            real(real64), intent(in) :: dt
            integer :: i
            do i = 1, size(particles)
                particles(i)%p = particles(i)%p + particles(i)%v * dt
            end do
        end subroutine update_pos

        !Update the velocity of the particles in a middle step
        subroutine update_vel(particles, dt)
            type(particle3d), intent(inout) :: particles(:)
            real(real64), intent(in) :: dt
            integer :: i
            do i = 1, size(particles)
                particles(i)%v = particles(i)%v + particles(i)%a * dt * real(0.5, kind=real64)
            end do
        end subroutine update_vel

        !Update the acceleration of the particles
        subroutine update_a(particles)
            type(particle3d), intent(inout) :: particles(:)
            integer :: i, j
            real(real64) :: r2, r3
            type(vector3d):: rji
            
            do i = 1, size(particles)
                do j = i+1, size(particles)
                    rji = particles(j)%p - particles(i)%p
                    r2 = normsquare(rji)
                    r3 = r2 * sqrt(r2)
                    particles(i)%a = particles(i)%a + particles(j)%m * rji / r3
                    particles(j)%a = particles(j)%a - particles(i)%m * rji / r3
                end do
            end do
        end subroutine update_a
end module particle