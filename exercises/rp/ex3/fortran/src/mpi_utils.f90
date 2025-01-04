module mpi_utils
   use geometry
   use particle
   use mpi
   implicit none

   integer :: mpi_particle_type
   integer :: mpi_vector3d_type

contains

   subroutine create_particle_mpi_type()
      integer :: blocklengths(3)
      integer :: types(3)
      integer(MPI_ADDRESS_KIND) :: displacements(3)
      type(particle3d) :: dummy
      integer(MPI_ADDRESS_KIND) :: base, addr
      integer :: ierr

      ! Define number of elements in each block
      blocklengths = [1, 3, 3]  ! mass, position (x,y,z), velocity (x,y,z)

      ! Define MPI types for each block
      types = [MPI_REAL8, MPI_REAL8, MPI_REAL8]

      ! Calculate displacements
      call MPI_GET_ADDRESS(dummy, base, ierr)
      call MPI_GET_ADDRESS(dummy%m, addr, ierr)
      displacements(1) = addr - base

      call MPI_GET_ADDRESS(dummy%p, addr, ierr)
      displacements(2) = addr - base

      call MPI_GET_ADDRESS(dummy%v, addr, ierr)
      displacements(3) = addr - base

      ! Create the MPI struct type
      call MPI_TYPE_CREATE_STRUCT(3, blocklengths, displacements, types, mpi_particle_type, ierr)
      call MPI_TYPE_COMMIT(mpi_particle_type, ierr)
   end subroutine create_particle_mpi_type

   subroutine free_particle_mpi_type()
      integer :: ierr
      call MPI_TYPE_FREE(mpi_particle_type, ierr)
   end subroutine free_particle_mpi_type


   subroutine create_vector_mpi_type()
      integer :: blocklengths(1), types(1)
      integer(MPI_ADDRESS_KIND) :: displacements(1), base, addr
      type(vector3d) :: dummy
      integer :: ierr

      ! vector3d has 3 reals
      blocklengths(1) = 3
      types(1)        = MPI_REAL8

      call MPI_GET_ADDRESS(dummy, base, ierr)
      call MPI_GET_ADDRESS(dummy%x, addr, ierr)  ! or dummy%y, dummy%z => same block
      displacements(1) = addr - base

      call MPI_TYPE_CREATE_STRUCT(1, blocklengths, displacements, types, mpi_vector3d_type, ierr)
      call MPI_TYPE_COMMIT(mpi_vector3d_type, ierr)
   end subroutine create_vector_mpi_type

   ! subroutine free_vector_mpi_type()
   !    integer :: ierr
   !    call MPI_TYPE_FREE(mpi_vector3d_type, ierr)
   ! end subroutine free_vector_mpi_type
   ! subroutine free_vector_mpi_type()
   !    integer :: ierr
   !    call MPI_TYPE_FREE(mpi_vector3d_type, ierr)
   ! end subroutine free_vector_mpi_type

end module mpi_utils
