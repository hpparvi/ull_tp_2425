module particle
  use geometry   
  use mpi
  implicit none

  
  type :: particle3d
      type(vector3d) :: p   ! position
      type(vector3d) :: v   ! velocity
      type(vector3d) :: a   ! acceleration
      real(kind=8) :: m     ! mass
  end type particle3d

  ! Variables para los tipos MPI
  integer :: mpi_type_particle3d, mpi_type_vector3d
  integer(MPI_ADDRESS_KIND), dimension(4) :: array_of_displacements
  integer, dimension(4) :: array_of_types         
  integer, dimension(4) :: array_of_blocklengths
  integer :: ierror

contains

  ! Subrutina para definir los tipos MPI
  subroutine define_mpi_types()
    ! Declaraciones de variables
    integer(MPI_ADDRESS_KIND), dimension(3) :: displs
    integer, dimension(3) :: types, blocklengths
    integer(MPI_ADDRESS_KIND), dimension(4) :: displs_particle
    integer, dimension(4) :: types_particle, blocklengths_particle
    
    ! Definimos el tipo MPI para los vectores
    displs = [0_MPI_ADDRESS_KIND, 1_MPI_ADDRESS_KIND * 8, 2_MPI_ADDRESS_KIND * 8]
    types = [MPI_REAL8, MPI_REAL8, MPI_REAL8]
    blocklengths = [1, 1, 1]

    call MPI_TYPE_CREATE_STRUCT(3, blocklengths, displs, types, mpi_type_vector3d, ierror)
    call MPI_TYPE_COMMIT(mpi_type_vector3d, ierror)

    ! Definimos el tipo MPI para las particulas
    displs_particle = [0_MPI_ADDRESS_KIND, 24_MPI_ADDRESS_KIND, 48_MPI_ADDRESS_KIND, 72_MPI_ADDRESS_KIND]
    types_particle = [mpi_type_vector3d, mpi_type_vector3d, mpi_type_vector3d, MPI_REAL8]
    blocklengths_particle = [1, 1, 1, 1]

    call MPI_TYPE_CREATE_STRUCT(4, blocklengths_particle, displs_particle, types_particle, mpi_type_particle3d, ierror)
    call MPI_TYPE_COMMIT(mpi_type_particle3d, ierror)

  end subroutine define_mpi_types

end module particle
