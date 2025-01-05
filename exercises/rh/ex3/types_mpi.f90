module types_mpi
  use, intrinsic :: iso_fortran_env
  use mpi_f08
  use geometry
  use particle
  implicit none
  
  ! Module-level variables for MPI datatypes
  TYPE(MPI_Datatype) :: mpi_point3d = MPI_DATATYPE_NULL
  TYPE(MPI_Datatype) :: mpi_vector3d = MPI_DATATYPE_NULL
  TYPE(MPI_Datatype) :: mpi_particle3d = MPI_DATATYPE_NULL
  
contains
  subroutine create_mpi_types()
    type(vector3d) :: dummy_vector
    type(point3d) :: dummy_point
    
    !integer, parameter :: num_blocks = 3 !things inside vector3d
    integer :: ierr
    
    ! variables to define mpi types  
    integer :: lengths(0:2)
    integer(kind=MPI_Address_Kind) :: displacements(0:2)
    type(MPI_Datatype) :: types(0:2)
    
    ! Create point3d type for mpi
    lengths = [1, 1, 1]
    displacements = [0 * SIZEOF(dummy_point%x), &
         & sizeof(dummy_point%x), &
         & sizeof(dummy_point%x) + sizeof(dummy_point%y)]
    types = [MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION]
    
    CALL MPI_Type_create_struct(3, lengths, displacements, types, mpi_point3d, ierr)
    if (ierr.ne.0) then
      print *, 'Error creating mpi type ', ierr
      call MPI_Abort(MPI_COMM_WORLD, ierr)
    end if
    CALL MPI_Type_commit(mpi_point3d, ierr)
    if (ierr.ne.0) then
      print *, 'Error creating mpi type ', ierr
       call MPI_Abort(MPI_COMM_WORLD, ierr)
    end if
    
    ! Create vector3d type for mpi
    displacements = [0 * SIZEOF(dummy_vector%x), &
         & sizeof(dummy_vector%x), &
         & sizeof(dummy_vector%x) + sizeof(dummy_vector%y)]
    
    CALL MPI_Type_create_struct(3, lengths, displacements, types, mpi_vector3d, ierr)
    if (ierr.ne.0) then
      print *, 'Error creating mpi type ', ierr
      call MPI_Abort(MPI_COMM_WORLD, ierr)
    end if
    CALL MPI_Type_commit(mpi_vector3d, ierr)
    if (ierr.ne.0) then
      print *, 'Error creating mpi type ', ierr
      call MPI_Abort(MPI_COMM_WORLD, ierr)
    end if

    ! Create particle3d type for mpi
    displacements = [ 0*sizeof(dummy_point), &
       & sizeof(dummy_point), &
       & sizeof(dummy_point) + sizeof(dummy_vector)]
    types = [ mpi_point3d, &
       & mpi_vector3d, &
       & mpi_DOUBLE_PRECISION]

    CALL MPI_Type_create_struct(3, lengths, displacements, types, mpi_particle3d, ierr)
    if (ierr.ne.0) then
      print *, 'Error creating mpi type ', ierr
      call MPI_Abort(MPI_COMM_WORLD, ierr)
    end if
    CALL MPI_Type_commit(mpi_particle3d, ierr)
    if (ierr.ne.0) then
        print *, 'Error creating mpi type ', ierr
        call MPI_Abort(MPI_COMM_WORLD, ierr)
    end if

  end subroutine create_mpi_types
 
  subroutine free_mpi_types()
        INTEGER :: ierr
        IF (mpi_particle3d /= MPI_DATATYPE_NULL) THEN
            CALL MPI_TYPE_FREE(mpi_particle3d, ierr)
        END IF
        
        IF (mpi_vector3d /= MPI_DATATYPE_NULL) THEN
            CALL MPI_TYPE_FREE(mpi_vector3d, ierr)
        END IF
        
        IF (mpi_point3d /= MPI_DATATYPE_NULL) THEN
            CALL MPI_TYPE_FREE(mpi_point3d, ierr)
        END IF
  end subroutine free_mpi_types
    

end module types_mpi
