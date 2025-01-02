module mpi_types
  use geometry !Importing the geometry module
  use particle !Importing the particle module
  use barnes   !Importing the barnes module, where the Barnes-Hut algorithm is defined
  use calcs    !Importing the calcs module, where subroutines to update properties of the particles are defined
  use data     !Importing the data module, where subroutines to read and save data are defined
  use iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  use mpi_f08  !Importing the mpi library to use MPI 
  implicit none

contains

  !Subroutine to create an MPI datatype for a 3D point
  subroutine create_MPI_POINT3D(MPI_POINT3D, ierr)
    type(MPI_Datatype), intent(out) :: MPI_POINT3D !MPI 3D point
    integer, intent(inout) :: ierr !Integer variable for error handling in MPI operations

    !Create an MPI 3D point datatype
    call MPI_Type_contiguous(3, MPI_DOUBLE_PRECISION, MPI_POINT3D, ierr)

    !Check for errors in creating the MPI 3D point datatype
    if (ierr .NE. MPI_SUCCESS) then
       print*, "Error in creating MPI datatype for 3D point: ", ierr !Print error message
       call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
    end if

    !Commit the MPI 3D point datatype to make it usable in MPI operations
    call MPI_Type_commit(MPI_POINT3D, ierr)

    !Check for errors in commiting the MPI 3D point datatype
    if (ierr .NE. MPI_SUCCESS) then
        print *, "Error in commiting MPI datatype for 3D point: ", ierr !Print error message
        call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
    end if

  end subroutine create_MPI_POINT3D
  
  !Subroutine to create an MPI datatype for a 3D vector
  subroutine create_MPI_VECTOR3D(MPI_VECTOR3D, ierr)
    type(MPI_Datatype), intent(out) :: MPI_VECTOR3D !MPI 3D vector
    integer, intent(inout) :: ierr !Integer variable for error handling in MPI operations

    !Create an MPI 3D vector datatype
    call MPI_Type_contiguous(3, MPI_DOUBLE_PRECISION, MPI_VECTOR3D, ierr)

    !Check for errors in creating the MPI 3D vector datatype
    if (ierr .NE. MPI_SUCCESS) then
       print*, "Error in creating MPI datatype for 3D vector: ", ierr !Print error message
       call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
    end if

    !Commit the MPI 3D vector datatype to make it usable in MPI operations
    call MPI_Type_commit(MPI_VECTOR3D, ierr)

    !Check for errors in commiting the MPI 3D vector datatype
    if (ierr .NE. MPI_SUCCESS) then
        print *, "Error in commiting MPI datatype for 3D vector: ", ierr !Print error message
        call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
    end if
    
  end subroutine create_MPI_VECTOR3D

  !Subroutine to create an MPI datatype for a 3D particle
  subroutine create_MPI_PARTICLE3D(MPI_PARTICLE3D, ierr)
    type(MPI_Datatype), intent(out) :: MPI_PARTICLE3D !MPI 3D particle
    type(MPI_Datatype) :: MPI_POINT3D, MPI_VECTOR3D   !MPI 3D point and MPI 3D vector to create MPI 3D particle
    type(vector3d) :: vector !3D vector
    type(point3d) :: point   !3D point
    
    integer, intent(inout) :: ierr !Integer variable for error handling in MPI operations
    integer, parameter :: block_count = 4 !Number of blocks to create for a 3D particle
    integer :: block_lengths(block_count) !Contain the length of each block
    integer(MPI_ADDRESS_KIND) :: displacements(block_count) !Contain the displacement for each block
    type(MPI_Datatype) :: block_types(block_count) !Contain the types of each block
    
    !Call subroutines to create MPI datatypes for 3D points and 3D vectors
    call create_MPI_POINT3D(MPI_POINT3D, ierr)
    call create_MPI_VECTOR3D(MPI_VECTOR3D, ierr)

    !Define block lengths for the MPI 3D particle structure
    block_lengths = [1, 1, 1, 1] ![position, velocity, acceleration, mass]

    !Define block types 
    block_types = [MPI_POINT3D, MPI_VECTOR3D, MPI_VECTOR3D, MPI_DOUBLE_PRECISION] ![position, velocity, acceleration, mass]

    !Calculate displacements for each block type
    displacements = [0_MPI_ADDRESS_KIND, sizeof(point), sizeof(point) + sizeof(vector), sizeof(point) + 2*sizeof(vector)] ![position, velocity, acceleration, mass]
    
    !Create an MPI 3D particle structure
    call MPI_Type_create_struct(block_count, block_lengths, displacements, block_types, MPI_PARTICLE3D, ierr)

    !Check for errors in creating the MPI 3D particle structure
    if (ierr .NE. MPI_SUCCESS) then
       print*, "Error in creating MPI datatype for 3D particle: ", ierr !Print error message
       call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
    end if

    !Commit the MPI 3D particle structure to make it usable in MPI operations
    call MPI_Type_commit(MPI_PARTICLE3D, ierr)

    !Check for errors in commiting the MPI 3D particle structure
    if (ierr .NE. MPI_SUCCESS) then
        print *, "Error in commiting MPI datatype for 3D particle: ", ierr !Print error message
        call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
    end if

  end subroutine create_MPI_PARTICLE3D
  
  !Subroutine to deallocate MPI datatypes
  subroutine deallocate_MPI_Types(MPI_POINT3D, MPI_VECTOR3D, MPI_PARTICLE3D, ierr)
    type(MPI_Datatype) :: MPI_POINT3D, MPI_VECTOR3D, MPI_PARTICLE3D   !MPI 3D point, 3D vector and 3D particle
    integer, intent(inout) :: ierr !Integer variable for error handling in MPI operations

    !Deallocate MPI 3D point datatype
    call MPI_Type_free(MPI_POINT3D, ierr)

    !Check for errors in deallocating the MPI 3D point
    if (ierr .NE. MPI_SUCCESS) then
       print*, "Error in deallocating MPI 3D point: ", ierr !Print error message
       call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
    end if

    !Deallocate MPI 3D vector datatype
    call MPI_Type_free(MPI_VECTOR3D, ierr)

    !Check for errors in deallocating the MPI 3D vector
    if (ierr .NE. MPI_SUCCESS) then
       print*, "Error in deallocating MPI 3D vector: ", ierr !Print error message
       call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
    end if

    !Deallocate MPI 3D particle datatype
    call MPI_Type_free(MPI_PARTICLE3D, ierr)

    !Check for errors in deallocating the MPI 3D particle
    if (ierr .NE. MPI_SUCCESS) then
       print*, "Error in deallocating MPI 3D particle: ", ierr !Print error message
       call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
    end if

  end subroutine deallocate_MPI_Types
    
end module mpi_types
