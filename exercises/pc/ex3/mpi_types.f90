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
  !Subroutine to create an MPI datatype for a 3D vector
  subroutine mpi_vector3d(mpi_vector3d)
    type(MPI_datatype), intent(out) :: mpi_vector3d
    integer(int64) :: ierr

    call MPI_Type_contiguous(3, MPI_REAL8, mpi_vector3d, ierr)

    !Check for errors in MPI_Type_contiguos
    if (ierr /= MPI_SUCCESS) then
       print*, "Error in creating MPI type for 3D vector:", ierr
       call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminates all MPI processes
    end if
    
    call MPI_Type_commit(mpi_vector3d, ierr)

    !Check for errors in MPI_Type_commit
    if (ierr /= MPI_SUCCESS) then
        print *, "Error in commiting MPI type for 3D vector: ", ierr
        call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminates all MPI processes
    end if
    
  end subroutine mpi_vector3d

  !Subroutine to create an MPI datatype for a 3D point
  subroutine mpi_point3d(mpi_point3d)
    type(MPI_datatype), intent(out) :: mpi_point3d
    integer(int64) :: ierr

    call MPI_Type_contiguos(3, MPI_REAL8, mpi_point3d, ierr)

    !Check for errors in MPI_Type_contiguos
    if (ierr /= MPI_SUCCESS) then
       print*, "Error in creating MPI type for 3D point:", ierr
       call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminates all MPI processes
    end if

    call MPI_Type_commit(mpi_point3d, ierr)

    !Check for errors in MPI_Type_commit
    if (ierr /= MPI_SUCCESS) then
        print *, "Error in commiting MPI type for 3D point: ", ierr
        call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminates all MPI processes
    end if

  end subroutine mpi_point3d
  
  !Subroutine to create an MPI datatype for particle
  subroutine mpi_particle(mpi_particle)
    
  end subroutine mpi_particle
  
  !Subroutine to deallocate MPI datatypes
  subroutine deallocate_mpi_types
    
    call MPI_Type_free(mpi_vector3d, ierr) !Deallocate MPI 3D vector
    call MPI_Type_free(mpi_point3, ierr)   !Deallocate MPI 3D point
    call MPI_Type_free(mpi_particle, ierr) !Deallocate MPI particle
    
  end subroutine deallocate_mpi_types
    
end module mpi_types
