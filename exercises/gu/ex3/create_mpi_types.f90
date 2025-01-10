module create_mpi_types
  use, intrinsic :: iso_fortran_env
  use mpi_f08
  use geometry
  use particle
  use tree
  implicit none

contains
  subroutine build_derived_vector(mesg_mpi_vector)
    type(vector3d) :: dummy_vector
    integer, parameter :: num_blocks = 3 !things inside vector3d
    integer :: ierr, i
    integer :: lengths(0:2)
    integer(kind=MPI_Address_Kind) :: displacements(0:2)
    type(MPI_Datatype) :: types(0:2)
    type(MPI_Datatype), intent(out) :: mesg_mpi_vector

    lengths = [1, 1, 1]
    displacements = [0*sizeof(dummy_vector%xx), &
         &sizeof(dummy_vector%xx), &
    &sizeof(dummy_vector%xx) + sizeof(dummy_vector%yy)]
    types = [MPI_Double_Precision, MPI_Double_Precision, MPI_Double_Precision]
    
    call MPI_Type_create_struct(3, lengths, displacements, types, mesg_mpi_vector, ierr)
    if (ierr.ne.0) then
       print *, 'Error in type create: ', ierr
       call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
    end if

    call MPI_TYPE_COMMIT(mesg_mpi_vector, ierr)
    if (ierr .ne. 0) then
       print *, 'Error in type commit: ', ierr
       call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
    end if

    return

  end subroutine build_derived_vector

  subroutine build_derived_point(mesg_mpi_point)
    type(point3d) :: dummy_point
    integer, parameter :: num_blocks = 3 !things inside vector3d
    integer :: ierr, i
    integer :: lengths(0:2)
    integer(kind=MPI_Address_Kind) :: displacements(0:2)
    type(MPI_Datatype) :: types(0:2)
    type(MPI_Datatype), intent(out) :: mesg_mpi_point

    lengths = [1, 1, 1]
    displacements = [0*sizeof(dummy_point%xx), &
         &sizeof(dummy_point%xx), &
    &sizeof(dummy_point%xx) + sizeof(dummy_point%yy)]
    types = [MPI_Double_Precision, MPI_Double_Precision, MPI_Double_Precision]
    
    call MPI_Type_create_struct(3, lengths, displacements, types, mesg_mpi_point, ierr)
    if (ierr.ne.0) then
       print *, 'Error in type create: ', ierr
       call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
    end if

    call MPI_TYPE_COMMIT(mesg_mpi_point, ierr)
    if (ierr .ne. 0) then
       print *, 'Error in type commit: ', ierr
       call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
    end if

    return

  end subroutine build_derived_point

  subroutine build_derived_particle(mesg_mpi_particle)
    type(particle3d) :: dummy_particle
    type(MPI_Datatype) :: MPI_Vector3d, MPI_point3d
    integer, parameter :: num_blocks = 3 !things inside vector3d
    integer :: ierr, i
    integer :: lengths(0:2)
    integer(kind=MPI_Address_Kind) :: displacements(0:2)
    type(MPI_Datatype) :: types(0:2)
    type(MPI_Datatype), intent(out) :: mesg_mpi_particle

    call build_derived_vector(MPI_vector3d)
    call build_derived_point(MPI_point3d)

    lengths = [1, 1, 1]
    displacements = [0*sizeof(dummy_particle%m), &
         &sizeof(dummy_particle%m), &
    &sizeof(dummy_particle%m) + sizeof(dummy_particle%p)]
    types = [MPI_Double_Precision, MPI_point3d, MPI_vector3d]
    
    call MPI_Type_create_struct(3, lengths, displacements, types, mesg_mpi_particle, ierr)
    if (ierr.ne.0) then
       print *, 'Error in type create: ', ierr
       call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
    end if

    call MPI_TYPE_COMMIT(mesg_mpi_particle, ierr)
    if (ierr .ne. 0) then
       print *, 'Error in type commit: ', ierr
       call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
    end if

    return

  end subroutine build_derived_particle

  subroutine build_derived_RANGE(mesg_mpi_range)
    type(range) :: dummy_range
    integer, parameter :: num_blocks = 2 !things inside range
    integer :: ierr, i
    integer :: lengths(0:1)
    integer(kind=MPI_Address_Kind) :: displacements(0:1)
    type(MPI_Datatype) :: types(0:1)
    type(MPI_Datatype), intent(out) :: mesg_mpi_range

    lengths = [3, 3]
    displacements = [0*sizeof(dummy_range%rmin), sizeof(dummy_range%rmin)]
    types = [MPI_Double_Precision, MPI_Double_Precision]
    
    call MPI_Type_create_struct(num_blocks, lengths, displacements, types, mesg_mpi_range, ierr)
    if (ierr.ne.0) then
       print *, 'Error in type create: ', ierr
       call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
    end if

    call MPI_TYPE_COMMIT(mesg_mpi_range, ierr)
    if (ierr .ne. 0) then
       print *, 'Error in type commit: ', ierr
       call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
    end if

    return

  end subroutine build_derived_range

end module create_mpi_types
