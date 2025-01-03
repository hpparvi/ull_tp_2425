module particles
  ! module that contains
  !   - data type: particle, with position, velocity and mass

  use geometry

  implicit none

  type :: particle
     type(point3d)     :: p  ! position
     type(vector3d)    :: v  ! velocity
     type(vector3d)    :: a  ! acceleration
     real  :: m              ! mass
  end type particle

end module particles



module mpi_particles
  use geometry
  use mpi_f08

  implicit none

contains
 subroutine create_mpi_datatypes(MPI_point3d, MPI_vector3d, MPI_particle)
  INTEGER :: ierr
  INTEGER :: lengths(0:2)
  INTEGER(KIND=MPI_ADDRESS_KIND) :: displacements(0:2)
  TYPE(MPI_Datatype) :: types(0:2)
  type(point3d) :: dummy_point
  type(vector3d) :: dummy_vector
  ! for mpi particle type
  INTEGER :: lengths_p(0:3)
  INTEGER(KIND=MPI_ADDRESS_KIND) :: displacements_p(0:3)
  TYPE(MPI_Datatype) :: types_p(0:3)

  type(MPI_Datatype), intent(out) ::  MPI_particle, MPI_vector3d, MPI_point3d

  ! Create point3d datatype for mpi
  lengths = [1, 1, 1]
  displacements = [0 * SIZEOF(dummy_point%x), &
       & sizeof(dummy_point%x), &
       & sizeof(dummy_point%x) + sizeof(dummy_point%y)]
  types = [MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION]
  CALL MPI_Type_create_struct(3, lengths, displacements, types, MPI_point3d, ierr)
  if (ierr.ne.0) then
     print *, 'Error creating mpi type ', ierr
     call MPI_Abort(MPI_COMM_WORLD, ierr)
  end if
  CALL MPI_Type_commit(MPI_point3d, ierr)
  if (ierr.ne.0) then
     print *, 'Error creating mpi type ', ierr
     call MPI_Abort(MPI_COMM_WORLD, ierr)
  end if
  ! Create vector3d datatype for mpi
  displacements = [0 * SIZEOF(dummy_vector%x), &
       & sizeof(dummy_vector%x), &
       & sizeof(dummy_vector%x) + sizeof(dummy_vector%y)]
  CALL MPI_Type_create_struct(3, lengths, displacements, types, MPI_vector3d, ierr)
  if (ierr.ne.0) then
     print *, 'Error creating mpi type ', ierr
     call MPI_Abort(MPI_COMM_WORLD, ierr)
  end if
  CALL MPI_Type_commit(MPI_vector3d, ierr)
  if (ierr.ne.0) then
     print *, 'Error creating mpi type ', ierr
     call MPI_Abort(MPI_COMM_WORLD, ierr)
  end if

  ! Create particle datatype for mpi
  lengths_p = [ 1, 1, 1, 1]
  displacements_p = [ 0*sizeof(dummy_point), &
       & sizeof(dummy_point), &
       & sizeof(dummy_point) + sizeof(dummy_vector), &
       & sizeof(dummy_point) + sizeof(dummy_vector) + sizeof(dummy_vector)]
   types_p = [ mpi_point3d, &
       & mpi_vector3d, &
       & mpi_vector3d, &
       & mpi_DOUBLE_PRECISION]

   CALL MPI_Type_create_struct(4, lengths_p, displacements_p, types_p, MPI_particle, ierr)
   if (ierr.ne.0) then
      print *, 'Error creating mpi type ', ierr
      call MPI_Abort(MPI_COMM_WORLD, ierr)
   end if
   CALL MPI_Type_commit(MPI_particle, ierr)
   if (ierr.ne.0) then
      print *, 'Error creating mpi type ', ierr
      call MPI_Abort(MPI_COMM_WORLD, ierr)
   end if

 end subroutine create_mpi_datatypes
end module mpi_particles
