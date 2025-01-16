module mpi_env
  use geometry  !Importing the geometry module
  use particle  !Importing the particle module
  use barnes    !Importing the barnes module, where the Barnes-Hut algorithm is defined
  use calcs     !Importing the calcs module, where subroutines to update properties of the particles are defined
  use data      !Importing the data module, where subroutines to read and save data are defined
  use mpi_types !Importing the mpi_types module, where subroutines to use MPI datatypes are defined
  use iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  use mpi_f08   !Importing the mpi library to use MPI 
  implicit none

contains
  
  !Subroutine to initialize the MPI environment
  subroutine Initialize_MPI(comsize, rank, ierr)
    integer, intent(inout) :: comsize, rank, ierr !Total number of processes, processes IDs and variable for error handling in MPI operations
    
    call MPI_INIT(ierr) !Initialize the MPI execution environment
    call MPI_COMM_SIZE(MPI_COMM_WORLD, comsize, ierr) !Get the number of MPI processes
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)    !Get the rank of the calling MPI process
    
  end subroutine Initialize_MPI

  !Subroutine to send an integer variable to all processes
  subroutine send_integer(i, ierr)
    integer(int64) :: i !Integer variable
    integer, intent(inout) :: ierr !Integer variable for error handling in MPI operations

    !Send the real variable
    call MPI_Bcast(i, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)

    !Check for errors in sending the integer variable
    if (ierr .NE. MPI_SUCCESS) then
        print *, "Error in broadcasting ", i, " to all processes: ", ierr !Print error message
        call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
     end if

   end subroutine send_integer

  !Subroutine to send a real variable to all processes
  subroutine send_real(r, ierr)
    real(real64) :: r !Real variable
    integer, intent(inout) :: ierr !Integer variable for error handling in MPI operations

    !Send the real variable
    call MPI_Bcast(r, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    !Check for errors in sending the real variable
    if (ierr .NE. MPI_SUCCESS) then
        print *, "Error in broadcasting ", r, " to all processes: ", ierr !Print error message
        call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
     end if

   end subroutine send_real

  !Subroutine to send simulation parameters to all processes
  subroutine send_data(n, dt, dt_out, t_end, theta, ierr)
    integer(int64), intent(in) :: n !Number of bodies
    real(real64), intent(in) :: dt, dt_out, t_end, theta !Time step, output time step and final time and parameter that determines the accuracy of the simulation
    integer, intent(inout) :: ierr !Integer variable for error handling in MPI operations

    call send_integer(n, ierr)   !Send the number of bodies
    call send_real(dt, ierr)     !Send the time step
    call send_real(dt_out, ierr) !Send the output final time
    call send_real(t_end, ierr)  !Send the final time
    call send_real(theta, ierr)  !Send theta
    
  end subroutine send_data

  !Subroutine to send particle data to all processes
  subroutine send_particles(p, MPI_PARTICLE3D, ierr)
    type(particle3d), intent(in) :: p(:) !Particles
    type(MPI_Datatype), intent(in) :: MPI_PARTICLE3D !MPI 3D particle
    integer, intent(inout) :: ierr !Integer variable for error handling in MPI operations

    !Send particle data
    call MPI_Bcast(p, size(p), MPI_PARTICLE3D, 0, MPI_COMM_WORLD, ierr)

    !Check for errors in sending particle data
    if (ierr .NE. MPI_SUCCESS) then
        print *, "Error in broadcasting the particles to all processes: ", ierr !Print error message
        call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
    end if
    
  end subroutine send_particles
  
  !Subroutine to determine how particles are distributed among processes
  subroutine particles_nodes(n, comsize, rank, ierr, i_start, i_end, n_local, n_extra, n_nodes, displacements)
    integer(int64), intent(in) :: n !Number of bodies
    integer(int64) :: i !Loop indexing variable
    integer, intent(inout) :: comsize, rank, ierr !Total number of processes, processes IDs and variable for error handling in MPI operations
    integer, intent(inout) :: i_start, i_end, n_local, n_extra !First and last index of the range of particles to process, number of particles per process and extra particles
    integer, allocatable, intent(inout) :: n_nodes(:), displacements(:) !Number of particles per process and displacement array


    !Calculate the number of particles per process
    n_local = n/comsize
    n_extra = mod(n, comsize)

    !Determine the start and end indices for each process
    if (rank < n_extra) then
       i_start = rank * (n_local + 1) + 1
       i_end = i_start + n_local
    else
       i_start = rank * n_local + n_extra + 1
       i_end = i_start + n_local - 1
    end if

    !Allocate memory for the number of particles per process and displacement array
    allocate(n_nodes(comsize), displacements(comsize))

    !Assign the number of particles per process
    do i = 1, comsize
       n_nodes(i) = n_local
    end do

    !Distribute extra particles among the first few processes
    do i = 1, n_extra
       n_nodes(i) = n_nodes(i) + 1
    end do

    !Compute displacements for gathering particle data from all processes
    displacements(1) = 0
    do i = 2, comsize
       displacements(i) = displacements(i-1) + n_nodes(i-1)
    end do

  end subroutine particles_nodes
  
  !Subroutine to gather particle data from all processes
  subroutine gather_particles(p, i_start, i_end, n_nodes, displacements, MPI_PARTICLE3D, rank, ierr)
    type(particle3d), intent(inout) :: p(:) !Particles
    type(MPI_Datatype), intent(in) :: MPI_PARTICLE3D !MPI 3D particle
    integer, intent(in) :: i_start, i_end, rank !First and last index of the range of particles to process and processor's ID
    integer, intent(inout) :: ierr, n_nodes(:), displacements(:) !Variable for error handling in MPI operations, number of particles per process and displacement array

    !Perform an all-to-all gather of particle data across all processes
    CALL MPI_Allgatherv(p(i_start:i_end), n_nodes(rank+1), MPI_PARTICLE3D, p, n_nodes, &
         displacements, MPI_PARTICLE3D, MPI_COMM_WORLD, ierr)

    !Check for errors in gathering particle data
    if (ierr .NE. MPI_SUCCESS) then
        print *, "Rank ", rank, ": Error in gathering the particles across all processes: ", ierr !Print error message
        call MPI_Abort(MPI_COMM_WORLD, ierr) !Terminate all MPI processes
    end if

  end subroutine gather_particles

  !Subroutine to finalize the MPI environment
  subroutine Finalize_MPI(ierr)
    integer, intent(inout) :: ierr !Integer variable for error handling in MPI operations

    !Finalize the MPI environment
    call MPI_FINALIZE(ierr)
    
  end subroutine Finalize_MPI
  
end module mpi_env

