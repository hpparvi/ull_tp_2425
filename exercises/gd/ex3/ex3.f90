PROGRAM ex3

  USE mpi_f08
  USE geometry
  USE ex3_deriv_types
  USE barnes_hut
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE      
  
  REAL(REAL64) :: dt, t_end, t, dt_out, t_out    
  character(len=25) :: input, output, orbit
  INTEGER :: comsize, ierr, workers, counts, counts_finale, i_start, i_end
  INTEGER, ALLOCATABLE :: sendcounts(:), displs(:)
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: pt_workers
  ! To test time of execution
  INTEGER :: count_begin, count_end, count_rate

  !! MPI variables
  TYPE(MPI_Datatype) :: mpi_point3d_type, mpi_vector3d_type, mpi_particle3d_type
  INTEGER, DIMENSION(3) :: block_lengths = [1, 1, 1]
  INTEGER(MPI_Address_kind) :: displacements(3)
  TYPE(MPI_Datatype), DIMENSION(3) :: old_types
  INTEGER(MPI_Address_kind) :: base_address, addr_p, addr_v, addr_m
  TYPE(particle3d) :: particle

  CALL MPI_init(ierr)
  CALL MPI_comm_size(MPI_COMM_WORLD, comsize, ierr)
  CALL MPI_comm_rank(MPI_COMM_WORLD, rank, ierr)
  workers = comsize - 1

  ALLOCATE(sendcounts(comsize), displs(comsize))
  
!!  Data reading
  
  input = "initial_conditions.dat"
  output = "output.dat"  
 
!!  We turn node 0 into the root and read the information from the file.
  root = 0
  IF (rank==root) THEN
    OPEN (1,file = input, status = 'old', action = 'read')
    READ(1,*) dt
    READ(1,*) dt_out
    READ(1,*) t_end
    READ(1,*) n
    ALLOCATE(pt(n))
    DO i = 1, n
      READ(1,*) pt(i)%m, pt(i)%p, pt(i)%v
    END DO
    CLOSE (UNIT=1) 
  
  END IF
  

!!  Send the information from the root to the rest of the nodes (workers).
  
  call MPI_Bcast(dt, 1, MPI_Double_precision, root, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(dt_out, 1, MPI_Double_precision, root, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(t_end, 1, MPI_Double_precision, root, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(n, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  ALLOCATE(a(n))  
  
  
!!  Create MPI types for point3d and vector3d.
  old_types = [MPI_Double_precision, MPI_Double_precision, MPI_Double_precision]  
  
!!  Define the point3d type.
  CALL MPI_Get_Address(particle%p, base_address, ierr)
  CALL CheckMPIError(ierr, "MPI_Get_Address 1")
  CALL MPI_Get_Address(particle%p%x, displacements(1), ierr)
  CALL CheckMPIError(ierr, "MPI_Get_Address 2")
  CALL MPI_Get_Address(particle%p%y, displacements(2), ierr)
  CALL CheckMPIError(ierr, "MPI_Get_Address 3")
  CALL MPI_Get_Address(particle%p%z, displacements(3), ierr)
  CALL CheckMPIError(ierr, "MPI_Get_Address 4")
  displacements = displacements - base_address
  CALL MPI_Type_create_struct(3, block_lengths, displacements, old_types, mpi_point3d_type, ierr)
  CALL CheckMPIError(ierr, "MPI_Type_create_struct 1")
  CALL MPI_Type_commit(MPI_point3d_type, ierr)
  CALL CheckMPIError(ierr, "MPI_Type_commit 1")

!!  Define the vector3d type.
  CALL MPI_Get_Address(particle%v, base_address, ierr)
  CALL CheckMPIError(ierr, "MPI_Get_Address 11")
  CALL MPI_Get_Address(particle%v%x, displacements(1), ierr)
  CALL CheckMPIError(ierr, "MPI_Get_Address 12")
  CALL MPI_Get_Address(particle%v%y, displacements(2), ierr)
  CALL CheckMPIError(ierr, "MPI_Get_Address 13")
  CALL MPI_Get_Address(particle%v%z, displacements(3), ierr)
  CALL CheckMPIError(ierr, "MPI_Get_Address 14")
  displacements = displacements - base_address
  CALL MPI_Type_create_struct(3, block_lengths, displacements, old_types, mpi_vector3d_type, ierr)
  CALL CheckMPIError(ierr, "MPI_Type_create_struct 11")
  CALL MPI_Type_commit(MPI_vector3d_type, ierr)
  CALL CheckMPIError(ierr, "MPI_Type_commit 11")

!!  Create an MPI type for particle3d. It depends on the two previous types.
  old_types = [mpi_point3d_type, mpi_vector3d_type, MPI_Double_precision]

  CALL MPI_Get_Address(particle, base_address, ierr)
  CALL CheckMPIError(ierr, "MPI_Get_Address 21")
  CALL MPI_Get_Address(particle%p, addr_p, ierr)
  CALL CheckMPIError(ierr, "MPI_Get_Address 22")
  CALL MPI_Get_Address(particle%v, addr_v, ierr)
  CALL CheckMPIError(ierr, "MPI_Get_Address 23")
  CALL MPI_Get_Address(particle%m, addr_m, ierr)
  CALL CheckMPIError(ierr, "MPI_Get_Address 24")
  displacements = [addr_p - base_address, addr_v - base_address, addr_m - base_address]

  CALL MPI_Type_create_struct(3, block_lengths, displacements, old_types, mpi_particle3d_type, ierr)
  CALL CheckMPIError(ierr, "MPI_Type_create_struct 21")
  CALL MPI_Type_commit(MPI_particle3d_type, ierr)  
  CALL CheckMPIError(ierr, "MPI_Type_commit 21")
  
!!  Move the particle information from the root to the workers
  IF (rank /= root) ALLOCATE(pt(n))
  CALL MPI_Bcast(pt, n, MPI_particle3d_type, root, MPI_COMM_WORLD, ierr)
  CALL CheckMPIError(ierr, "MPI_Bcast 1")
  
!!  Start measuring time here
  IF (rank == root) CALL SYSTEM_CLOCK(count_begin, count_rate)

  
!!  Initialization of head node
  ALLOCATE(head)
  CALL Calculate_ranges(head)
  head%type = 0
  CALL Nullify_Pointers(head)
  
!!  Create the initial tree
  DO i = 1,n
    CALL Find_Cell(head,temp_cell,pt(i)%p)
    CALL Place_Cell(temp_cell,pt(i)%p,i)
  END DO
  CALL Borrar_empty_leaves(head)
  CALL Calculate_masses(head)

!! Calculate the limits and the size to be sent.
  counts = CEILING(REAL(n)/REAL(workers))
  counts_finale = n - counts * (workers - 1)
  
!! Calculate initial acceleration
!! Calculate the start and end values of the particle loop for each worker node.
!! a message will be show with the number of cores.

  a = vector3d(0,0,0)

  IF (rank == root) THEN
     PRINT'(A,I2,A)', 'Parallelized process with ', comsize, ' cores' 
     i_start = 0
     i_end = 0
  ELSEIF (rank /= workers) THEN
     i_start=1+(rank-1)*counts
     i_end = (rank)*counts
  ELSE
     i_start=1+(rank-1)*counts
     i_end = (rank-1)*counts + counts_finale
  END IF
  
!!  Calculate the parameters needed to send the information  
  DO i = 1, comsize
     IF (i.EQ.1) THEN
       sendcounts(i) = 0
       displs(i) = 0
     ELSEIF (i /= workers+1) THEN
       sendcounts(i) = counts
       displs(i) = counts*(i - 2)
     ELSE
       sendcounts(i) = counts_finale
       displs(i) = counts*(i - 2)
     END IF
  END DO
  ALLOCATE(pt_workers(sendcounts(rank+1)))
  
  
 
  CALL Calculate_forces(head,pt,a,i_start, i_end)
!!  Wait til everyone is done
  !CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  !CALL CheckMPIError(ierr, "MPI_Barrier 1")



!!  Open the file only for the root.
  IF (rank.EQ.root) THEN
    OPEN(3, file = output, status = 'replace', action = 'write')
  END IF
  
  t = 0.0
  t_out = 0.0

!! Main loop
  DO 
    t = t + dt
    t_out = t_out + dt
    
!!  Send each worker its particular group of particles.
    CALL MPI_Scatterv(pt, sendcounts, displs, mpi_particle3d_type, &
         & pt_workers, sendcounts(rank+1), mpi_particle3d_type, root, MPI_COMM_WORLD,ierr)
    CALL CheckMPIError(ierr, "MPI_Scatterv 1")
    pt_workers%v = pt_workers%v + a(i_start:i_end) * (dt/2)
    pt_workers%p = pt_workers%p + pt_workers%v * dt
    
    CALL MPI_Allgatherv(pt_workers, sendcounts(rank+1), mpi_particle3d_type,&
         & pt, sendcounts, displs, mpi_particle3d_type, MPI_COMM_WORLD,ierr)
    CALL CheckMPIError(ierr, "MPI_Allgatherv 1")
    
!!  The positions have changed, so the tree needs to be reset.
    CALL Borrar_tree(head)
    CALL Calculate_ranges(head)
    head%type = 0
    CALL Nullify_Pointers(head)
    DO i = 1,n
      CALL Find_Cell(head,temp_cell,pt(i)%p)
      CALL Place_Cell(temp_cell,pt(i)%p,i)
    END DO
    CALL Borrar_empty_leaves(head)
    CALL Calculate_masses(head)
    a = vector3d(0,0,0)
    CALL Calculate_forces(head,pt,a,i_start, i_end)

!!  Send agaion to each worker its particular group of particles.
    CALL MPI_scatterv(pt, sendcounts, displs, mpi_particle3d_type, &
         & pt_workers, sendcounts(rank+1), mpi_particle3d_type, root, MPI_COMM_WORLD,ierr)
    CALL CheckMPIError(ierr, "MPI_Scatterv 2")
    pt_workers%v = pt_workers%v + a(i_start:i_end) * (dt/2)
    
    CALL MPI_allgatherv(pt_workers, sendcounts(rank+1), mpi_particle3d_type,&
         & pt, sendcounts, displs, mpi_particle3d_type, MPI_COMM_WORLD,ierr)
    CALL CheckMPIError(ierr, "MPI_Allgatherv 2")
         
    !CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    !CALL CheckMPIError(ierr, "MPI_Barrier 2")

    IF ((t_out >= dt_out) .AND. rank == root) THEN
      write (3, '(3F12.4)', advance='no') t
      DO i = 1,n
        IF(i<n) THEN
          write (3, '(6F16.6)', advance='no') pt(i)%p
        ELSE
          write (3, '(6F16.6)') pt(i)%p
        END IF
      END DO
      t_out = 0.0 ! reset of output parameter 
    END IF
!!  Exit conditions
    IF (t>= t_end) THEN
      EXIT
    END IF
  END DO
  
!!  End measuring time here
  IF (rank.EQ.root) THEN
    CLOSE (3)
    CALL SYSTEM_CLOCK(count_end)
    PRINT'(A,F7.3,A)', "The elapsed time is ", &
       REAL(count_end-count_begin) / REAL(count_rate), &
       " seconds"  
  END IF

!!  Deallocate a finish program
  DEALLOCATE(pt)
  DEALLOCATE(a)
  DEALLOCATE(pt_workers)
  DEALLOCATE(sendcounts)
  DEALLOCATE(displs)
  CALL MPI_Type_free(mpi_vector3d_type,ierr)
  CALL CheckMPIError(ierr, "MPI_Type_free 1")
  CALL MPI_Type_free(mpi_point3d_type,ierr)
  CALL CheckMPIError(ierr, "MPI_Type_free 2")
  CALL MPI_Type_free(mpi_particle3d_type,ierr)
  CALL CheckMPIError(ierr, "MPI_Type_free 3")
  PRINT *, "Proceso ", rank, " liberando memoria y finalizando."
  CALL MPI_Finalize(ierr)
  CALL CheckMPIError(ierr, "MPI_Finalize")
  PRINT *, "Proceso ", rank, " finalizó correctamente."
  
CONTAINS
!!  Subroutine to display the MPI error.
  SUBROUTINE CheckMPIError(ierr, func_name)
    USE mpi_f08
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ierr
    CHARACTER(*), INTENT(IN) :: func_name
  
    IF (ierr /= MPI_SUCCESS) THEN
      PRINT *, "Error en la función MPI:", TRIM(func_name), "Código de error:", ierr
      CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
    END IF
  END SUBROUTINE CheckMPIError
END PROGRAM ex3
