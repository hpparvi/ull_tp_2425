PROGRAM ex3
  USE geometry
  USE particle
  USE bh
  USE mpi
  IMPLICIT NONE

  INTEGER :: start_count, end_count, count_rate, n, n_local, i, ierr, rank, size
  REAL :: dt, t_end, t, dt_out, t_out, execution_time
  REAL(kind=kind(1.0d0)) :: theta
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: par, local_par
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a, local_a
  TYPE(cell), POINTER :: head_ptr

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

  ! Initialization and input
  IF (rank == 0) THEN
    print*, "Enter value of the timestep, dt: "
    read*, dt
    print*, "Enter value of the output time, dt_out: "
    read*, dt_out
    print*, "Enter value of the final time, t_end: "
    read*, t_end
    print*, "Enter value of the number of particles, n: "
    read*, n
    print*, "Enter value of theta: "
    read*, theta
  END IF

  CALL MPI_BCAST(dt, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(dt_out, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(t_end, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(theta, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

  n_local = n / size
  IF (MOD(n, size) /= 0 .AND. rank == 0) THEN
    print*, "Error: Number of particles not divisible by number of processes."
    CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  END IF

  ALLOCATE(local_par(n_local))
  ALLOCATE(local_a(n_local))
  ALLOCATE(a(n))

  IF (rank == 0) THEN
    ALLOCATE(par(n))
    open(unit=10, file="IC.txt", status='old')
    DO i = 1, n
      read(10, *) par(i)%m, par(i)%p, par(i)%v
    END DO
    close(unit=10)
  END IF

  CALL MPI_SCATTER(par, n_local, MPI_TYPE_PARTICLE3D, & 
                   local_par, n_local, MPI_TYPE_PARTICLE3D, & 
                   0, MPI_COMM_WORLD, ierr)

  local_a = vector3d(0.0, 0.0, 0.0)

  CALL system_clock(count=start_count, count_rate=count_rate)

  t = 0.0
  t_out = 0.0

  ! Main loop
  DO WHILE (t <= t_end)

    ! Half-step velocity update
    DO i = 1, n_local
      local_par(i)%v = local_par(i)%v + mulrv(REAL(dt / 2, KIND=8), local_a(i))
    END DO

    ! Position update
    DO i = 1, n_local
      local_par(i)%p = local_par(i)%p + mulrv(REAL(dt, KIND=8), local_par(i)%v)
    END DO

    ! Rebuild tree and calculate local forces
    ALLOCATE(head)
    CALL Calculate_ranges(head, local_par)
    head%type = 0
    CALL Nullify_Pointers(head)

    DO i = 1, n_local
      CALL Find_Cell(head, temp_cell, local_par(i))
      CALL Place_Cell(temp_cell, local_par(i), i)
    END DO

    CALL Borrar_empty_leaves(head)
    CALL Calculate_masses(head, local_par)

    local_a = vector3d(0.0, 0.0, 0.0)

    ! Llamada a la subrutina con theta
    CALL Calculate_forces(head, local_par, n_local, theta)

    ! Gather forces from all processes
    CALL MPI_ALLGATHER(local_a, n_local, MPI_TYPE_VECTOR3D, & 
                       a, n_local, MPI_TYPE_VECTOR3D, & 
                       MPI_COMM_WORLD, ierr)

    ! Half-step velocity update with global forces
    DO i = 1, n_local
      local_par(i)%v = local_par(i)%v + mulrv(REAL(dt / 2, KIND=8), a((rank * n_local) + i))
    END DO

    t = t + dt
    t_out = t_out + dt

    ! Output results
    IF (t_out >= dt_out .AND. rank == 0) THEN
      open(unit=11, file="output.dat", status='replace')
      DO i = 1, n
        WRITE(11, *) i, par(i)%p
      END DO
      close(unit=11)
      t_out = 0.0
    END IF

  END DO

  CALL system_clock(count=end_count)
  execution_time = REAL(end_count - start_count) / REAL(count_rate)

  IF (rank == 0) THEN
    print*, "Total execution time:", execution_time, "seconds"
  END IF

  CALL MPI_FINALIZE(ierr)

END PROGRAM ex3
