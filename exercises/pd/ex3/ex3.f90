program ex3
    use mpi
    use barneshut
    use geometry
    use particle
    implicit none

    real(8) :: dt, t_end, t, dt_out, t_out, exec_time
    integer :: start, end, count_rate, i, n, rank, size, ierr
    integer :: particles_per_proc, start_idx, end_idx
    type(particle3d), dimension(:), allocatable :: particles
    real(8), dimension(:), allocatable :: buffer_send, buffer_recv, gathered_particles
    type(cell), pointer :: head, temp_cell
    integer :: particle_size

    ! Inicializar MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

    ! Leer las entradas
    if (rank == 0) then
        print*, "Timestep value:"
        read*, dt
        print*, "Output time value:"
        read*, dt_out
        print*, "Value of the final time:"
        read*, t_end
        print*, "Number of particles:"
        read*, n
    end if

    ! Distribuir n entre los procesos
    call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    particles_per_proc = n / size

    ! Determinar el tamaño de cada partícula como arreglo lineal
    particle_size = 10  ! 1 masa + 3 posición + 3 velocidad + 3 aceleración

    ! Preparar buffers
    if (rank == 0) then
        allocate(buffer_send(n * particle_size))
        allocate(particles(n))

        open(unit=4, file="initial_conditions.txt", status="old")  
        do i = 1, n
            read(4,*) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
                     particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
        end do
        close(unit=4)

        ! Copiar datos al buffer de envío
        do i = 1, n
            buffer_send((i-1)*particle_size+1:i*particle_size) = [particles(i)%m, &
                                                                  real(particles(i)%p%x, kind=8), &
                                                                  real(particles(i)%p%y, kind=8), &
                                                                  real(particles(i)%p%z, kind=8), &
                                                                  real(particles(i)%v%x, kind=8), &
                                                                  real(particles(i)%v%y, kind=8), &
                                                                  real(particles(i)%v%z, kind=8), &
                                                                  0.0_8, 0.0_8, 0.0_8]
        end do
    end if

    allocate(buffer_recv(particles_per_proc * particle_size))
    ! Solo asigna 'particles' si no está asignado aún, de lo contrario, desasigna antes
    if (.not. allocated(particles)) then
    	allocate(particles(particles_per_proc))  ! Asignar solo si no está asignado
    else
    	deallocate(particles)  ! Desasignar primero si ya está asignado
    	allocate(particles(particles_per_proc))  ! Volver a asignar
    end if

    
    ! Allocate a separate buffer for gathered data
    if (rank == 0) then
        allocate(gathered_particles(n * particle_size))
    end if

    ! Distribuir datos con MPI_Scatter
    call MPI_Scatter(buffer_send, particles_per_proc * particle_size, MPI_REAL8, &
                     buffer_recv, particles_per_proc * particle_size, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    ! Reconstruir partículas locales
    do i = 1, particles_per_proc
        particles(i)%m = buffer_recv((i-1)*particle_size+1)
        particles(i)%p%x = buffer_recv((i-1)*particle_size+2)
        particles(i)%p%y = buffer_recv((i-1)*particle_size+3)
        particles(i)%p%z = buffer_recv((i-1)*particle_size+4)
        particles(i)%v%x = buffer_recv((i-1)*particle_size+5)
        particles(i)%v%y = buffer_recv((i-1)*particle_size+6)
        particles(i)%v%z = buffer_recv((i-1)*particle_size+7)
        particles(i)%a = vector3d(0.0, 0.0, 0.0)
    end do

    ! Inicializar el árbol y los cálculos
    allocate(head)
    call Calculate_Ranges(head, particles)
    head%type = 0
    call Nullify_Pointers(head)

    ! Crear el árbol inicial
    do i = 1, n 
        call Find_Cell(head, temp_cell, particles(i))
        call Place_Cell(temp_cell, particles(i), i)
    end do

    call Borrar_empty_leaves(head)
    call Calculate_masses(head, particles)

    ! Calcular aceleraciones iniciales
    do i = 1, n
        particles(i)%a = vector3d(0.0, 0.0, 0.0)
    end do
    call Calculate_forces(head, particles)

    ! Guardar resultados:
    open(unit=5, file="output.dat", status='replace')

    ! Bucle principal
    t_out = 0.0
    t = 0.0
    do while (t <= t_end)
        do i = 1, n 
            particles(i)%v = vsum(particles(i)%v, mulvr(dt / 2, particles(i)%a))
            particles(i)%p = sumvp(particles(i)%p, mulvr(dt, particles(i)%v))
        end do

        ! Borrar y reinicializar el árbol
        call Borrar_tree(head)
        call Calculate_Ranges(head, particles)
        head%type = 0
        call Nullify_Pointers(head)

        do i = 1, n 
            call Find_Cell(head, temp_cell, particles(i))
            call Place_Cell(temp_cell, particles(i), i)
        end do

        call Borrar_empty_leaves(head)
        call Calculate_masses(head, particles)

        do i = 1, n
            particles(i)%a = vector3d(0.0, 0.0, 0.0)
        end do
        call Calculate_forces(head, particles)

        do i = 1, n 
            particles(i)%v = vsum(particles(i)%v, mulvr(dt / 2, particles(i)%a))
        end do

        t_out = t_out + dt

        ! Guardar resultados
        if (t_out >= dt_out) then
            call MPI_Gather(buffer_recv, particles_per_proc * particle_size, MPI_REAL8, &
                            gathered_particles, particles_per_proc * particle_size, MPI_REAL8, &
                            0, MPI_COMM_WORLD, ierr)
            if (rank == 0) then
                do i = 1, n 
                    write(5,*) t, i, gathered_particles((i-1)*particle_size+1)
                end do
            end if
            t_out = 0.0
        end if
        t = t + dt
    end do

    close(unit=5)

    ! Calcular tiempo de ejecución
    call system_clock(count=end)
    exec_time = (end - start) / count_rate

    if (rank == 0) then
        print*, "Execution time = ", exec_time, " s"
    end if

    ! Finalizar MPI
    call MPI_FINALIZE(ierr)

end program ex3
