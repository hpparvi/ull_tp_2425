program ex2
    !$use omp_lib
    use barneshut
    use geometry
    use particle
    implicit none

    real(8):: dt, t_end, t, dt_out, t_out, exec_time, r2,r3
    type(particle3d), dimension(:), allocatable:: particles 
    integer:: start,end, count_rate,i,n
    type(cell), pointer :: head, temp_cell
  

 ! Read inputs
    print*, "Timestep value:"
    read*, dt
    print*, "Output time value:"
    read*, dt_out
    print*, "Value of the final time:"
    read*, t_end
    print*, "Number of particles:"
    read*, n

    allocate(particles(n))
   
    ! We open the file where initial conditions are written
    open(unit=4, file="initial_conditions.txt", status="old")  

    do i=1,n 
        read(4,*) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
        particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
    end do
    
    close(unit=4)

    
 
    ! Calculate the starting time
    call system_clock(count=start, count_rate=count_rate)

    !! Inicialización head node
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(head)
    call Calculate_Ranges(head,particles)
    head%type = 0
    call Nullify_Pointers(head)

    !! Creación del árbol inicial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,n 
        call Find_Cell(head,temp_cell,particles(i))
        call Place_Cell(temp_cell, particles(i),i)
    end do

    call Borrar_empty_leaves(head)
    call Calculate_masses(head, particles)

    
    !! Calcular aceleraciones iniciales
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel shared(particles, head, n, dt)
    do i = 1, n
        particles(i)%a = vector3d(0.0, 0.0, 0.0)
    end do
    call Calculate_forces(head, particles)
    !$omp end parallel

    !!Save results:
    open(unit=5, file="output.dat", status='replace')

    !! Bucle principal
!!!!!!!!!!!!!!!!!!

    t_out = 0.0
    t = 0.0
    do while (t <= t_end)
        !$omp parallel shared(particles, head, n, dt)
        do i=1,n 
            particles(i)%v = vsum(particles(i)%v, mulvr(dt / 2, particles(i)%a))
            particles(i)%p = sumvp(particles(i)%p , mulvr(dt, particles(i)%v))
        end do
        !! Las posiciones han cambiado, por lo que tenemos que borrar
!! y reinicializar el árbol
        
        !$omp single
        call Borrar_tree(head)
        call Calculate_Ranges(head,particles)
        head%type = 0
        call Nullify_Pointers(head)

        do i=1,n 
            call Find_Cell(head,temp_cell, particles(i))
            call Place_Cell(temp_cell, particles(i), i)
        end do

        call Borrar_empty_leaves(head)
        call Calculate_masses(head, particles)
        !$omp end single
        do i = 1, n
            particles(i)%a = vector3d(0.0, 0.0, 0.0)
        end do
        call Calculate_forces(head, particles)

        do i=1,n 
            particles(i)%v = vsum(particles(i)%v, mulvr(dt / 2, particles(i)%a))
        end do
        !$omp end parallel

        t_out = t_out + dt

        if (t_out >= dt_out) then
            do i = 1,n 
                write(5,*) t, i, particles(i)%p
            end do
            t_out = 0.0
        end if
        t = t + dt
    end do

    close(unit=5)

    !Calculate final time
    call system_clock(count=end)

    !Execution time:
    exec_time = (end - start) / count_rate

    print*, "Execution time = ", exec_time," s"

end program ex2
