program initialconditions

    implicit none

    integer :: i, N
    real(8) :: mass, rx, ry, rz, vx, vy, vz, theta, phi, velocidad_magnitud, r
    real(8), dimension(:), allocatable :: seed

    open(unit=4, file='initial_conditions.txt', status='replace')   !Aquí guardaremos la información

    print*, "Number of bodies?"
    read*, N

    print*, "Magnitude of velocity?"
    read*, velocidad_magnitud

    mass = 1.0d0 / N

    do i = 1, N
        ! Generar phi entre 0 y 2*pi
        call random_number(phi)
        phi = phi * 2.0d0 * 3.141592653589793d0

        ! Generar cos(theta) uniformemente entre -1 y 1
        call random_number(theta)
        theta = acos(1.0d0 - 2.0d0 * theta)

        ! Convertir a coordenadas cartesianas para posición
        rx = sin(theta) * cos(phi)
        ry = sin(theta) * sin(phi)
        rz = cos(theta)

        ! Calcular la distancia radial desde el centro (origen)
        r = sqrt(rx**2 + ry**2 + rz**2)

        ! Normalizar las componentes de la posición y escalar por la magnitud de la velocidad
        vx = -velocidad_magnitud * rx / r
        vy = -velocidad_magnitud * ry / r
        vz = -velocidad_magnitud * rz / r

        ! Escribir las condiciones iniciales
        write(4, '(F10.3, 3F16.8, 3F16.8)') mass, rx, ry, rz, vx, vy, vz

        print*, 'Velocidad de la partícula ', i, ' = ', sqrt(vx**2 + vy**2 + vz**2)

    end do

    close(unit=4)
    close(unit=3)

end program initialconditions
