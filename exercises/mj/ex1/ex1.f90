PROGRAM leapfrog
    use geometry
    use particle
    implicit none

    INTEGER :: i, j, k, n, num_steps
    REAL(kind=8) :: dt, t_end, dt_out, t_out
    REAL(kind=8) :: r2, r3
    type(vector3d) :: rji
    type(particle3d), allocatable :: particles(:) !particulas
    type(vector3d), allocatable :: a(:) !aceleracion

    ! Leer parámetros de tiempo y número de partículas
    PRINT*, "Introduce dt"
    READ(*,*) dt
    PRINT*, "Introduce dtout"
    READ(*,*) dt_out
    PRINT*, "Introduce tend"
    READ(*,*) t_end
    PRINT*, "Introduce numero de particulas"
    READ(*,*) n

    ! Determinar el número de pasos de tiempo
    num_steps = INT(t_end / dt)

    ! Asignación de memoria para partículas
    ALLOCATE(particles(n))
    ALLOCATE(a(n))

    ! Leer los datos de masa, posición y velocidad inicial de cada partícula
    DO i = 1, n
        READ(*,*) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
                   particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
    END DO

    ! Inicializar aceleraciones a cero
    DO i = 1, n
        particles(i)%a = vector3d(0.0, 0.0, 0.0)
    END DO

    ! Cálculo inicial de la aceleración
    DO i = 1, n
        DO j = i + 1, n
            rji = vector3d(particles(j)%p%x - particles(i)%p%x, &
                           particles(j)%p%y - particles(i)%p%y, &
                           particles(j)%p%z - particles(i)%p%z)
            r2 = rji%x**2 + rji%y**2 + rji%z**2
            r3 = r2 * sqrt(r2)

            particles(i)%a = particles(i)%a + mulvr(rji, particles(j)%m / r3)
            particles(j)%a = particles(j)%a - mulvr(rji, particles(i)%m / r3)
        END DO
    END DO

    ! Integración de tiempo y actualización de posiciones/velocidades
    t_out = 0.0
    DO i = 1, num_steps
        ! Actualizar velocidades
        DO j = 1, n
            particles(j)%v = particles(j)%v + mulvr(particles(j)%a, dt / 2.0)
        END DO

        ! Actualizar posiciones
        DO j = 1, n
            particles(j)%p = particles(j)%p + mulvr(particles(j)%v, dt)
        END DO

        ! Recalcular aceleración
        DO j = 1, n
            particles(j)%a = vector3d(0.0, 0.0, 0.0)  ! Inicializar a cero en cada paso
        END DO

        DO j = 1, n
            DO k = j + 1, n
                rji = vector3d(particles(k)%p%x - particles(j)%p%x, &
                               particles(k)%p%y - particles(j)%p%y, &
                               particles(k)%p%z - particles(j)%p%z)
                r2 = rji%x**2 + rji%y**2 + rji%z**2
                r3 = r2 * sqrt(r2)

                particles(j)%a = particles(j)%a + mulvr(rji, particles(k)%m / r3)
                particles(k)%a = particles(k)%a - mulvr(rji, particles(j)%m / r3)
            END DO
        END DO

        ! Actualizar velocidades 
        DO j = 1, n
            particles(j)%v = particles(j)%v + mulvr(particles(j)%a, dt / 2.0)
        END DO

        ! Imprimir posiciones
        t_out = t_out + dt
        IF (t_out >= dt_out) THEN
            DO j = 1, n
                PRINT*, particles(j)%p%x, particles(j)%p%y, particles(j)%p%z
            END DO
            t_out = 0.0
        END IF
    END DO

END PROGRAM leapfrog
