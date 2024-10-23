PROGRAM leapfrog
  use particle
  IMPLICIT NONE

  
  INTEGER(int64) :: i,j,k
  INTEGER(int64) :: n                                !número de cuerpos
  REAL(real64) :: dt, t_end, t, dt_out, t_out         !paso temporal, ...
  REAL(real64) :: rs, r2, r3                          !posicion, r al cuadrado y al cubo
  type(vector3d) :: rji                       !posición entre dos cuerpos
  
  type(particle3d), allocatable :: p(:)         !particulas p tienen (%p, v, m)
  type(vector3d), allocatable :: a(:)
  
  character(len=30) :: datos
  character(len=30) :: orbitas
  INTEGER :: io_status  ! Variable porque me daba errores

  datos = 'data_input.dat'
  orbitas = 'data_output.dat'
  OPEN(10, file=datos, status='old', action='read', iostat=io_status)


  ! Leer los parámetros básicos y depurar
  READ(10, *, IOSTAT=io_status) dt
  print *, "El valor leido del salto temporal (dt) es ", dt

  READ(10, *, IOSTAT=io_status) dt_out
  print *, "El valor leido del tiempo de impresion de resultados (dt_out) es", dt_out

  READ(10, *, IOSTAT=io_status) t_end
  print *, "El valor leido del tiempo total de la simulacion (t_end) es ", t_end

  READ(10, *, IOSTAT=io_status) n
  print *, "El valor leido del numero de cuerpos (n) es ", n

  ! Asignar memoria para las partículas
  ALLOCATE(p(n))
  allocate(a(n))

  ! Leer las posiciones, velocidades y masas de las partículas
  DO i = 1, n
     READ(10, *, IOSTAT=io_status) p(i)%p
     print *, "El valor leido de la posicion de la particula", i, "es:", p(i)%p

     READ(10, *, IOSTAT=io_status) p(i)%v
     print *, "El valor leido de la velocidad de la particula", i, "es:", p(i)%v

     READ(10, *, IOSTAT=io_status) p(i)%m
     print *, "El valor leido de la masa de la particula", i, "es:", p(i)%m

  END DO

  CLOSE(10)


  do i = 1, n
     a(i) = vector3d(0.0, 0.0, 0.0)
  end do

  open(11, file = orbitas, status = 'old', action = 'write', iostat = io_status)
  
  DO i = 1, n
     DO j = i+1, n
        rji = distance(p(i)%p, p(j)%p)
        r2 = rji%x**2+rji%y**2+rji%z**2
        r3 = r2 * sqrt(r2)
        a(i) = a(i) + ((p(j)%m * rji) / r3)
        a(j) = a(j) - ((p(i)%m * rji) / r3)
     END DO
  END DO

  t = 0.0                                                  !new
  t_out = 0.0
  DO while (t <= t_end)                                    !changed bc fortran doesnt allow non-entire variables at do-loops
     do i = 1, n                                           !added
     p(i)%v = p(i)%v + (a(i) * (dt/2))
     p(i)%p = p(i)%p + (p(i)%v * dt)
     a(i) = vector3d(0.0, 0.0, 0.0)
     end do
     DO i = 1, n
        DO j = i+1, n
           rji = distance(p(i)%p, p(j)%p)
           r2 = rji%x**2+rji%y**2+rji%z**2
           r3 = r2 * sqrt(r2)
           a(i) = a(i) + ((p(j)%m * rji) / r3)
           a(j) = a(j) - ((p(i)%m * rji) / r3)
        END DO
     END DO
     do i = 1, n                                            !added
      p(i)%v = p(i)%v + (a(i) * (dt/2))
     end do
   
     t_out = t_out + dt
      IF (t_out >= dt_out) THEN
         DO i = 1,n
            write (11, '(3F12.6)') p(i)%p%x, p(i)%p%y, p(i)%p%z
         END DO

         t_out = 0.0
      END IF
      t = t + dt
   END DO

   close(11)
  
END PROGRAM leapfrog
