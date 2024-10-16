PROGRAM leapfrog
  use particl
  IMPLICIT NONE

  
  INTEGER :: i,j,k
  INTEGER :: n                                !número de cuerpos
  REAL :: dt, t_end, t, dt_out, t_out         !paso temporal, ...
  REAL :: rs, r2, r3                          !posicion sun??, r al cuadrado y al cubo
  type(vector3d) :: rji                       !posición entre dos cuerpos

  type(particle), allocatable :: p(:)         !particulas p tienen (%p, v, a, m)

  print*, 'Por favor, introduzca el salto temporal: '
  READ*, dt
  print*, 'Por favor, introduza el dt_out: '
  READ*, dt_out
  print*, 'Por favor, introduzca el tiempo final: '
  READ*, t_end
  print*, 'Por favor, introduzca el numero de cuerpos: '
  READ*, n

  
  ALLOCATE(p(n))

  DO i = 1, n
     print*, 'Por favor, introduzca la posicion (tres numeros)',i,':'
     read*, p(i)%p
  END DO

  DO i = 1, n
     print*, 'Por favor, introduzca la velocidad (tres numeros)',i,':'
     read*,  p(i)%v
  END DO

  DO i = 1, n
     print*, 'Por favor, introduzca la masa del cuerpo (un numero)',i,':'
     read*, p(i)%m
  END DO

  do i = 1, n
     p(i)%a = vector3d(0.0, 0.0, 0.0)
  end do
  
  DO i = 1, n
     DO j = i+1, n
        rji = distance(p(j)%p, p(i)%p)
        r2 = modulus(rji)
        r3 = r2 * SQRT(r2)
        p(i)%a = p(i)%a .vsv. (p(j)%m .rpv. rji) .ver. r3
        p(j)%a = p(j)%a .vrv. (p(i)%m .rpv. rji) .ver. r3
     END DO
  END DO

  t = 0.0                                                  !new
  t_out = 0.0
  DO while (t <= t_end)                                    !changed bc fortran doesnt allow non-entire variables at do-loops
     do i = 1, n                                           !added
     p(i)%v = p(i)%v .vsv. (p(i)%a .vpr. (dt/2))
     p(i)%p = p(i)%p .psv. (p(i)%v .vpr. dt)
     p(i)%a = vector3d(0.0, 0.0, 0.0)
     end do
     DO i = 1, n
        DO j = i+1, n
        rji = distance(p(j)%p, p(i)%p)
        r2 = modulus(rji)
        r3 = r2 * SQRT(r2)
        p(i)%a = p(i)%a .vsv. (p(j)%m .rpv. rji) .ver. r3
        p(j)%a = p(j)%a .vrv. (p(i)%m .rpv. rji) .ver. r3
        END DO
     END DO
     do i = 1, n                                            !added
      p(i)%v = p(i)%v .vsv. (p(i)%a .vpr. (dt/2))
     end do
   
     t_out = t_out + dt
      IF (t_out >= dt_out) THEN
         DO i = 1,n
            PRINT*, p(i)%p
         END DO

         t_out = 0.0
      END IF
  END DO
  
END PROGRAM leapfrog
