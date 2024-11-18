PROGRAM particulas
  USE EX2_deriv_types
  USE geometry
  IMPLICIT NONE
  
  INTEGER :: I,N
  INTEGER :: values(1:8), k
  INTEGER, dimension(:), allocatable :: seed
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: pt   !Particles
  TYPE(vector3d)  :: prueba    !Accelerations
  REAL :: mass, dt, dt_out, t_end
  CHARACTER(len=25) :: output
  
  output =  "initial_conditions.dat"
  
  CALL date_and_time(values=values)
  CALL random_seed(size=k)
  
  ALLOCATE(seed(1:k))
  SEED(:) = values(8)
  CALL random_seed(put=seed)

  
  OPEN(1, file = output, status = 'replace', action = 'write')
  PRINT*, "Valor de dt?"
  READ*, dt
  WRITE(1,'(F12.6)') dt
  PRINT*, "Valor de dt_out?"
  READ*, dt_out
  WRITE(1,'(F12.6)') dt_out
  PRINT*, "Valor de t_end?"
  READ*, t_end
  WRITE(1,'(F12.6)') t_end
  PRINT*, "Number of bodies?"
  READ*, N
  WRITE(1,'(I5)') N
  mass = 1.0 / N
  ALLOCATE(pt(N))
  DO I= 1,N
    CALL random_number(pt(I)%p%x)
    DO
      CALL random_number(pt(I)%p%y)
      IF ((pt(I)%p%x**2 + pt(I)%p%y**2) .LE. 1) EXIT
    END DO
    DO
      CALL random_number(pt(I)%p%z)
      IF ((pt(I)%p%x**2 + pt(I)%p%y**2 + pt(I)%p%z**2) .LE. 1) EXIT
    END DO
  pt(I)%m = mass
  pt(I)%v = vector3d(0,0,0)
  WRITE(1,'(7F12.6)') pt(i)%m, pt(i)%p, pt(i)%v
  END DO
  
  
END PROGRAM particulas