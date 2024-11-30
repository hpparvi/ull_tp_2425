PROGRAM particulas
  USE EX2_deriv_types
  USE geometry
  IMPLICIT NONE
  
  INTEGER :: I,N
  INTEGER :: values(1:8), k
  INTEGER, dimension(:), allocatable :: seed
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: pt   !Particles
  TYPE(vector3d)  :: prueba    !Accelerations
  TYPE(point3d)  :: c_o_m    !Accelerations
  REAL(REAL64) :: mass, dt, dt_out, t_end, vel
  REAL(REAL64) :: G = 0.45 ! Factor debido a que no toda la masa está en el centro de masa
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
  c_o_m = point3d(0,0,0)
  ALLOCATE(pt(N))
  DO I= 1,N
    DO
      CALL random_number(pt(I)%p%x)
      pt(I)%p%x = 2.0 * (pt(I)%p%x - 0.5)
      CALL random_number(pt(I)%p%y)
      pt(I)%p%y = 2.0 * (pt(I)%p%y - 0.5)
      CALL random_number(pt(I)%p%z)
      pt(I)%p%z = 2.0 * (pt(I)%p%z - 0.5)
      IF ((pt(I)%p%x**2 + pt(I)%p%y**2 + pt(I)%p%z**2) .LE. 1) EXIT
    END DO
    c_o_m%x = c_o_m%x + pt(I)%p%x/N
    c_o_m%y = c_o_m%y + pt(I)%p%y/N
    c_o_m%z = c_o_m%z + pt(I)%p%z/N
    CALL random_number(pt(I)%v%x)
    CALL random_number(pt(I)%v%y)
    CALL random_number(pt(I)%v%z)
    pt(I)%m = mass
    !pt(I)%v%x = 2.0 * (pt(I)%v%x - 0.5)/(N*N)
    !pt(I)%v%y = 2.0 * (pt(I)%v%y - 0.5)/(N*N)
    !pt(I)%v%z = 2.0 * (pt(I)%v%z - 0.5)/(N*N)
  END DO
  DO I= 1,N
    vel = sqrt(1/distance(pt(I)%p,c_o_m))*G
    pt(I)%v = vel*orthv(c_o_m - pt(I)%p,vector3d(0,0,1))
    WRITE(1,'(7E14.6)') pt(i)%m, pt(i)%p, pt(i)%v
  END DO
  
END PROGRAM particulas