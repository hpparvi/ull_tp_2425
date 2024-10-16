PROGRAM leapfrog

  use geometry
  use particles

  IMPLICIT NONE

  INTEGER :: i,j,k, i_t
  INTEGER :: n, n_t
  double precision :: dt, t_end, t, dt_out, t_out
  double precision :: rs, r2, r3
  
  type(particle), dimension(:), allocatable :: partics
!  type(vector3d), dimension(:), allocatable :: a_v
  type(vector3d)                            :: rji_v
  double precision                          :: r3_2
  
  double precision, DIMENSION(:), ALLOCATABLE :: m
  double precision, DIMENSION(:,:), ALLOCATABLE :: r,v,a
  double precision, DIMENSION(3) :: rji
  
  READ*, dt
  READ*, dt_out
  READ*, t_end
  READ*, n

  n_t = t_end/dt
  
  allocate(partics(n))
!  allocate(a_v(n))
  
  ALLOCATE(m(n))
  ALLOCATE(r(n,3))
  ALLOCATE(v(n,3))
  ALLOCATE(a(n,3))
  
  DO i = 1, n
     READ*, m(i), r(i,:),v(i,:)
     print*, "Write again for checks"
     READ*, partics(i)%m, partics(i)%p, partics(i)%v
  END DO

  do i = 1,n
     partics(i)%a = vector3d(0,0,0)
  end do
  
  DO i = 1,n
     DO j = i+1,n
        rji_v = vecpp(partics(j)%p, partics(i)%p)
        r3_2 = norm(rji_v)**3
        partics(i)%a = partics(i)%a + (partics(j)%m * rji_v / r3_2)
        partics(j)%a = partics(j)%a - (partics(i)%m * rji_v / r3_2)
        
        rji = r(j,:) - r(i,:)
        r2 = SUM(rji**2)
        r3 = r2 * SQRT(r2)
        a(i,:) = a(i,:) + m(j) * rji / r3
        a(j,:) = a(j,:) - m(i) * rji / r3
     END DO
  END DO
  
  t_out = 0.0
  
  DO i_t = 0, n_t
     t = i_t * dt
     partics(:)%v = partics(:)%v + (partics(:)%a * (dt/2))
     partics(:)%p = partics(:)%p + (partics(:)%v * dt)
     
     v = v + a * dt/2
     r = r + v * dt

     a = 0.0
     DO i = 1,n
        DO j = i+1,n
           rji_v = vecpp(partics(j)%p, partics(i)%p)
           r3_2 = norm(rji_v)**3
           partics(i)%a = partics(i)%a + (partics(j)%m * rji_v / r3_2)
           partics(j)%a = partics(j)%a - (partics(i)%m * rji_v / r3_2)
        
           rji = r(j,:) - r(i,:)
           r2 = SUM(rji**2)
           r3 = r2 * SQRT(r2)
           a(i,:) = a(i,:) + m(j) * rji / r3
           a(j,:) = a(j,:) - m(i) * rji / r3
        END DO
     END DO
     
     partics(:)%v = partics(:)%v + (partics(:)%a * (dt/2))
     
     v = v + a * dt/2
     
     t_out = t_out + dt
     IF (t_out >= dt_out) THEN
        DO i = 1,n
           print*, "Time = ", t
           PRINT*, r(i,:)
           print*, partics(i)%p
        END DO
        t_out = 0.0
     END IF
     
  END DO
  
END PROGRAM leapfrog
