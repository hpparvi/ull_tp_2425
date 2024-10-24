PROGRAM leapfrog
    USE geometry
    USE particle
    IMPLICIT NONE
    INTEGER :: i,j,k
    INTEGER :: n
    REAL :: dt, t_end, t, dt_out, t_out
    REAL :: rs, r2, r3
    REAL, DIMENSION(:), ALLOCATABLE :: m
    REAL, DIMENSION(:,:), ALLOCATABLE :: r,v,a
    REAL, DIMENSION(3) :: rji
    READ*, dt
    READ*, dt_out
    READ*, t_end
    READ*, n
    ALLOCATE(m(n))
    ALLOCATE(r(n,3))
    ALLOCATE(v(n,3))
    ALLOCATE(a(n,3))
    DO i = 1, n
        READ*, m(i), r(i,:),v(i,:)
    END DO
    a = 0.0
    DO i = 1,n
        DO j = i+1,n
            rji = r(j,:) - r(i,:)
            r2 = SUM(rji**2)
            r3 = r2 * SQRT(r2)
            a(i,:) = a(i,:) + m(j) * rji / r3
            a(j,:) = a(j,:) - m(i) * rji / r3
        END DO
    END DO
    t_out = 0.0
    DO t = 0.0, t_end, dt
        v = v + a * dt/2
        r = r + v * dt
        a = 0.0
        DO i = 1,n
            DO j = i+1,n
                rji = r(j,:) - r(i,:)
                r2 = SUM(rji**2)
                r3 = r2 * SQRT(r2)
                a(i,:) = a(i,:) + m(j) * rji / r3
                a(j,:) = a(j,:) - m(i) * rji / r3
            END DO
        END DO
        v = v + a * dt/2
        t_out = t_out + dt
        IF (t_out >= dt_out) THEN
            DO i = 1,n
                PRINT*, r(i,:)
            END DO
            t_out = 0.0
        END IF
    END DO
END PROGRAM leapfrog
