program particulas     !This program initializes a random mass, positions and velocity for the number of particles desired by the user. The mass and position are
                       !multiplied by 100 in order to avoid having a lot of particles in a really small space.

  use iso_fortran_env  !This module ensures all variables are defined as 64-bit.
  IMPLICIT NONE
  
  INTEGER(int64) :: I,N
  INTEGER :: values(1:8), k
  INTEGER, dimension(:), allocatable :: seed
  REAL(real64) :: mass,rx,ry,rz
  real(real64) :: vx, vy, vz
  character(len=30) :: data
  real(real64) :: a, b, c, d, e, f

  data = 'data_input.dat'

  open(10, file=data, status='replace', action='write')
  
  CALL date_and_time(values=values)
  CALL random_seed(size=k)
  ALLOCATE(seed(1:k))
  SEED(:) = values(8)
  CALL random_seed(put=seed)

  
  PRINT*, "Number of bodies?"
  READ*, N
  WRITE(10, *) 0.001    ! dt
  WRITE(10, *) 0.1      ! dt_out
  WRITE(10, *) 100.0    ! t_end
  WRITE(10, *) n        ! Número de cuerpos
  
  DO I= 1,N
     ! Random mass, from 0 to 1 and multiplied by 100
     CALL random_number(mass)
     mass = mass * 100

     ! Random positions, from 0 to 1 and multiplied by a factor 100 and by a random +/- 1
     CALL random_number(rx)
     call random_number(a)
     if (a <= 0.5) then
        a = -a
     end if
     rx = rx*100 * a

     CALL random_number(ry)
     call random_number(b)
     if (b <= 0.5) then
        b = -b
     end if
     ry = ry*100 * b
        
     CALL random_number(rz)
     call random_number(c)
     if (c <= 0.5) then
        c = -c
     end if
     rz = rz*100 * c

     ! Random velocities, from 0 to 1 and multiplied by a random +/- 1
     CALL random_number(vz)
     call random_number(d)
     if (d <= 0.5) then
        d = -d
     end if
     vz = vz * d
     
     CALL random_number(vy)
     call random_number(e)
     if (e <= 0.5) then
        e = -e
     end if
     vy = vy * e

     CALL random_number(vx)
     call random_number(f)
     if (f <= 0.5) then
        f = -f
     end if
     vx = vx * f

     ! Write the data input
     WRITE(10,'(F5.2, 6F11.3)') mass, rx,ry,rz, vx,vy,vz

  END DO

  close(10)
  
END program particulas
