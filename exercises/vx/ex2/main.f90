PROGRAM tree
  use iso_fortran_env
  use geometry
  use definitions
  use barneshut
  IMPLICIT NONE
!  INTEGER :: i,j,k,n
!  REAL :: dt, t_end, t, dt_out, t_out, rs, r2, r3
!  REAL(real64), PARAMETER :: theta = 1
!  REAL(real64), DIMENSION(:), ALLOCATABLE :: m
!  REAL(real64), DIMENSION(:,:), ALLOCATABLE :: r,v,a
!  REAL(real64), DIMENSION(3) :: rji
!  type(particle3d), allocatable :: p(:)               !Particles of the problem.
!  type(vector3d), allocatable :: a(:)                 !Their accelerations.
  type(vector3d), allocatable :: a(:)
  type(particle3d), allocatable :: p(:) 
  TYPE (CELL), POINTER :: head, temp_cell
  type(vector3d) :: rji

  
  !! Lectura de datos
!!!!!!!!!!!!!!!!!!!
  
  !Input and output files.
  character(len=30) :: data                          
  character(len=30) :: orbits

  data = 'data_input_b.dat'                            
  orbits = 'data_output.dat'

  !Read the needed parameters.
  !Read the needed parameters.
  OPEN(10, file=data, status='old', action='read')

  READ(10, *) dt
  PRINT *, "The read value of the time step (dt) is ", dt

  READ(10, *) dt_out
  PRINT *, "The read value of the output time step (dt_out) is", dt_out

  READ(10, *) t_end
  PRINT *, "The read value of the total simulation time (t_end) is ", t_end

  READ(10, *) n
  PRINT *, "The read value of the number of bodies (n) is ", n

  ALLOCATE(a(n))
  ALLOCATE(p(n))

  !Read positions, velocities, and masses of particles.
  DO i = 1, n
     READ(10, *) p(i)%m, p(i)%p%x, p(i)%p%y, p(i)%p%z, &
                                      p(i)%v%x, p(i)%v%y, p(i)%v%z
     PRINT *, "Particle", i, ":"
     PRINT *, "  Position:", p(i)%p
     PRINT *, "  Velocity:", p(i)%v
     PRINT *, "  Mass:    ", p(i)%m
  END DO

  CLOSE(10)

  
  
  !! InicializaciÂ´on head node
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(head)
  CALL Calculate_ranges(head, p)
  head%type = 0
  CALL Nullify_Pointers(head)


  !! Creación del árbol inicial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1, n
     CALL Find_Cell(head, temp_cell, p(i))
     CALL Place_Cell(temp_cell, p(i), i)
  END DO
  CALL Borrar_empty_leaves(head)
  CALL Calculate_masses(head, p)
  
  !! Calcular aceleraciones iniciales
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1, n
     a(i) = vector3d(0.0, 0.0, 0.0)
  END DO

  open(11, file = orbits, status = 'old', action = 'write')
  CALL Calculate_forces(head, n, p, rji, a)
  
  !! Bucle principal
!!!!!!!!!!!!!!!!!!
  t_out = 0.0
  t = 0.0
  DO WHILE (t <= t_end)
     do i = 1,n
        p(i)%v = p(i)%v + a(i) * (dt / 2)
     end do
     do i = 1, n
        p(i)%p = p(i)%p + p(i)%v * dt
     end do

     
     !! Las posiciones han cambiado, por lo que tenemos que borrar
     !! y reinicializar el árbol
     
     CALL Borrar_tree(head)
     CALL Calculate_ranges(head, p)
     head%type = 0
     CALL Nullify_Pointers(head)
     DO i = 1, n
        CALL Find_Cell(head, temp_cell, p(i))
        CALL Place_Cell(temp_cell, p(i), i)
     END DO
     CALL Borrar_empty_leaves(head)
     CALL Calculate_masses(head, p)
     
     DO i = 1, n
        a(i) = vector3d(0.0, 0.0, 0.0)
     END DO
     CALL Calculate_forces(head, n, p, rji, a)
     do i=1,n
        p(i)%v = p(i)%v + a(i) * (dt / 2)
     end do
     t_out = t_out + dt
     IF (t_out >= dt_out) THEN
        WRITE(11, '(E12.2)', ADVANCE='NO') t
        DO i = 1, n
           WRITE(11, '(3X, E12.6, 3X, E12.6, 3X, E12.6)', ADVANCE='NO') &
               p(i)%p%x, p(i)%p%y, p(i)%p%z
        END DO
        WRITE(11, *)
        t_out = 0.0
     END IF
     
     t = t + dt
  END DO
END PROGRAM tree
