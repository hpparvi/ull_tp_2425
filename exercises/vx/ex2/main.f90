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
  
  TYPE (CELL), POINTER :: head, temp_cell

  
  !! Lectura de datos
!!!!!!!!!!!!!!!!!!!
  
  !Input and output files.
  character(len=30) :: data                          
  character(len=30) :: orbits

  data = 'data_input.dat'                            
  orbits = 'data_output.dat'

  !Read the needed parameters.
  OPEN(10, file=data, status='old', action='read')

  READ(10, *) dt
  print *, "The read value of the time step (dt) is ", dt

  READ(10, *) dt_out
  print *, "The read value of the output time step (dt_out) is", dt_out

  READ(10, *) t_end
  print *, "The read value of the total simulation time (t_end) is ", t_end

  READ(10, *) n
  print *, "The read value of the number of bodies (n) is ", n
  

  ALLOCATE(a(n,3))

    !Read positions, velocities, and masses of particles.
  DO i = 1, n
     READ(10, *) p(i)%p
     print *, "The read value of the position of particle", i, "is:", p(i)%p

     READ(10, *) p(i)%v
     print *, "The read value of the velocity of particle", i, "is:", p(i)%v

     READ(10, *) p(i)%m
     print *, "The read value of the mass of particle", i, "is:", p(i)%m
  END DO

  CLOSE(10)

  
  
  !! InicializaciÂ´on head node
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(head)
  CALL Calculate_ranges(head)
  head%type = 0
  CALL Nullify_Pointers(head)


  !! Creación del árbol inicial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1, n
     CALL Find_Cell(head, temp_cell, r(i, :))
     CALL Place_Cell(temp_cell, r(i, :), i)
  END DO
  CALL Borrar_empty_leaves(head)
  CALL Calculate_masses(head)
  
  !! Calcular aceleraciones iniciales
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  a = 0.0
  CALL Calculate_forces(head)
  
  !! Bucle principal
!!!!!!!!!!!!!!!!!!
  t_out = 0.0
  t = 0.0
  DO WHILE (t <= t_end)
     v = v + a * dt / 2
     r = r + v * dt
     
     !! Las posiciones han cambiado, por lo que tenemos que borrar
     !! y reinicializar el árbol
     CALL Borrar_tree(head)
     CALL Calculate_ranges(head)
     head%type = 0
     CALL Nullify_Pointers(head)
     DO i = 1, n
        CALL Find_Cell(head, temp_cell, r(i, :))
        CALL Place_Cell(temp_cell, r(i, :), i)
     END DO
     CALL Borrar_empty_leaves(head)
     CALL Calculate_masses(head)
     
     a = 0.0
     CALL Calculate_forces(head)
     v = v + a * dt / 2
     
     t_out = t_out + dt
     IF (t_out >= dt_out) THEN
        DO i = 1, 10
           PRINT *, r(i, :)
        END DO
        PRINT *, "-----------------------------------"
        PRINT *, ""
        t_out = 0.0
     END IF
     
     t = t + dt
  END DO
END PROGRAM tree
