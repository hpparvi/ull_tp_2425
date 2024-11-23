PROGRAM ex2

  USE geometry
  USE ex2_deriv_types
  USE barnes_hut
  USE ISO_FORTRAN_ENV
  
  IMPLICIT NONE      
  
  REAL(REAL64) :: dt, t_end, t, dt_out, t_out    
  character(len=25) :: input, output, orbit
  
!! Lectura de datos
!!!!!!!!!!!!!!!!!!!
  input = "initial_conditions.dat"
  output = "output.dat"
  
  OPEN (1,file = input, status = 'old', action = 'read')
  

  READ(1,*) dt
  READ(1,*) dt_out
  READ(1,*) t_end
  READ(1,*) n
  ALLOCATE(pt(n))
  ALLOCATE(a(n))
  DO i = 1, n
    READ(1,*) pt(i)%m, pt(i)%p, pt(i)%v
  END DO
  
!! Inicialización head node
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(head)
  CALL Calculate_ranges(head)
  head%type = 0
  CALL Nullify_Pointers(head)
  
  
!! Creaci´on del ´arbol inicial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO i = 1,n
    CALL Find_Cell(head,temp_cell,pt(i)%p)
    CALL Place_Cell(temp_cell,pt(i)%p,i)
  END DO
  CALL Borrar_empty_leaves(head)
  CALL Calculate_masses(head)


!! Calcular aceleraciones iniciales
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  a = vector3d(0,0,0)
  CALL Calculate_forces(head)

!! Bucle principal
!!!!!!!!!!!!!!!!!!
  
  OPEN(3, file = output, status = 'replace', action = 'write')
  t = 0.0
  t_out = 0.0
  DO 
    t = t + dt
    t_out = t_out + dt
    pt%v = pt%v + a * (dt/2)
    pt%p = pt%p + pt%v * dt
    !! Las posiciones han cambiado, por lo que tenemos que borrar
    !! y reinicializar el ´arbol

    CALL Borrar_tree(head)
    CALL Calculate_ranges(head)
    head%type = 0
    CALL Nullify_Pointers(head)
    DO i = 1,n
      CALL Find_Cell(head,temp_cell,pt(i)%p)
      CALL Place_Cell(temp_cell,pt(i)%p,i)
    END DO
    CALL Borrar_empty_leaves(head)
    CALL Calculate_masses(head)
    a = vector3d(0,0,0)
    CALL Calculate_forces(head)
    pt%v = pt%v + a * (dt/2)
    IF (t_out >= dt_out) THEN
      write (3, '(3F12.4)', advance='no') t
      DO i = 1,n
        IF(i<n) THEN
          write (3, '(3F12.6)', advance='no') pt(i)%p
        ELSE
          write (3, '(3F12.6)') pt(i)%p
        END IF
      END DO
      t_out = 0.0 ! reset of output parameter 
    END IF
    ! Exit conditions
    IF (t>= t_end) THEN
      EXIT
    END IF
  END DO
  
 
END PROGRAM ex2