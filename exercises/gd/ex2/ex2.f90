PROGRAM tree

  USE geometry
  USE particle
  USE ISO_FORTRAN_ENV
  
  IMPLICIT NONE
  INTEGER :: i,j,k,n
  REAL(REAL64) :: dt, t_end, t, dt_out, t_out, rs, r2, r3
  REAL(REAL64), PARAMETER :: theta = 1
  !REAL, DIMENSION(:), ALLOCATABLE :: m
  !REAL, DIMENSION(:,:), ALLOCATABLE :: r,v,a                            
  character(len=25) :: input, output, orbit
  
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: pt   !Particles
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a    !Accelerations
  TYPE(vector3d) :: r
  !REAL, DIMENSION(3) :: rji
  
  TYPE RANGE
    TYPE(point3d)  :: min,max
    !REAL, DIMENSION(3) :: min,max
  END TYPE RANGE
  
  TYPE CPtr
    TYPE(CELL), POINTER :: ptr
  END TYPE CPtr
  
  TYPE CELL
    TYPE (RANGE) :: range
    !REAL, DIMENSION(3) :: part
    INTEGER :: pos
    INTEGER :: type !! 0 = no particle; 1 = particle; 2 = conglomerado
    !REAL :: mass
    TYPE(particle3d) :: pt 
    
    !REAL, DIMENSION(3) :: c_o_m
    TYPE(vector3d)  :: c_o_m
    TYPE (CPtr), DIMENSION(2,2,2) :: subcell
  END TYPE CELL
  
  TYPE (CELL), POINTER :: head, temp_cell
  
!! Lectura de datos
!!!!!!!!!!!!!!!!!!!
  input = "initial_conditions.dat"
  output = "output.dat"
  orbit = "orbit.dat"
  
  OPEN (1,file = input, status = 'old', action = 'read')
  

  READ(1,*) dt
  READ(1,*) dt_out
  READ(1,*) t_end
  READ(1,*) n
  !ALLOCATE(m(n))
  !ALLOCATE(r(n,3))
  !ALLOCATE(v(n,3))
  ALLOCATE(pt(n))
  ALLOCATE(a(n))
  DO i = 1, n
    READ(1,*) pt(i)%m, pt(i)%p, pt(i)%v
    !READ*, m(i), r(i,:),v(i,:)
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
  
  OPEN(2, file = output, status = 'replace', action = 'write')
  OPEN(3, file = orbit, status = 'replace', action = 'write')
  t = 0.0
  t_out = 0.0
  DO 
  !PRINT*, 'Entra do', t
    t = t + dt
    t_out = t_out + dt
    !v = v + a * dt/2
    !r = r + v * dt
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
    !v = v + a * dt/2
    !PRINT*, 'aceleraciones', a
    pt%v = pt%v + a * (dt/2)
    IF (t_out >= dt_out) THEN
      write (3, '(3F12.4)', advance='no') t
      DO i = 1,n
          write (2, '(3F12.4)', advance='no') t
          write (2, '(I5)', advance='no') i
          write (2, '(3F12.6)', advance='no') pt(i)%p
          write (2, '(3F12.6)', advance='no') pt(i)%v
          write (2, '(3F12.6)') pt(i)%m
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
  
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_Ranges !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calcula los rangos de las part´ıculas en la
!! matriz r en las 3 dimensiones y lo pone en la
!! variable apuntada por goal
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Calculate_Ranges(goal)
    TYPE(CELL),POINTER :: goal
    !REAL, DIMENSION(3) :: mins,maxs,medios
    TYPE(point3d) :: mins,maxs,medios
    REAL(REAL64) :: span
    mins%x = MINVAL(pt%p%x)
    mins%y = MINVAL(pt%p%y)
    mins%z = MINVAL(pt%p%z)
    maxs%x = MAXVAL(pt%p%x)
    maxs%y = MAXVAL(pt%p%y)
    maxs%z = MAXVAL(pt%p%z)
  
    !mins = MINVAL(pt%p,DIM=1)
    !maxs = MAXVAL(pt%p,DIM=1)
      ! Al calcular span le sumo un 10% para que las
      ! particulas no caigan justo en el borde
    span = MAX(maxs%x - mins%x, maxs%y - mins%y, maxs%z - mins%z) * 1.1
    medios%x = (maxs%x + mins%x) / 2.0
    medios%y = (maxs%y + mins%y) / 2.0
    medios%z = (maxs%z + mins%z) / 2.0
    goal%range%min%x = medios%x - span/2.0
    goal%range%min%y = medios%y - span/2.0
    goal%range%min%z = medios%z - span/2.0
    goal%range%max%x = medios%x + span/2.0
    goal%range%max%y = medios%y + span/2.0
    goal%range%max%z = medios%z + span/2.0
  END SUBROUTINE Calculate_Ranges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Find_Cell !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Encuentra la celda donde colocaremos la particula.
!! Si la celda que estamos considerando no tiene
!! particula o tiene una particula, es esta celda donde
!! colocaremos la particula.
!! Si la celda que estamos considerando es un "conglomerado",
!! buscamos con la funci´on BELONGS a que subcelda de las 8
!! posibles pertenece y con esta subcelda llamamos de nuevo
!! a Find_Cell
!!
!! NOTA: Cuando se crea una celda "conglomerado" se crean las
!! 8 subceldas, por lo que podemos asumir que siempre existen
!! las 8. Las celdas vac´ıas se borran al final del todo, cuando
!! todo el ´arbol ha sido ya creado.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  RECURSIVE SUBROUTINE Find_Cell(root,goal,part)
    !REAL, DIMENSION(3) :: part
    TYPE(point3d) :: part
    TYPE(CELL),POINTER :: root,goal,temp
    INTEGER :: i,j,k
      SELECT CASE (root%type)
        CASE (2)
          out: DO i = 1,2
            DO j = 1,2
              DO k = 1,2
                IF (Belongs(part,root%subcell(i,j,k)%ptr)) THEN
                  CALL Find_Cell(root%subcell(i,j,k)%ptr,temp,part)
                  goal => temp
                  EXIT out
              END IF
            END DO
          END DO
        END DO out
      CASE DEFAULT
        goal => root
    END SELECT
  END SUBROUTINE Find_Cell
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Place_Cell !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Se ejecuta tras Find_Cell, en la celda que
!! esa funci´on nos devuelve, por lo que siempre
!! es una celda de tipo 0 (sin particula) o de tipo 1
!! (con una particula). En el caso de que es una celda
!! de tipo 1 habra que subdividir la celda y poner en
!! su lugar las dos particulas (la que originalmente
!! estaba, y la nueva).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE Place_Cell(goal,part,n)
    TYPE(CELL),POINTER :: goal,temp
    !REAL, DIMENSION(3) :: part
    TYPE(point3d) :: part
    INTEGER :: n
    
    SELECT CASE (goal%type)
      CASE (0)
        goal%type = 1
        goal%pt%p = part
        goal%pos = n
      CASE (1)
        CALL Crear_Subcells(goal)
        CALL Find_Cell(goal,temp,part)
        CALL Place_Cell(temp,part,n)
      CASE DEFAULT
        print*,"SHOULD NOT BE HERE. ERROR!" 
      END SELECT
  END SUBROUTINE Place_Cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Crear_Subcells !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Esta funcion se llama desde Place_Cell y
!! solo se llama cuando ya hay una particula
!! en la celda, con lo que la tenemos que
!! subdividir. Lo que hace es crear 8 subceldas
!! que "cuelgan" de goal y la particula que
!! estaba en goal la pone en la subcelda que
!! corresponda de la 8 nuevas creadas.
!!
!! Para crear las subceldas utilizar las funciones
!! CALCULAR_RANGE, BELONGS y NULLIFY_POINTERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Crear_Subcells(goal)
    TYPE(CELL), POINTER :: goal
    !REAL,DIMENSION(3) :: part
    TYPE(point3d) :: part
    INTEGER :: i,j,k,n
    INTEGER, DIMENSION(3) :: octant
    part = goal%pt%p
    goal%type=2
    DO i = 1,2
      DO j = 1,2
        DO k = 1,2
          octant = (/i,j,k/)
          ALLOCATE(goal%subcell(i,j,k)%ptr)
          goal%subcell(i,j,k)%ptr%range%min = Calcular_Range (0,goal,octant)
          goal%subcell(i,j,k)%ptr%range%max = Calcular_Range (1,goal,octant)
          IF (Belongs(part,goal%subcell(i,j,k)%ptr)) THEN
            goal%subcell(i,j,k)%ptr%pt%p = part
            goal%subcell(i,j,k)%ptr%type = 1
            goal%subcell(i,j,k)%ptr%pos = goal%pos
          ELSE
            goal%subcell(i,j,k)%ptr%type = 0
          END IF
        CALL Nullify_Pointers(goal%subcell(i,j,k)%ptr)
        END DO
      END DO
    END DO
  END SUBROUTINE Crear_Subcells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Nullify_Pointers !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Simplemente me NULLIFYca los punteros de
!! las 8 subceldas de la celda "goal"
!!
!! Se utiliza en el bucle principal y por
!! CREAR_SUBCELLS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Nullify_Pointers(goal)
    TYPE(CELL), POINTER :: goal
    INTEGER :: i,j,k
    DO i = 1,2
      DO j = 1,2
        DO k = 1,2
          NULLIFY(goal%subcell(i,j,k)%ptr)
        END DO
      END DO
    END DO
  END SUBROUTINE Nullify_Pointers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Belongs !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Devuelve TRUE si la particula "part" est´a
!! dentro del rango de la celda "goal"
!!
!! Utilizada por FIND_CELL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION Belongs (part,goal)
    !REAL, DIMENSION(3) :: part
    TYPE(point3d) :: part
    TYPE(CELL), POINTER :: goal
    LOGICAL :: Belongs
    IF (part%x > goal%range%min%x .AND. &
        part%x <= goal%range%max%x .AND. &
        part%y > goal%range%min%y .AND. &
        part%y <= goal%range%max%y .AND. &
        part%z > goal%range%min%z .AND. &
        part%z <= goal%range%max%z) THEN
      Belongs = .TRUE.
    ELSE
      Belongs = .FALSE.
    END IF
  END FUNCTION Belongs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calcular_Range !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Dado un octante "otctant" (1,1,1, 1,1,2 ... 2,2,2),
!! calcula sus rangos en base a los rangos de
!! "goal". Si "what" = 0 calcula los minimos. Si what=1
!! calcula los maximos.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION Calcular_Range (what,goal,octant)
    INTEGER :: what,n
    TYPE(CELL), POINTER :: goal
    INTEGER, DIMENSION(3) :: octant
    !REAL, DIMENSION(3) :: Calcular_Range, valor_medio
    TYPE(point3d) :: Calcular_Range, valor_medio
    
    valor_medio%x = (goal%range%min%x + goal%range%max%x) / 2.0
    valor_medio%y = (goal%range%min%y + goal%range%max%y) / 2.0
    valor_medio%z = (goal%range%min%z + goal%range%max%z) / 2.0
    
    !SELECT CASE (what)
    !CASE (0)
    !  WHERE (octant == 1)
    !    Calcular_Range = goal%range%min
    !  ELSEWHERE
    !    Calcular_Range = valor_medio
    !  ENDWHERE
    !CASE (1)
    !  WHERE (octant == 1)
    !    Calcular_Range = valor_medio
    !  ELSEWHERE
    !    Calcular_Range = goal%range%max
    !  ENDWHERE
    !END SELECT
    SELECT CASE (what)
    CASE (0)
        IF (octant(1) == 1) THEN
            Calcular_Range%x = goal%range%min%x
        ELSE
            Calcular_Range%x = valor_medio%x
        END IF

        IF (octant(2) == 1) THEN
            Calcular_Range%y = goal%range%min%y
        ELSE
            Calcular_Range%y = valor_medio%y
        END IF

        IF (octant(3) == 1) THEN
            Calcular_Range%z = goal%range%min%z
        ELSE
            Calcular_Range%z = valor_medio%z
        END IF

    CASE (1)
        ! Asignación condicional por componente para el caso what = 1
        IF (octant(1) == 1) THEN
            Calcular_Range%x = valor_medio%x
        ELSE
            Calcular_Range%x = goal%range%max%x
        END IF

        IF (octant(2) == 1) THEN
            Calcular_Range%y = valor_medio%y
        ELSE
            Calcular_Range%y = goal%range%max%y
        END IF

        IF (octant(3) == 1) THEN
            Calcular_Range%z = valor_medio%z
        ELSE
            Calcular_Range%z = goal%range%max%z
        END IF
    END SELECT
  END FUNCTION Calcular_Range


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Borrar_empty_leaves !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Se llama una vez completado el arbol para
!! borrar (DEALLOCATE) las celdas vac´ıas (i.e.
!! sin part´ıcula).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE Borrar_empty_leaves(goal)
    TYPE(CELL),POINTER :: goal
    INTEGER :: i,j,k
    IF (ASSOCIATED(goal%subcell(1,1,1)%ptr)) THEN
      DO i = 1,2
        DO j = 1,2
          DO k = 1,2
            CALL Borrar_empty_leaves(goal%subcell(i,j,k)%ptr)
            IF (goal%subcell(i,j,k)%ptr%type == 0) THEN
              DEALLOCATE (goal%subcell(i,j,k)%ptr)
            END IF
          END DO
        END DO
      END DO
    END IF
  END SUBROUTINE Borrar_empty_leaves

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Borrar_tree !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Borra el arbol completo, excepto la "head".
!!
!! El arbol se ha de regenerar continuamente,
!! por lo que tenemos que borrar el antiguo
!! para evitar "memory leaks".
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE Borrar_tree(goal)
  TYPE(CELL),POINTER :: goal
  INTEGER :: i,j,k
    DO i = 1,2
      DO j = 1,2
        DO k = 1,2
          IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
            CALL Borrar_tree(goal%subcell(i,j,k)%ptr)
            DEALLOCATE (goal%subcell(i,j,k)%ptr)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE Borrar_tree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_masses !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Nos calcula para todas las celdas que cuelgan
!! de "goal" su masa y su center-of-mass.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE Calculate_masses(goal)
    TYPE(CELL),POINTER :: goal
    INTEGER :: i,j,k
    REAL(REAL64) :: mass
    !REAL, DIMENSION(3) :: c_o_m
    ! Tratamos el centro de masa como el vector que une el punto 0,0,0 con la posición del centro de masa.
    TYPE(vector3d) :: c_o_m
    TYPE(point3d) :: V3d_0 = point3d(0,0,0)

    !goal%mass = 0
    !goal%c_o_m = 0
    goal%pt%m = 0
    goal%c_o_m = vector3d(0,0,0)
    SELECT CASE (goal%type)
    CASE (1)
      !goal%mass = m(goal%pos)
      !goal%c_o_m = r(goal%pos,:)
      goal%pt%m = pt(goal%pos)%m
      goal%c_o_m = V3d_0 - pt(goal%pos)%p
    CASE (2)
      DO i = 1,2
        DO j = 1,2
          DO k = 1,2
            IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
              CALL Calculate_masses(goal%subcell(i,j,k)%ptr)
              !mass = goal%mass
              mass = goal%pt%m
              !goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
              goal%pt%m = goal%pt%m + goal%subcell(i,j,k)%ptr%pt%m
              !goal%c_o_m = (mass * goal%c_o_m + &
              !goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m) / goal%mass
              goal%c_o_m = ((goal%c_o_m * mass) + &
              goal%subcell(i,j,k)%ptr%pt%m * goal%subcell(i,j,k)%ptr%c_o_m) / goal%pt%m
            END IF
          END DO
        END DO
      END DO
    END SELECT
  END SUBROUTINE Calculate_masses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_forces !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calcula las fuerzas de todas las particulas contra "head".
!! Se sirve de la funcion Calculate_forces_aux que es la
!! que en realidad hace los calculos para cada particula
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Calculate_forces(head)
    TYPE(CELL),POINTER :: head
    INTEGER :: i,j,k,start,end
    
    DO i = 1,n
      CALL Calculate_forces_aux(i,head)
    END DO
    
  END SUBROUTINE Calculate_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_forces_aux !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Dada una particula "goal" calcula las fuerzas
!! sobre ella de la celda "tree". Si "tree" es una
!! celda que contiene una sola particula el caso
!! es sencillo, pues se tratan de dos particulas.
!!
!! Si "tree" es una celda conglomerado, hay que ver primero
!! si l/D < theta. Es decir si el lado de la celda (l)
!! dividido entre la distancia de la particula goal
!! al center_of_mass de la celda tree (D) es menor que theta.
!! En caso de que asi sea, tratamos a la celda como una
!! sola particula. En caso de que no se menor que theta,
!! entonces tenemos que considerar todas las subceldas
!! de tree y para cada una de ellas llamar recursivamente
!! a Calculate_forces_aux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE Calculate_forces_aux(goal,tree)
    TYPE(CELL),POINTER :: tree
    INTEGER :: i,j,k,goal
    REAL(REAL64) :: l,D
    !para las fuerzas, convertinos el centro de masa en un punto
    TYPE(point3d) :: c_o_m
    TYPE(point3d) :: V3d_0 = point3d(0,0,0)
      !PRINT*, 'goal;' ,goal ,'pos;' ,tree%pos, 'type;' ,tree%type
    SELECT CASE (tree%type)
    CASE (1)
      IF (goal .NE. tree%pos) THEN
        !rji = tree%c_o_m - r(goal,:) - 
        !r2 = SUM(rji**2)
        !r3 = r2 * SQRT(r2)
        !a(goal,:) = a(goal,:) + m(tree%pos) * rji / r3
        c_o_m = tree%c_o_m - V3d_0
        !PRINT*, c_o_m , 'Centro de masas'
        r = pt(goal)%p - c_o_m 
        !PRINT*, r , 'r'
        a(goal) = a(goal) + (pt(tree%pos)%m *r)/ distance(pt(goal)%p,c_o_m)**3 !accel of particle i due to j
        !PRINT*, a(goal), 'ACELERACION:', goal
      END IF
    CASE (2)
          !! El rango tiene el mismo span en las 3 dimensiones
          !! por lo que podemos considerar una dimension cualquiera
          !! para calcular el lado de la celda (en este caso la
          !! dimension 1)
      c_o_m = tree%c_o_m - V3d_0
      l = tree%range%max%x - tree%range%min%x
      !rji = tree%c_o_m - r(goal,:)
      !r2 = SUM(rji**2)
      !D = SQRT(r2)
      D = distance(pt(goal)%p,c_o_m)
      IF (l/D < theta) THEN
        !! Si conglomerado, tenemos que ver si se cumple l/D < @
        !r3 = r2 * D
        !a(goal,:) = a(goal,:) + tree%mass * rji / r3
        r = pt(goal)%p - c_o_m 
        a(goal) = a(goal) + (pt(tree%pos)%m *r)/ distance(pt(goal)%p,c_o_m)**3 !accel of particle i due to j
        !PRINT*, a(goal), "CASE 2"
      ELSE
        DO i = 1,2
          DO j = 1,2
            DO k = 1,2
              IF (ASSOCIATED(tree%subcell(i,j,k)%ptr)) THEN
                !PRINT*, i,j,k ,'LISTEN'
                CALL Calculate_forces_aux(goal,tree%subcell(i,j,k)%ptr)
              END IF
            END DO
          END DO
        END DO
      END IF
    END SELECT
  END SUBROUTINE Calculate_forces_aux

! Subroutine to calculate acelerations 
 SUBROUTINE calc_acc
 ! Init acceleration
  a = vector3d(0,0,0)
  DO i = 1,n
    DO j = i+1,n ! j > i to avoid duplications
      r = pt(i)%p - pt(j)%p  ! We need the vector thar separate the particla i and j
      a(i) = a(i) + (pt(j)%m *r)/ distance(pt(i)%p,pt(j)%p)**3 !accel of particle i due to j
      a(j) = a(j) - (pt(i)%m *r)/ distance(pt(i)%p,pt(j)%p)**3 !accel of particle j due to i
    END DO
  END DO
 END SUBROUTINE calc_acc
 
END PROGRAM tree