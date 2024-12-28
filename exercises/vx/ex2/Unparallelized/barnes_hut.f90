module barneshut       !This module implements the Barnes-Hut algorithm.
  
  use iso_fortran_env  !This module ensures all variables are defined as 64-bit.
  use geometry         !This module defines 3D vector and point operations for vector3d and point3d types.
  use definitions      !This module defines key data structures for N-body simulations.
  implicit none

  INTEGER(int64) :: n
  REAL(real64) :: dt, t_end, t, dt_out, t_out, rs
  REAL(real64), PARAMETER :: theta = 50

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
  SUBROUTINE Calculate_Ranges(goal, p)
    TYPE(CELL),POINTER :: goal
    REAL(real64), DIMENSION(3) :: mins,maxs,medios
    REAL(real64) :: span
    INTEGER(INT64) :: i
    type(particle3d), allocatable :: p(:)
    
    REAL(real64), ALLOCATABLE, DIMENSION(:,:) :: positions
    ALLOCATE(positions(SIZE(p), 3))
    
    do i = 1, size(p)
       positions(i, :) = [p(i)%p%x, p(i)%p%y, p(i)%p%z]
    end do
    
    mins = MINVAL(positions,DIM=1)
    maxs = MAXVAL(positions,DIM=1)
    ! Al calcular span le sumo un 10% para que las
    ! particulas no caigan justo en el borde
    span = MAXVAL(maxs - mins) * 1.1
    medios = (maxs + mins) / 2.0
    goal%range%min = medios - span/2.0
    goal%range%max = medios + span/2.0
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
  RECURSIVE SUBROUTINE Find_Cell(root,goal,p)
    type(particle3d), intent(in) :: p
    TYPE(CELL),POINTER :: root,goal,temp
    INTEGER(int64) :: i,j,k
    SELECT CASE (root%type)
    CASE (2)
       out: DO i = 1,2
          DO j = 1,2
             DO k = 1,2
                IF (Belongs(p,root%subcell(i,j,k)%ptr)) THEN
                   CALL Find_Cell(root%subcell(i,j,k)%ptr,temp,p)
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
  RECURSIVE SUBROUTINE Place_Cell(goal,p,n)
    TYPE(CELL),POINTER :: goal,temp
    !REAL(real64), DIMENSION(3) :: part
    type(particle3d), intent(in) :: p
    INTEGER(int64) :: n
    SELECT CASE (goal%type)
    CASE (0)
       goal%type = 1
       goal%p = p
       goal%pos = n
    CASE (1)
       CALL Crear_Subcells(goal)
       CALL Find_Cell(goal,temp,p)
       CALL Place_Cell(temp,p,n)
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
    type(particle3d) :: p
    INTEGER(int64) :: i,j,k,n
    INTEGER(int64), DIMENSION(3) :: octant
    p = goal%p
    goal%type=2
    DO i = 1,2
       DO j = 1,2
          DO k = 1,2
             octant = (/i,j,k/)
             ALLOCATE(goal%subcell(i,j,k)%ptr)
             goal%subcell(i,j,k)%ptr%range%min = Calcular_Range (int(0, int64),goal,octant)
             goal%subcell(i,j,k)%ptr%range%max = Calcular_Range (int(1, int64),goal,octant)
             IF (Belongs(p,goal%subcell(i,j,k)%ptr)) THEN
                goal%subcell(i,j,k)%ptr%p = p
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
    INTEGER(int64) :: i,j,k
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
  FUNCTION Belongs (p,goal)
    type(particle3d) :: p
    TYPE(CELL), POINTER :: goal
    LOGICAL :: Belongs
    IF (p%p%x >= goal%range%min(1) .AND. &
         p%p%x < goal%range%max(1) .AND. &
         p%p%y >= goal%range%min(2) .AND. &
         p%p%y < goal%range%max(2) .AND. &
         p%p%z >= goal%range%min(3) .AND. &
         p%p%z < goal%range%max(3)) THEN
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
    INTEGER(int64) :: what,n
    TYPE(CELL), POINTER :: goal
    INTEGER(int64), DIMENSION(3) :: octant
    REAL(real64), DIMENSION(3) :: Calcular_Range, valor_medio
    valor_medio = (goal%range%min + goal%range%max) / 2.0
    SELECT CASE (what)
    CASE (0)
       WHERE (octant == 1)
          Calcular_Range = goal%range%min
       ELSEWHERE
          Calcular_Range = valor_medio
       ENDWHERE
    CASE (1)
       WHERE (octant == 1)
          Calcular_Range = valor_medio
       ELSEWHERE
          Calcular_Range = goal%range%max
       ENDWHERE
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
    INTEGER(int64) :: i,j,k
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
    INTEGER(int64) :: i,j,k
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
  RECURSIVE SUBROUTINE Calculate_masses(goal, p)
    TYPE(CELL),POINTER :: goal
    INTEGER(int64) :: i,j,k
    REAL(real64) :: mass
    type(point3d), allocatable :: c_o_m
    type(particle3d), intent(inout) :: p(:)
    goal%mass = 0
    goal%c_o_m%x = 0
    goal%c_o_m%y = 0
    goal%c_o_m%z = 0
    SELECT CASE (goal%type)
    CASE (1)
       goal%mass = goal%p%m
       goal%c_o_m%x = goal%p%p%x
       goal%c_o_m%y = goal%p%p%y
       goal%c_o_m%z = goal%p%p%z
    CASE (2)
       DO i = 1,2
          DO j = 1,2
             DO k = 1,2
                IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                   CALL Calculate_masses(goal%subcell(i,j,k)%ptr, p)
                   mass = goal%mass
                   goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
                   goal%c_o_m = DEBUGGERPRO9000B((mass * goal%c_o_m + &
                        goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m) / goal%mass)
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
  SUBROUTINE Calculate_forces(head, n, p, rji, a)
    TYPE(CELL),POINTER :: head
    integer(int64), intent(in) :: n
    INTEGER(int64) :: i,j,k,start,end
    type(particle3d), intent(inout) :: p(:)
    type(vector3d), intent(inout) :: rji
    type(vector3d), intent(inout) :: a(:)
    DO i = 1,n
       CALL Calculate_forces_aux(i,head, p, rji, a)
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
  RECURSIVE SUBROUTINE Calculate_forces_aux(goal,tree, p, rji, a)
    TYPE(CELL),POINTER :: tree
    INTEGER(int64) :: i,j,k,goal
    REAL(real64) :: l,D
    type(vector3d), intent(inout) :: rji
    real(real64) :: r2, r3
    type(particle3d), intent(inout) :: p(:)
    type(vector3d), intent(inout) :: a(:)

    
    SELECT CASE (tree%type)
    CASE (1)
       IF (goal .NE. tree%pos) THEN
          rji = tree%c_o_m - p(goal)%p
          r2 = rji%x**2 + rji%y**2 + rji%z**2
          r3 = r2 * SQRT(r2)
          a(goal) = a(goal) + (p(tree%pos)%m * rji / r3)
       END IF
    CASE (2)
       !! El rango tiene el mismo span en las 3 dimensiones
       !! por lo que podemos considerar una dimension cualquiera
       !! para calcular el lado de la celda (en este caso la
       !! dimension 1)
       l = tree%range%max(1) - tree%range%min(1)
       rji = tree%c_o_m - p(goal)%p
       r2 = rji%x**2 + rji%y**2 + rji%z**2
       D = SQRT(r2)
       IF (l/D < theta) THEN
          !! Si conglomerado, tenemos que ver si se cumple l/D < @
          r3 = r2 * D
          a(goal) = a(goal) + (tree%mass * rji / r3)
       ELSE
          DO i = 1,2
             DO j = 1,2
                DO k = 1,2
                   IF (ASSOCIATED(tree%subcell(i,j,k)%ptr)) THEN
                      CALL Calculate_forces_aux(goal,tree%subcell(i,j,k)%ptr, p, rji, a)
                   END IF
                END DO
             END DO
          END DO
       END IF
    END SELECT
  END SUBROUTINE Calculate_forces_aux
end module barneshut
