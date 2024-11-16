PROGRAM tree
   IMPLICIT NONE
   INTEGER :: i,j,k,n
   REAL :: dt, t_end, t, dt_out, t_out, rs, r2, r3
   REAL, PARAMETER :: theta = 1
   REAL, DIMENSION(:), ALLOCATABLE :: m
   REAL, DIMENSION(:,:), ALLOCATABLE :: r,v,a
   REAL, DIMENSION(3) :: rji
   TYPE RANGE
      REAL, DIMENSION(3) :: min,max
   END TYPE RANGE
   TYPE CPtr
      TYPE(CELL), POINTER :: ptr
   END TYPE CPtr
   TYPE CELL
      TYPE (RANGE) :: range
      REAL, DIMENSION(3) :: part
      INTEGER :: pos
      INTEGER :: type !! 0 = no particle; 1 = particle; 2 = conglomerado
      REAL :: mass
      REAL, DIMENSION(3) :: c_o_m
      TYPE (CPtr), DIMENSION(2,2,2) :: subcell
   END TYPE CELL
   TYPE (CELL), POINTER :: head, temp_cell
   !! Lectura de datos
   !!!!!!!!!!!!!!!!!!!
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
   !! Inicialización head node
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   ALLOCATE(head)
   CALL Calculate_ranges(head)
   head%type = 0
   CALL Nullify_Pointers(head)
   !! Creación del árbol inicial
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   DO i = 1,n
      CALL Find_Cell(head,temp_cell,r(i,:))
      CALL Place_Cell(temp_cell,r(i,:),i)
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
   DO t = 0.0, t_end, dt
      v = v + a * dt/2
      r = r + v * dt
      !! Las posiciones han cambiado, por lo que tenemos que borrar
      !! y reinicializar el árbol
      CALL Borrar_tree(head)
      CALL Calculate_ranges(head)
      head%type = 0
      CALL Nullify_Pointers(head)
      DO i = 1,n
         CALL Find_Cell(head,temp_cell,r(i,:))
         CALL Place_Cell(temp_cell,r(i,:),i)
      END DO
      CALL Borrar_empty_leaves(head)
      CALL Calculate_masses(head)

      a = 0.0
      CALL Calculate_forces(head)
      v = v + a * dt/2
      t_out = t_out + dt
      IF (t_out >= dt_out) THEN
         DO i = 1,10
            PRINT*, r(i,:)
         END DO
         PRINT*, "-----------------------------------"
         PRINT*, ""
         t_out = 0.0
      END IF
   END DO
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_Ranges
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calcula los rangos de las partı́culas en la
!! matriz r en las 3 dimensiones y lo pone en la
!! variable apuntada por goal
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE Calculate_Ranges(goal)
      TYPE(CELL),POINTER :: goal
      REAL, DIMENSION(3) :: mins,maxs,medios
      REAL :: span
      mins = MINVAL(r,DIM=1)
      maxs = MAXVAL(r,DIM=1)
! Al calcular span le sumo un 10% para que las
! particulas no caigan justo en el borde
      span = MAXVAL(maxs - mins) * 1.1
      medios = (maxs + mins) / 2.0
      goal%range%min = medios - span/2.0
      goal%range%max = medios + span/2.0
   END SUBROUTINE Calculate_Ranges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Find_Cell
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Encuentra la celda donde colocaremos la particula.
!! Si la celda que estamos considerando no tiene
!! particula o tiene una particula, es esta celda donde
!! colocaremos la particula.
!! Si la celda que estamos considerando es un "conglomerado",
!! buscamos con la función BELONGS a que subcelda de las 8
!! posibles pertenece y con esta subcelda llamamos de nuevo
!! a Find_Cell
!!
!! NOTA: Cuando se crea una celda "conglomerado" se crean las
!! 8 subceldas, por lo que podemos asumir que siempre existen
!! las 8. Las celdas vacı́as se borran al final del todo, cuando
!! todo el árbol ha sido ya creado.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   RECURSIVE SUBROUTINE Find_Cell(root,goal,part)
      REAL, DIMENSION(3) :: part
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
!! Place_Cell
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Se ejecuta tras Find_Cell, en la celda que
!! esa función nos devuelve, por lo que siempre
!! es una celda de tipo 0 (sin particula) o de tipo 1
!! (con una particula). En el caso de que es una celda
!! de tipo 1 habra que subdividir la celda y poner en
!! su lugar las dos particulas (la que originalmente
!! estaba, y la nueva).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   RECURSIVE SUBROUTINE Place_Cell(goal,part,n)
      TYPE(CELL),POINTER :: goal,temp
      REAL, DIMENSION(3) :: part
      INTEGER :: n
      SELECT CASE (goal%type)
       CASE (0)
         goal%type = 1
         goal%part = part
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
!! Crear_Subcells
!!
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
      REAL,DIMENSION(3) :: part
      INTEGER :: i,j,k,n
      INTEGER, DIMENSION(3) :: octant
      part = goal%part
      goal%type=2
      DO i = 1,2
         DO j = 1,2
            DO k = 1,2
               octant = (/i,j,k/)
               ALLOCATE(goal%subcell(i,j,k)%ptr)
               goal%subcell(i,j,k)%ptr%range%min = Calcular_Range (0,goal,octant)
               goal%subcell(i,j,k)%ptr%range%max = Calcular_Range (1,goal,octant)

               IF (Belongs(part,goal%subcell(i,j,k)%ptr)) THEN
                  goal%subcell(i,j,k)%ptr%part = part
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
   !! Nullify_Pointers
   !!
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
   !! Belongs
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! Devuelve TRUE si la particula "part" está
   !! dentro del rango de la celda "goal"
   !!
   !! Utilizada por FIND_CELL
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   FUNCTION Belongs (part,goal)
      REAL, DIMENSION(3) :: part
      TYPE(CELL), POINTER :: goal
      LOGICAL :: Belongs
      IF (part(1) >= goal%range%min(1) .AND. &
         part(1) <= goal%range%max(1) .AND. &
         part(2) >= goal%range%min(2) .AND. &
         part(2) <= goal%range%max(2) .AND. &
         part(3) >= goal%range%min(3) .AND. &
         part(3) <= goal%range%max(3)) THEN
         Belongs = .TRUE.
      ELSE
         Belongs = .FALSE.
      END IF
   END FUNCTION Belongs
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Calcular_Range
   !!
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
      REAL, DIMENSION(3) :: Calcular_Range, valor_medio
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
   !! Borrar_empty_leaves
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! Se llama una vez completado el arbol para
   !! borrar (DEALLOCATE) las celdas vacı́as (i.e.
   !! sin partı́cula).
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
   !! Borrar_tree
   !!
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
!! Calculate_masses
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Nos calcula para todas las celdas que cuelgan
!! de "goal" su masa y su center-of-mass.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   RECURSIVE SUBROUTINE Calculate_masses(goal)
      TYPE(CELL),POINTER :: goal
      INTEGER :: i,j,k
      REAL :: mass
      REAL, DIMENSION(3) :: c_o_m
      goal%mass = 0
      goal%c_o_m = 0
      SELECT CASE (goal%type)
       CASE (1)
         goal%mass = m(goal%pos)
         goal%c_o_m = r(goal%pos,:)
       CASE (2)
         DO i = 1,2
            DO j = 1,2
               DO k = 1,2
                  IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                     CALL Calculate_masses(goal%subcell(i,j,k)%ptr)
                     mass = goal%mass
                     goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
                     goal%c_o_m = (mass * goal%c_o_m + &
                        goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m)
                  END IF
               END DO
            END DO
         END DO
      END SELECT
   END SUBROUTINE Calculate_masses
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_forces
!!
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
!! Calculate_forces_aux
!!
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
      REAL :: l,D
      SELECT CASE (tree%type)
       CASE (1)
         IF (goal .NE. tree%pos) THEN
            rji = tree%c_o_m - r(goal,:)
            r2 = SUM(rji**2)
            r3 = r2 * SQRT(r2)
            a(goal,:) = a(goal,:) + m(tree%pos) * rji / r3
         END IF
       CASE (2)
!! El rango tiene el mismo span en las 3 dimensiones
!! por lo que podemos considerar una dimension cualquiera
!! para calcular el lado de la celda (en este caso la
!! dimension 1)
         l = tree%range%max(1) - tree%range%min(1)
         rji = tree%c_o_m - r(goal,:)
         r2 = SUM(rji**2)
         D = SQRT(r2)
         IF (l/D < theta) THEN
!! Si conglomerado, tenemos que ver si se cumple l/D < @
            r3 = r2 * D
            a(goal,:) = a(goal,:) + tree%mass * rji / r3
         ELSE
            DO i = 1,2
               DO j = 1,2
                  DO k = 1,2
                     IF (ASSOCIATED(tree%subcell(i,j,k)%ptr)) THEN
                        CALL Calculate_forces_aux(goal,tree%subcell(i,j,k)%ptr)
                     END IF
                  END DO
               END DO
            END DO
         END IF
      END SELECT
   END SUBROUTINE Calculate_forces_aux
END PROGRAM tree
