MODULE barnes
  USE geometry !Importing the geometry module
  USE particle !Importing the particle module
  USE iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  IMPLICIT NONE

  !Define integer variables for loop indexing
  INTEGER(INT64) :: i, j, k

  !Define real variables for squared and cubed distances
  REAL(REAL64) :: r2, r3

  REAL(REAL64), PARAMETER :: theta = 1

  TYPE RANGE
     REAL(REAL64), DIMENSION(3) :: min, max
  END TYPE RANGE
  
  TYPE CPtr
     TYPE(CELL), POINTER :: ptr
  END TYPE CPtr
  
  TYPE CELL
     TYPE(RANGE) :: range
     TYPE(particle3d) :: part !Particles
     INTEGER(INT64) :: pos
     INTEGER(INT64) :: type !! 0 = no particle; 1 = particle; 2 = conglomerado
     TYPE(vector3d) :: c_o_m
     TYPE (CPtr), DIMENSION(2, 2, 2) :: subcell
  END TYPE CELL

CONTAINS
  
  SUBROUTINE Calculate_Ranges(goal, p)
    TYPE(CELL), POINTER :: goal
    REAL(REAL64), DIMENSION(3) :: mins, maxs, medios
    REAL(REAL64) :: span
    TYPE(particle3d), INTENT(IN) :: p(:) !Particles

    
    mins = [MINVAL([p(:)%p%x]), MINVAL([p(:)%p%y]), MINVAL([p(:)%p%z])]
    maxs = [MAXVAL([p(:)%p%x]), MAXVAL([p(:)%p%y]), MAXVAL([p(:)%p%z])]
    span = MAXVAL(maxs - mins) * 1.1
    
    medios = (maxs + mins) * 0.5
    goal%range%min = medios - span * 0.5
    goal%range%max = medios + span * 0.5
    
  END SUBROUTINE Calculate_Ranges

  RECURSIVE SUBROUTINE Find_Cell(root, goal, part)
    TYPE(CELL), POINTER :: root, goal, temp
    TYPE(particle3d), INTENT(INOUT) :: part !Particles

   
    SELECT CASE (root%type)
    CASE (2)
       out: DO i = 1, 2
          DO j = 1, 2
             DO k = 1, 2
                IF (Belongs(part, root%subcell(i, j, k)%ptr)) THEN
                   
                   CALL Find_Cell(root%subcell(i, j, k)%ptr, temp, part)
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

  RECURSIVE SUBROUTINE Place_Cell(goal, part, n)
    TYPE(CELL), POINTER :: goal, temp
    TYPE(particle3d), INTENT(INOUT) :: part !Particles
    INTEGER(INT64), INTENT(IN) :: n !Number of bodies
    
    SELECT CASE (goal%type)
    CASE (0)
       goal%type = 1
       goal%part = part
       goal%pos = n
       
    CASE (1)
       CALL Crear_Subcells(goal, part, n)
       CALL Find_Cell(goal, temp, part)
       CALL Place_Cell(temp, part, n)
       
    CASE DEFAULT
       print*,"SHOULD NOT BE HERE. ERROR!"
    END SELECT
    
  END SUBROUTINE Place_Cell

  SUBROUTINE Crear_Subcells(goal, part, n)
    TYPE(CELL), POINTER :: goal
    TYPE(particle3d), INTENT(INOUT) :: part !Particles
    INTEGER(INT64), INTENT(IN) :: n !Number of bodies
    INTEGER, DIMENSION(3) :: octant
    
    part = goal%part
    goal%type = 2
    
    DO i = 1, 2
       DO j = 1, 2
          DO k = 1, 2
             octant = (/i, j, k/)
             ALLOCATE(goal%subcell(i, j, k)%ptr)
             
             goal%subcell(i, j, k)%ptr%range%min = Calcular_Range(0, goal, octant)
             goal%subcell(i, j, k)%ptr%range%max = Calcular_Range(1 ,goal, octant)
             
             IF (Belongs(part, goal%subcell(i, j, k)%ptr)) THEN
                goal%subcell(i, j, k)%ptr%part = part
                goal%subcell(i, j, k)%ptr%type = 1
                goal%subcell(i, j, k)%ptr%pos = goal%pos
                
             ELSE
                goal%subcell(i, j, k)%ptr%type = 0
             END IF
             
             CALL Nullify_Pointers(goal%subcell(i, j, k)%ptr)
             
          END DO
       END DO
    END DO
  END SUBROUTINE Crear_Subcells

  SUBROUTINE Nullify_Pointers(goal)
    TYPE(CELL), POINTER :: goal

    DO i = 1, 2
       DO j = 1, 2
          DO k = 1, 2
             NULLIFY(goal%subcell(i, j, k)%ptr)
          END DO
       END DO
    END DO
  END SUBROUTINE Nullify_Pointers

  FUNCTION Belongs (part, goal)
    TYPE(particle3d), INTENT(INOUT) :: part !Particles
    TYPE(CELL), POINTER :: goal
    LOGICAL :: Belongs
    
    IF ((part%p%x >= goal%range%min(1) .AND. part%p%x <= goal%range%max(1)) .AND. &
         (part%p%y >= goal%range%min(2) .AND. part%p%y <= goal%range%max(2)) .AND. &
         (part%p%z >= goal%range%min(3) .AND. part%p%z <= goal%range%max(3))) THEN
       Belongs = .TRUE.
    ELSE
       Belongs = .FALSE.
    END IF

  END FUNCTION Belongs

  FUNCTION Calcular_Range (what, goal, octant)
    INTEGER(INT32) :: what
    TYPE(CELL), POINTER :: goal
    INTEGER(INT32), DIMENSION(3) :: octant
    REAL(REAL64), DIMENSION(3) :: Calcular_Range, valor_medio
    
    valor_medio = (goal%range%min + goal%range%max) * 0.5
    
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

  RECURSIVE SUBROUTINE Borrar_empty_leaves(goal)
    TYPE(CELL), POINTER :: goal

    IF (ASSOCIATED(goal%subcell(1, 1, 1)%ptr)) THEN
       DO i = 1, 2
          DO j = 1, 2
             DO k = 1 ,2
                CALL Borrar_empty_leaves(goal%subcell(i, j, k)%ptr)
                
                IF (goal%subcell(i, j, k)%ptr%type == 0) THEN
                   DEALLOCATE (goal%subcell(i, j, k)%ptr)
                END IF
             END DO
          END DO
       END DO
    END IF
    
  END SUBROUTINE Borrar_empty_leaves

  RECURSIVE SUBROUTINE Borrar_tree(goal)
    TYPE(CELL), POINTER :: goal

    DO i = 1, 2
       DO j = 1, 2
          DO k = 1, 2
             IF (ASSOCIATED(goal%subcell(i, j, k)%ptr)) THEN
                
                CALL Borrar_tree(goal%subcell(i, j, k)%ptr)
                DEALLOCATE (goal%subcell(i, j, k)%ptr)
                
             END IF
          END DO
       END DO
    END DO
  END SUBROUTINE Borrar_tree

  RECURSIVE SUBROUTINE Calculate_masses(goal, p)
    TYPE(CELL), POINTER :: goal
    TYPE(vector3d) :: c_o_m
    TYPE(particle3d), INTENT(INOUT) :: p(:) !Particles
    
    goal%part%m = 0
    goal%c_o_m = vector3d(0, 0, 0)
    
    SELECT CASE (goal%type)
    CASE (1)
       goal%part%m = p(goal%pos)%m
       goal%c_o_m = point_to_vector(p(goal%pos)%p)
       
    CASE (2)
       DO i = 1, 2
          DO j = 1, 2
             DO k = 1, 2
                IF (ASSOCIATED(goal%subcell(i, j, k)%ptr)) THEN
                   CALL Calculate_masses(goal%subcell(i, j, k)%ptr, p)
                   
                   goal%part%m = goal%part%m + goal%subcell(i, j, k)%ptr%part%m
                   goal%c_o_m = (goal%part%m * goal%c_o_m + &
                        goal%subcell(i, j, k)%ptr%part%m * goal%subcell(i, j, k)%ptr%c_o_m)/goal%part%m
                END IF
             END DO
          END DO
       END DO
    END SELECT
  END SUBROUTINE Calculate_masses

  SUBROUTINE Calculate_forces(head, n, p, rji)
    TYPE(CELL), POINTER :: head
    INTEGER(INT64), INTENT(IN) :: n !Number of bodies
    TYPE(particle3d), INTENT(INOUT) :: p(:) !Particles
    TYPE(vector3d), INTENT(INOUT) :: rji !Vector from one particle to another
    
    DO i = 1, n
       CALL Calculate_forces_aux(i, head, p, rji)
    END DO
  END SUBROUTINE Calculate_forces

  RECURSIVE SUBROUTINE Calculate_forces_aux(goal, tree, p, rji)
    INTEGER(INT64) :: goal
    TYPE(CELL), POINTER :: tree
    TYPE(particle3d), INTENT(INOUT) :: p(:) !Particles
    TYPE(vector3d), INTENT(INOUT) :: rji !Vector from one particle to another
    REAL(REAL64) :: l, D
    
    SELECT CASE (tree%type)
    CASE (1)
       IF (goal .NE. tree%pos) THEN
          rji = tree%c_o_m - point_to_vector(p(goal)%p)
          r2 = (NORM(rji))**2  !Square of the vector
          r3 = r2 * SQRT(r2) !Cube of the vector
          p(goal)%a = p(goal)%a + p(tree%pos)%m * rji / r3
       END IF
       
    CASE (2)
       l = tree%range%max(1) - tree%range%min(1)
       rji = tree%c_o_m - point_to_vector(p(goal)%p)
       D = NORM(rji)
       IF (l/D < theta) THEN
          r3 = D**3
          p(goal)%a = p(goal)%a + tree%part%m * rji / r3
          
       ELSE
          DO i = 1, 2
             DO j = 1, 2
                DO k = 1, 2
                   IF (ASSOCIATED(tree%subcell(i, j, k)%ptr)) THEN
                      CALL Calculate_forces_aux(goal, tree%subcell(i, j, k)%ptr, p, rji)
                   END IF
                END DO
             END DO
          END DO
       END IF
    END SELECT
  END SUBROUTINE Calculate_forces_aux

END MODULE barnes
