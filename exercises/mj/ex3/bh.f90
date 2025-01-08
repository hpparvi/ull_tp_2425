MODULE bh

  USE geometry
  USE particle
  
  IMPLICIT NONE
   

   TYPE RANGE
    REAL, DIMENSION(3) :: min,max
   END TYPE RANGE
  
   TYPE CPtr
    TYPE(CELL), POINTER :: ptr
   END TYPE CPtr
  
   TYPE CELL
    TYPE (RANGE) :: range
    TYPE(particle3d):: par
    INTEGER :: pos
    INTEGER :: type !! 0 = no particle; 1 = particle; 2 = conglomerado
    REAL :: mass
    TYPE(point3d) :: c_o_m
    TYPE (CPtr), DIMENSION(2,2,2) :: subcell
   END TYPE CELL
  
   TYPE (CELL), POINTER :: head, temp_cell
   
    CONTAINS 
        SUBROUTINE Calculate_Ranges(goal, par)
            TYPE(CELL),POINTER :: goal
            TYPE(particle3d), DIMENSION(:) :: par
            REAL, DIMENSION(3) :: mins,maxs,medios
            REAL :: span
        
        mins(1) = MINVAL(par%p%x)
        mins(2) = MINVAL(par%p%y)
		mins(3) = MINVAL(par%p%z)
		maxs(1) = MAXVAL(par%p%x)
		maxs(2) = MAXVAL(par%p%y)
		maxs(3) = MAXVAL(par%p%z)
        ! Al calcular span le sumo un 10% para que las
        ! particulas no caigan justo en el borde
        span = MAXVAL(maxs - mins) * 1.1
        medios = (maxs + mins) / 2.0
        goal%range%min = medios - span/2.0
        goal%range%max = medios + span/2.0
    END SUBROUTINE Calculate_Ranges
   
    RECURSIVE SUBROUTINE Find_Cell(root,goal,par)
      TYPE(particle3d):: par
      TYPE(CELL),POINTER :: root,goal,temp
      INTEGER :: i,j,k
      
      SELECT CASE (root%type)
        CASE (2)
          out: DO i = 1,2
            DO j = 1,2
              DO k = 1,2
                IF (Belongs(par,root%subcell(i,j,k)%ptr)) THEN
                  CALL Find_Cell(root%subcell(i,j,k)%ptr,temp,par)
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
    
    
    RECURSIVE SUBROUTINE Place_Cell(goal,par,n)
      TYPE(CELL),POINTER :: goal,temp
      TYPE(particle3d):: par
      INTEGER :: n
      SELECT CASE (goal%type)
        CASE (0)
          goal%type = 1
          goal%par = par
          goal%pos = n
        CASE (1)
          CALL Crear_Subcells(goal)
          CALL Find_Cell(goal,temp,par)
          CALL Place_Cell(temp,par,n)
        CASE DEFAULT
          print*,"SHOULD NOT BE HERE. ERROR!"
      END SELECT
    END SUBROUTINE Place_Cell
 
    
    SUBROUTINE Crear_Subcells(goal)
      TYPE(CELL), POINTER :: goal
      TYPE(particle3d):: par
      INTEGER :: i,j,k,n
      INTEGER, DIMENSION(3) :: octant
      par = goal%par
      goal%type=2
      
      DO i = 1,2
        DO j = 1,2
          DO k = 1,2
            octant = (/i,j,k/)
            ALLOCATE(goal%subcell(i,j,k)%ptr)
            goal%subcell(i,j,k)%ptr%range%min = Calcular_Range (0,goal,octant)
            goal%subcell(i,j,k)%ptr%range%max = Calcular_Range (1,goal,octant)
            IF (Belongs(par,goal%subcell(i,j,k)%ptr)) THEN
              goal%subcell(i,j,k)%ptr%par = par
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
    
    FUNCTION Belongs (par,goal)
     TYPE(particle3d):: par
     TYPE(CELL), POINTER :: goal
     LOGICAL :: Belongs
     IF (par%p%x >= goal%range%min(1) .AND. par%p%x < goal%range%max(1) .AND. &
		    par%p%y >= goal%range%min(2) .AND. par%p%y < goal%range%max(2) .AND. &
		    par%p%z >= goal%range%min(3) .AND. par%p%z < goal%range%max(3)) THEN
       Belongs = .TRUE.
     ELSE
       Belongs = .FALSE.
     END IF
    END FUNCTION Belongs
   
   
    FUNCTION Calcular_Range (what,goal,octant)
      INTEGER :: what,n
      TYPE(CELL), POINTER :: goal
      INTEGER, DIMENSION(3) :: octant
      REAL, DIMENSION(3) :: Calcular_Range, valor_medio
      
      valor_medio = (goal%range%min + goal%range%max) / 2
      
      SELECT CASE (what)
      CASE (0)
        WHERE (octant == 1)
          Calcular_Range = goal%range%min
        ELSEWHERE
          Calcular_Range = valor_medio
        END WHERE
      CASE (1)
        WHERE (octant == 1)
          Calcular_Range = valor_medio
        ELSEWHERE
          Calcular_Range = goal%range%max
        END WHERE
      END SELECT
    END FUNCTION Calcular_Range
    
  
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
    
  
    RECURSIVE SUBROUTINE Calculate_masses(goal, par)
      TYPE(CELL),POINTER :: goal
      INTEGER :: i,j,k
      REAL :: mass
      TYPE(point3d) :: c_o_m ! Center of mass
      TYPE(particle3d), DIMENSION(:) :: par
      
      goal%mass = 0
      goal%c_o_m = point3d(0,0,0)
      
      SELECT CASE (goal%type)
        CASE (1)
          goal%mass = par(goal%pos)%m
          goal%c_o_m = vector_to_point(par(goal%pos))%p
        CASE (2)
          DO i = 1,2
            DO j = 1,2
              DO k = 1,2
                IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                  CALL Calculate_masses(goal%subcell(i,j,k)%ptr, par)
                  mass = goal%mass
                  goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
                  goal%c_o_m%x = (mass * goal%c_o_m%x + goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m%x) / goal%mass
                  goal%c_o_m%y = (mass * goal%c_o_m%y + goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m%y) / goal%mass
                  goal%c_o_m%z = (mass * goal%c_o_m%z + goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m%z) / goal%mass
                END IF
              END DO
            END DO
        END DO
      END SELECT
    END SUBROUTINE Calculate_masses
    
   
    SUBROUTINE Calculate_forces(head, par, n,theta)
      TYPE(CELL),POINTER :: head
      INTEGER :: i,j,k,n
      REAL(kind=kind(1.0d0)), intent(in) :: theta
      TYPE(particle3d), DIMENSION(:) :: par
      
      
      n = size(par)
   
      DO i = 1,n
        CALL Calculate_forces_aux(i,head,par,theta)
      END DO
      
    END SUBROUTINE Calculate_forces
    
    
    RECURSIVE SUBROUTINE Calculate_forces_aux(goal,tree, par, theta)
      TYPE(CELL),POINTER :: tree
      INTEGER :: i,j,k,goal
      REAL :: l,D,r2,r3
      REAL(kind=kind(1.0d0)), intent(in) :: theta
      TYPE(particle3d), DIMENSION(:) :: par
      TYPE(vector3d)::rji
      
      SELECT CASE (tree%type)
        CASE (1)
          IF (goal .NE. tree%pos) THEN
            rji = point_to_vector(tree%c_o_m) - par(goal)%p
            r2 = rji%x**2+rji%y**2+rji%z**2
            r3 = r2 * sqrt(r2) + 10**(-7)
            par(goal)%a = par(goal)%a + par(tree%pos)%m * divvr(rji,REAL(r3, KIND=8))
          END IF
        CASE (2)
          l = tree%range%max(1) - tree%range%min(1)
          rji = point_to_vector(tree%c_o_m) - par(goal)%p
          r2 = rji%x**2+rji%y**2+rji%z**2
          D = sqrt(r2)
          IF (l/D < theta) THEN
            !! Si conglomerado, tenemos que ver si se cumple l/D < @
            r3 = r2 * D + 10**(-7)
            par(goal)%a = par(goal)%a + divvr(mulrv(REAL(tree%mass,KIND=8),rji),REAL(r3, KIND=8))

          ELSE
            DO i = 1,2
              DO j = 1,2
                DO k = 1,2
                  IF (ASSOCIATED(tree%subcell(i,j,k)%ptr)) THEN
                    CALL Calculate_forces_aux(goal,tree%subcell(i,j,k)%ptr,par,theta)
                  END IF
                END DO
              END DO
            END DO
          END IF
      END SELECT
    END SUBROUTINE Calculate_forces_aux
END MODULE bh
