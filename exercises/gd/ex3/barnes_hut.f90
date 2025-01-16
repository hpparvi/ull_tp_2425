MODULE barnes_hut

  USE geometry
  USE ex3_deriv_types
  USE ISO_FORTRAN_ENV
  !$ use omp_lib
  IMPLICIT NONE
  INTEGER :: i,j,k,n
  INTEGER :: rank, root
  REAL(REAL64), PARAMETER :: theta = 1
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: pt   !Particles
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a    !Accelerations
  TYPE(vector3d) :: r
  TYPE (CELL), POINTER :: head, temp_cell
  
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_Ranges !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calculates the ranges of the particles in the
!! matrix `r` across the 3 dimensions and stores them in the
!! variable pointed to by `goal`.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Calculate_Ranges(goal)
    TYPE(CELL),POINTER :: goal
    TYPE(point3d) :: mins,maxs,medios
    REAL(REAL64) :: span
      ! The minimum and maximum are calculated coordinate by coordinate.
    mins%x = MINVAL(pt%p%x)
    mins%y = MINVAL(pt%p%y)
    mins%z = MINVAL(pt%p%z)
    maxs%x = MAXVAL(pt%p%x)
    maxs%y = MAXVAL(pt%p%y)
    maxs%z = MAXVAL(pt%p%z)
      ! When calculating the span, it add 10% so that the
      ! particles do not fall exactly on the edge.
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
!! Finds the cell where we will place the particle.
!! If the cell we are considering does not have a particle
!! or has a particle, it is this cell where we will place
!! the particle.
!! If the cell we are considering is a "cluster",
!! we use the BELONGS function to find which of the 8
!! possible subcells it belongs to, and with this subcell,
!! we call Find_Cell again.
!!
!! NOTE: When a "cluster" cell is created, the 8 subcells are
!! also created, so we can assume that all 8 always exist.
!! Empty cells are deleted at the very end, when
!! the entire tree has already been created.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  RECURSIVE SUBROUTINE Find_Cell(root,goal,part)
  
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
!! It is executed after Find_Cell, in the cell that
!! the function returns, so it is always a cell of type 0
!! (without a particle) or of type 1 (with a particle).
!! In the case of a type 1 cell, the cell must be subdivided
!! and both particles (the one that was originally there,
!! and the new one) must be placed in it.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE Place_Cell(goal,part,n)
    TYPE(CELL),POINTER :: goal,temp
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
        print*,"SHOULD NOT BE HERE. ERROR!", rank
      END SELECT
  END SUBROUTINE Place_Cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Crear_Subcells !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This function is called from Place_Cell and
!! it is only called when there is already a particle
!! in the cell, so it needs to be subdivided.
!! What it does is create 8 subcells
!! that "hang" from goal, and the particle that
!! was in goal is placed in the corresponding subcell
!! among the 8 newly created ones.
!!
!! To create the subcells, use the functions
!! CALCULATE_RANGE, BELONGS, and NULLIFY_POINTERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Crear_Subcells(goal)
    TYPE(CELL), POINTER :: goal
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
!! It simply NULLIFYes the pointers of
!! the 8 subcells of the "goal" cell.
!!
!! It is used in the main loop and by
!! CREATE_SUBCELLS.
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
!! Returns TRUE if the particle "part" is
!! within the range of the "goal" cell.
!!
!! Used by FIND_CELL.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION Belongs (part,goal)
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
!! Given an octant (1,1,1, 1,1,2 ... 2,2,2),
!! it calculates its ranges based on the ranges of
!! "goal". If "what" = 0, it calculates the minimums.
!! If "what" = 1, it calculates the maximums.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION Calcular_Range (what,goal,octant)
    INTEGER :: what,n
    TYPE(CELL), POINTER :: goal
    INTEGER, DIMENSION(3) :: octant
    TYPE(point3d) :: Calcular_Range, valor_medio
    !! The average, maximum, and minimum values are calculated coordinate by coordinate.    
    valor_medio%x = (goal%range%min%x + goal%range%max%x) / 2.0
    valor_medio%y = (goal%range%min%y + goal%range%max%y) / 2.0
    valor_medio%z = (goal%range%min%z + goal%range%max%z) / 2.0
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
!! It is called once the tree is completed to
!! DEALLOCATE the empty cells (i.e.
!! cells without a particle).
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
!! Deletes the entire tree, except for the "head".
!!
!! The tree needs to be regenerated continuously,
!! so we must delete the old one
!! to avoid "memory leaks".
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
!! It calculates for all the cells hanging
!! from "goal" their mass and their center-of-mass.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE Calculate_masses(goal)
    TYPE(CELL),POINTER :: goal
    INTEGER :: i,j,k
    REAL(REAL64) :: mass
    !! We treat the center of mass as the vector that connects the point (0,0,0) 
    !! with the position of the center of mass.
    TYPE(vector3d) :: c_o_m
    TYPE(point3d) :: V3d_0 = point3d(0,0,0)

    goal%pt%m = 0
    goal%c_o_m = vector3d(0,0,0)
    SELECT CASE (goal%type)
    CASE (1)
      goal%pt%m = pt(goal%pos)%m
      goal%c_o_m = V3d_0 - pt(goal%pos)%p
    CASE (2)
      DO i = 1,2
        DO j = 1,2
          DO k = 1,2
            IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
              CALL Calculate_masses(goal%subcell(i,j,k)%ptr)
              mass = goal%pt%m
              goal%pt%m = goal%pt%m + goal%subcell(i,j,k)%ptr%pt%m
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
!! Calculates the forces of all the particles against "head".
!! It uses the function Calculate_forces_aux, which is the
!! one that actually performs the calculations for each particle.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Calculate_forces(head,pt,a,i_start, i_end)
    TYPE(CELL),POINTER, intent(in) :: head
    type(vector3d), intent(inout) :: a(:)
    INTEGER :: i,j,k,i_start, i_end
    type(particle3d), intent(in) :: pt(:)
    !!  Process does not execute when there are no particles in the node
    IF(i_end /= 0) THEN
      DO i = i_start, i_end
        CALL Calculate_forces_aux(i,pt,a,head)
      END DO
    END IF
    
  END SUBROUTINE Calculate_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_forces_aux !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Given a particle "goal", it calculates the forces
!! on it from the "tree" cell. If "tree" is a
!! cell containing a single particle, the case is simple
!! because there are only two particles involved.
!!
!! If "tree" is a cluster cell, we first need to check
!! if l/D < theta. That is, if the cell side length (l)
!! divided by the distance from the "goal" particle to the
!! center_of_mass of the "tree" cell (D) is smaller than theta.
!! If this is the case, we treat the cell as a single particle.
!! If it is not smaller than theta, then we must consider all
!! subcells of "tree" and recursively call Calculate_forces_aux
!! for each of them.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE Calculate_forces_aux(goal,pt,a,tree)
    type(particle3d), intent(in) :: pt(:)
    type(vector3d), intent(inout) :: a(:)
    TYPE(CELL),POINTER, intent(in) :: tree
    INTEGER :: i,j,k,goal
    REAL(REAL64) :: l,D
    !! For the forces, we treat the center of mass as a point.    
    TYPE(vector3d) :: r
    TYPE(point3d) :: c_o_m
    TYPE(point3d) :: V3d_0 = point3d(0,0,0)
    SELECT CASE (tree%type)
    CASE (1)
      IF (goal .NE. tree%pos) THEN
        c_o_m =  V3d_0 + tree%c_o_m 
        r = pt(goal)%p - c_o_m 
        a(goal) = a(goal) + (pt(tree%pos)%m *r)/ distance(pt(goal)%p,c_o_m)**3 !accel of particle i due to j
      END IF
    CASE (2)
      c_o_m = V3d_0 + tree%c_o_m 
      l = tree%range%max%x - tree%range%min%x
      D = distance(pt(goal)%p,c_o_m)
      IF (l/D < theta) THEN
        r = pt(goal)%p - c_o_m 
        a(goal) = a(goal) + (tree%pt%m *r)/ distance(pt(goal)%p,c_o_m)**3 !accel of particle i due to j
      ELSE
        DO i = 1,2
          DO j = 1,2
            DO k = 1,2
              IF (ASSOCIATED(tree%subcell(i,j,k)%ptr)) THEN
                CALL Calculate_forces_aux(goal,pt,a,tree%subcell(i,j,k)%ptr)
              END IF
            END DO
          END DO
        END DO
      END IF
    END SELECT
  END SUBROUTINE Calculate_forces_aux
 
END MODULE barnes_hut
