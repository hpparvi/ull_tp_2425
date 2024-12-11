MODULE tree
  USE iso_fortran_env
  use mpi_f08
  USE geometry
  USE particle
  IMPLICIT NONE

  TYPE RANGE
     REAL(real64), DIMENSION(3) :: rmin, rmax
  END TYPE RANGE

  type CPtr
     type(CELL), pointer :: ptr
  end type CPtr

  TYPE CELL
     TYPE (RANGE) :: range
     type(particle3d) :: part
     INTEGER :: pos
     INTEGER :: TYPE !! 0 = no particle; 1 = particle; 2 = conglomerate
     type(point3d) :: center_of_mass
     real(real64) :: mass
     TYPE (CPtr), DIMENSION(2, 2, 2) :: subcell
  END TYPE CELL

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Calculate_Ranges !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Calculates the ranges of the particle positions
  !! in all 3 dimensions and then stores it in
  !! the variable "goal" is pointing to.
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calculate_Ranges(parts, goal)
    TYPE(CELL), POINTER, intent(inout) :: goal
    type(particle3d), dimension(:), intent(in) :: parts
    real(real64), dimension(3) :: mins, maxs, medios
    REAL(real64) :: span
    mins = (/minval(parts%p%xx, dim=1), minval(parts%p%yy, dim=1), &
         & minval(parts%p%zz, dim=1)/)

    maxs = (/maxval(parts%p%xx, dim=1), maxval(parts%p%yy, dim=1), &
         & maxval(parts%p%zz, dim=1)/)
    
    ! Al calcular span le sumo un 10% para que las
    ! particulas no caigan justo en el borde
    span = MAXVAL(maxs - mins) * 1.1
    medios = (maxs + mins) / 2.0
    goal%range%rmin = medios - span/2.0
    goal%range%rmax = medios + span/2.0
  END SUBROUTINE Calculate_Ranges

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Belongs !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Returns .TRUE. if the particle is inside
  !! the goal cell
  !!
  !! Used by Find_Cell
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical FUNCTION Belongs (part, goal)
    type(particle3d), intent(in) :: part
    TYPE(CELL), POINTER, intent(inout) :: goal

    IF (part%p%xx .GE. goal%range%rmin(1) .AND. &
         part%p%xx .LT. goal%range%rmax(1) .AND. &
         part%p%yy .GE. goal%range%rmin(2) .AND. &
         part%p%yy .LT. goal%range%rmax(2) .AND. &
         part%p%zz .GE. goal%range%rmin(3) .AND. &
         part%p%zz .LT. goal%range%rmax(3)) THEN
       Belongs = .TRUE.
    ELSE
       Belongs = .FALSE.
    END IF
  END FUNCTION Belongs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Find_Cell !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Finds the cell where the particle will go.
  !! If the cell under consideration is empty or has
  !! one particle, the particle goes there.
  !! If the cell under consideration is a conglomerate,
  !! we use the BELONGS function to find in which subcell
  !! the particle goes and calls Find_Cell on that subcell.
  !!
  !! NOTE: When a "conglomerate" (type 2) cell is created, 8
  !! subcells are created in it, so we can assume they all exist.
  !! The empty cells are deleted at the end, when the whole tree
  !! has been created.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE Find_Cell(root, goal, part)
    type(particle3d), intent(in) :: part
    TYPE(CELL), POINTER, intent(in) :: root
    TYPE(CELL), POINTER, intent(inout) :: goal
    TYPE(CELL), POINTER :: temp
    INTEGER :: i, j, k
    SELECT CASE (root%TYPE)
    CASE (2)
       out: DO i = 1, 2
          DO j = 1, 2
             DO k = 1, 2
                IF (Belongs(part, root%subcell(i,j,k)%ptr)) THEN
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
  !! Nullify_Pointers !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Nullifies the pointers for the 8 subcells
  !! in the "goal" cell.
  !!
  !! Used in the main program and by
  !! Create_Subcells.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Nullify_Pointers(goal)
    TYPE(CELL), POINTER :: goal
    INTEGER :: i,j,k
    DO i = 1, 2
       DO j = 1, 2
          DO k = 1, 2
             NULLIFY(goal%subcell(i, j, k)%ptr)
          END DO
       END DO
    END DO
  END SUBROUTINE Nullify_Pointers

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Calculate_OctRange !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Calculates the ranges for a given octant ((1,1,1),
  !! (1,1,2), ..., (2,2,2)) using the ranges given by
  !! "goal". If what = 0 it returns the minima, and
  !! for what = 1 it returns the maxima.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION Calculate_OctRange (what, goal, octant)
    INTEGER, intent(in) :: what
    TYPE(CELL), POINTER, intent(inout) :: goal
    INTEGER, DIMENSION(3), intent(in) :: octant
    REAL(real64), DIMENSION(3) :: Calculate_OctRange, valor_medio
    valor_medio = (goal%range%rmin + goal%range%rmax) / 2.0

    SELECT CASE (what)
    CASE (0)
       WHERE (octant == 1)
          Calculate_OctRange = goal%range%rmin
       ELSEWHERE
          Calculate_OctRange = valor_medio
       ENDWHERE
    CASE (1)
       WHERE (octant == 1)
          Calculate_OctRange = valor_medio
       ELSEWHERE
          Calculate_OctRange = goal%range%rmax
       ENDWHERE
    END SELECT
  END FUNCTION Calculate_OctRange

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Create_Subcells !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Calls Calculate_OctRange, Belongs, Nullify_Pointers
  !!
  !! Used when we are trying to place a particle
  !! into a cell that already has another in it.
  !! Creates 8 subcells under the goal cell,
  !! and places the particle inside goal into the
  !! corresponding subcell.
  !!
  !! Called by Place_Cell
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Create_Subcells(goal)
    TYPE(CELL), POINTER, intent(inout) :: goal
    type(particle3d) :: part
    INTEGER :: i, j, k
    INTEGER, DIMENSION(3) :: octant
    part = goal%part
    goal%TYPE=2
    DO i = 1,2
       DO j = 1,2
          DO k = 1,2
             octant = (/i,j,k/)
             ALLOCATE(goal%subcell(i,j,k)%ptr)
             goal%subcell(i,j,k)%ptr%range%rmin = Calculate_OctRange (0,goal,octant)
             goal%subcell(i,j,k)%ptr%range%rmax = Calculate_OctRange (1,goal,octant)
             IF (Belongs(part,goal%subcell(i,j,k)%ptr)) THEN
                goal%subcell(i,j,k)%ptr%part = part
                goal%subcell(i,j,k)%ptr%TYPE = 1
                goal%subcell(i,j,k)%ptr%pos = goal%pos
             ELSE
                goal%subcell(i,j,k)%ptr%TYPE = 0
             END IF
             CALL Nullify_Pointers(goal%subcell(i,j,k)%ptr)
          END DO
       END DO
    END DO
  END SUBROUTINE Create_Subcells

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Place_Cell !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Executed on the cell returned by Find_Cell,
  !! so a cell of type 0 (empty) or 1 (single particle).
  !! If it is a type 1 cell, it will be subdivided and
  !! each particle stored in the corresponding subcell.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE Place_Cell(goal, part, n)
    TYPE(CELL), POINTER, intent(inout) :: goal
    type(cell), pointer :: temp
    type(particle3d), intent(in) :: part
    INTEGER, intent(in) :: n
    SELECT CASE (goal%type)
       
    CASE (0)
       goal%TYPE = 1
       goal%part = part
       goal%pos = n
       
    CASE (1)
       CALL Create_Subcells(goal)
       CALL Find_Cell(goal, temp, part)
       CALL Place_Cell(temp, part, n)
       
    CASE DEFAULT
       PRINT*,"SHOULD NOT BE HERE. ERROR!"
    END SELECT
  END SUBROUTINE Place_Cell

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Delete_empty_leaves !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Called once the tree is completed in order
  !! to deallocate the empty cells.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE Delete_empty_leaves(goal)
    TYPE(CELL), POINTER, intent(inout) :: goal
    INTEGER :: i,j,k
    IF (ASSOCIATED(goal%subcell(1,1,1)%ptr)) THEN
       DO i = 1,2
          DO j = 1,2
             DO k = 1,2
                CALL Delete_empty_leaves(goal%subcell(i,j,k)%ptr)
                IF (goal%subcell(i,j,k)%ptr%TYPE == 0) THEN
                   DEALLOCATE (goal%subcell(i,j,k)%ptr)
                END IF
             END DO
          END DO
       END DO
    END IF
  END SUBROUTINE Delete_empty_leaves

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Delete_tree !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Deletes the whole tree, except the head.
  !!
  !! The tree must be regenerated continuously,
  !! so we must delete the old one in order to
  !! avoid memory leaks.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE Delete_tree(goal)
    TYPE(CELL), POINTER, intent(inout) :: goal
    INTEGER :: i, j, k
    DO i = 1, 2
       DO j = 1, 2
          DO k = 1, 2
             IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                CALL Delete_tree(goal%subcell(i,j,k)%ptr)
                DEALLOCATE (goal%subcell(i,j,k)%ptr)
             END IF
          END DO
       END DO
    END DO
  END SUBROUTINE Delete_tree

  !! And now with the physics:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Calculate_masses !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Calculates the mass and center of mass of each
  !! subcell inside the goal cell.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE Calculate_masses(goal)
    TYPE(CELL), POINTER, intent(inout) :: goal
    INTEGER :: i, j, k
    REAL(real64) :: mass
    type(vector3d) :: c_o_m
    type(point3d) :: zero = point3d(0.,0.,0.)
    goal%mass = 0.
    goal%center_of_mass = point3d(0.,0.,0.)
    SELECT CASE (goal%TYPE)
    CASE (1)
       goal%mass = goal%part%m
       goal%center_of_mass = goal%part%p
    CASE (2)
       DO i = 1, 2
          DO j = 1, 2
             DO k = 1, 2
                IF (ASSOCIATED(goal%subcell(i, j, k)%ptr)) THEN
                   CALL Calculate_masses(goal%subcell(i, j, k)%ptr)
                   mass = goal%mass
                   goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
                   
                   c_o_m = (mass * (goal%center_of_mass - zero) + &
                        goal%subcell(i,j,k)%ptr%mass * (goal%subcell(i,j,k)%ptr%center_of_mass - zero)) &
                        &/ goal%mass

                   goal%center_of_mass = point3d(c_o_m%xx, c_o_m%yy, c_o_m%zz) ! this is what we spanish call a chapuza
                END IF
             END DO
          END DO
       END DO
    END SELECT
  END SUBROUTINE Calculate_masses

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Calculate_forces_aux !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Given a goal particle (goal is the particle ID),
  !! calculates the force the "tree" cell exerts on it.
  !!
  !! If "tree" is a type 1 cell it's a two-particle
  !! case, which is simple.
  !!
  !! If "tree" is a type 2 cell, we have to check whether
  !! l/D < theta, that is, whether the side of the cell (l)
  !! divided by the distance from the goal particle to the
  !! center_of_mass of the tree cell (D) is less than theta.
  !! If it is, we treat the cell as a single particle. If it
  !! is not, we consider each subcell of the tree and call
  !! this function for each of them.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE Calculate_forces_aux(goal, tree, parts, aa, theta)
    TYPE(CELL), POINTER, intent(in) :: tree
    type(particle3d), intent(in) :: parts(:)
    integer, intent(in) :: goal
    type(vector3d), intent(inout) :: aa(:)
    real(real64), intent(in) :: theta
    
    INTEGER :: i, j, k
    REAL(real64) :: l, D, r2
    type(vector3d) :: rji
    
    SELECT CASE (tree%TYPE)
      
    CASE (1)
      
       IF (goal .NE. tree%pos) THEN !checks it's not the particle w itself
          rji = normalize(tree%part%p - parts(goal)%p)
          r2 = distance(parts(goal)%p, tree%part%p)**2
          aa(goal) = aa(goal) + tree%part%m * rji / r2
       END IF
       
    CASE (2)
       !! Span is the same in all 3 dimensions so we can use any
       !! to calculate l
       
       l = tree%range%rmax(1) - tree%range%rmin(1)
       rji = normalize(tree%center_of_mass - parts(goal)%p)
       r2 = distance(tree%center_of_mass, parts(goal)%p)**2
       D = distance(tree%center_of_mass, parts(goal)%p)
       IF (l/D < theta) THEN
          !! Si conglomerado, tenemos que ver si se cumple l/D < theta
          aa(goal) = aa(goal) + tree%mass * rji / r2
       ELSE
          DO i = 1,2
             DO j = 1,2
                DO k = 1,2
                   IF (ASSOCIATED(tree%subcell(i,j,k)%ptr)) THEN
                      CALL Calculate_forces_aux(goal,tree%subcell(i,j,k)%ptr, parts, aa, theta)
                   END IF
                END DO
             END DO
          END DO
       END IF
    END SELECT
  END SUBROUTINE Calculate_forces_aux


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Calculate_forces !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Calcula las fuerzas de todas las particulas contra "head".
  !! Se sirve de la funcion Calculate_forces_aux que es la
  !! que en realidad hace los calculos para cada particula
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calculate_forces(head, aa, parts, theta, istart, iend)
    TYPE(CELL), POINTER, intent(in) :: head
    type(vector3d), intent(inout) :: aa(:)
    type(particle3d), intent(in) :: parts(:)
    real(real64), intent(in) :: theta
    integer, intent(in) :: istart, iend
    INTEGER :: i, n
    n = size(parts)
    
    DO i = istart, iend
       CALL Calculate_forces_aux(i, head, parts, aa, theta)
    END DO

  END SUBROUTINE Calculate_forces

  
END MODULE tree
