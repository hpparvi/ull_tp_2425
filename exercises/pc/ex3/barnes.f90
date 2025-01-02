MODULE barnes
  USE geometry !Importing the geometry module
  USE particle !Importing the particle module
  USE iso_fortran_env !Importing iso_fortran_env to specify the number of bits for variables 
  IMPLICIT NONE

  !Define a RANGE type to hold the minimum and maximum range in 3D space
  TYPE RANGE
     TYPE(vector3d) :: min, max
  END TYPE RANGE
  
  !Define a CPtr type to hold a pointer to a CELL type
  TYPE CPtr
     TYPE(CELL), POINTER :: ptr
  END TYPE CPtr
  
  !Define a CELL type to store data about a cell
  TYPE CELL
     TYPE(RANGE) :: range
     TYPE(particle3d) :: part !Particle information for the cell
     INTEGER(INT64) :: pos, type !In type: 0 = no particle; 1 = particle; 2 = conglomerate (group of particles)
     TYPE(vector3d) :: c_o_m !Center of mass of the particles in the cell
     REAL(REAL64) :: mass !Mass of the particles in the cell
     TYPE(CPtr), DIMENSION(2, 2, 2) :: subcell !Subcells of the current cell
  END TYPE CELL

CONTAINS
  
  !Subroutine to calculate the ranges of the particles in the 3D position matrix and store them in goal
  SUBROUTINE Calculate_Ranges(goal, p)
    TYPE(CELL), POINTER :: goal !Target cell
    TYPE(vector3d) :: mins, maxs, medios !Vectors to store the minimum, maximum and mean values
    REAL(REAL64) :: span !Range of the particle positions
    TYPE(particle3d), INTENT(IN) :: p(:) !Particles
    
    !Calculate the minimum and maximum values for x, y, z coordinates of the particles
    mins = vector3d(MINVAL([p(:)%p%x]), MINVAL([p(:)%p%y]), MINVAL([p(:)%p%z]))
    maxs = vector3d(MAXVAL([p(:)%p%x]), MAXVAL([p(:)%p%y]), MAXVAL([p(:)%p%z]))
    
    !Calculate the span
    span = MAXVAL([maxs%x - mins%x, maxs%y - mins%y, maxs%z - mins%z]) * 1.1

    !Calculate the mean value for x, y, z coordinates of the particles
    medios = (maxs + mins) * REAL(0.5, REAL64)
   
    !Set the minimum and maximum ranges for the goal cell
    goal%range%min = medios - vector3d(span, span, span) * REAL(0.5, REAL64)
    goal%range%max = medios + vector3d(span, span, span) * REAL(0.5, REAL64)

  END SUBROUTINE Calculate_Ranges

  !Recursive subroutine to find the appropriate subcell for placing the particle
  RECURSIVE SUBROUTINE Find_Cell(root, goal, part)
    TYPE(CELL), POINTER :: root, goal, temp !Root node of the tree, target and temporary cell
    TYPE(particle3d), INTENT(IN) :: part !Particle
    INTEGER(INT64) :: i, j, k !Loop indexing variables

    SELECT CASE (root%type)
    CASE (2) !Conglomerate cell
       out: DO i = 1, 2
          DO j = 1, 2
             DO k = 1, 2
                IF (Belongs(part, root%subcell(i, j, k)%ptr)) THEN !Check which subcell the particle belongs to
                   CALL Find_Cell(root%subcell(i, j, k)%ptr, temp, part) !Recurse into the subcell
                   goal => temp !Set target cell to temporary cell
                   EXIT out
                END IF
             END DO
          END DO
       END DO out
       
    CASE DEFAULT !No particle or one particle in the cell
       goal => root !Set the target cell to root node
    END SELECT
    
  END SUBROUTINE Find_Cell

  !Recursive subroutine to place the particle into the appropriate cell
  RECURSIVE SUBROUTINE Place_Cell(goal, part, n)
    TYPE(CELL), POINTER :: goal, temp !Target and temporary cell
    TYPE(particle3d), INTENT(IN) :: part !Particle
    INTEGER(INT64) :: n !Number of bodies
    
    SELECT CASE (goal%type)
    CASE (0) !No particle in the cell
       goal%type = 1
       goal%part = part
       goal%pos = n
       
    CASE (1) !One particle in the cell
       !Subdivide the cell to fit the particle
       CALL Crear_Subcells(goal)
       CALL Find_Cell(goal, temp, part)
       CALL Place_Cell(temp, part, n)
       
    CASE DEFAULT
       PRINT*, "ERROR: Unexpected condition in Place_Cell!"
    END SELECT
    
  END SUBROUTINE Place_Cell

  !Subroutine to subdivide the cell when a particle is placed into it
  SUBROUTINE Crear_Subcells(goal)
    TYPE(CELL), POINTER :: goal !Target cell
    TYPE(particle3d) :: part !Particle
    INTEGER(INT64), DIMENSION(3) :: octant
    INTEGER(INT64) :: i, j, k !Loop indexing variables
    
    part = goal%part
    goal%type = 2 !Mark the cell as a conglomerate
    
    !Loop over the subcells and initialize them
    DO i = 1, 2
       DO j = 1, 2
          DO k = 1, 2
             octant = (/i, j, k/)
             ALLOCATE(goal%subcell(i, j, k)%ptr)
             
             !Set the range for each subcell
             goal%subcell(i, j, k)%ptr%range%min = Calcular_Range(INT(0, INT64), goal, octant)
             goal%subcell(i, j, k)%ptr%range%max = Calcular_Range(INT(1, INT64), goal, octant)
             
             !Check if the particle belongs to this subcell and assign it
             IF (Belongs(part, goal%subcell(i, j, k)%ptr)) THEN
                goal%subcell(i, j, k)%ptr%part = part
                goal%subcell(i, j, k)%ptr%type = 1
                goal%subcell(i, j, k)%ptr%pos = goal%pos
                
             ELSE
                goal%subcell(i, j, k)%ptr%type = 0 !Empty subcell
             END IF
             
             !Nullify pointers for clean-up
             CALL Nullify_Pointers(goal%subcell(i, j, k)%ptr)
          END DO
       END DO
    END DO
    
  END SUBROUTINE Crear_Subcells

  !Subroutine to nullify all pointers in the subcells
  SUBROUTINE Nullify_Pointers(goal) 
    TYPE(CELL), POINTER :: goal !Target cell
    INTEGER(INT64) :: i, j, k !Loop indexing variables

    DO i = 1, 2
       DO j = 1, 2
          DO k = 1, 2
             NULLIFY(goal%subcell(i, j, k)%ptr) !Nullify each pointer in subcells
          END DO
       END DO
    END DO
    
  END SUBROUTINE Nullify_Pointers

  !Function to check if a particle belongs to a given cell based on its range
  FUNCTION Belongs (part, goal)
    TYPE(particle3d) :: part !Particle
    TYPE(CELL), POINTER :: goal !Target cell
    LOGICAL :: Belongs
    
    !Check if the particle is within the bounds of the cell's range
    IF ((part%p%x .GE. goal%range%min%x) .AND. &
        (part%p%x .LT. goal%range%max%x) .AND. &
        (part%p%y .GE. goal%range%min%y) .AND. &
        (part%p%y .LT. goal%range%max%y) .AND. &
        (part%p%z .GE. goal%range%min%z) .AND. &
        (part%p%z .LT. goal%range%max%z)) THEN
       Belongs = .TRUE.
       
    ELSE
       Belongs = .FALSE.
    END IF

  END FUNCTION Belongs

  !Function to calculate the range for a subcell
  FUNCTION Calcular_Range(what, goal, octant)
    INTEGER(INT64) :: what !Integer variable that determines the range for a subcell
    TYPE(CELL), POINTER :: goal !Target cell
    INTEGER(INT64), DIMENSION(3) :: octant !Subcells
    TYPE(vector3d) :: Calcular_Range, valor_medio

    !Calculate the mean value of the range
    valor_medio = (goal%range%min + goal%range%max) * REAL(0.5, REAL64)
    
    SELECT CASE (what)
    CASE (0) !Calculate the range for a subcell with either the minimum or the mean value
       !If octant == 1, Calcular_Range = goal%range%min (minimum). If not, Calcular_Range = valor_medio (mean)
       Calcular_Range%x = MERGE(goal%range%min%x, valor_medio%x, octant(1) == 1)
       Calcular_Range%y = MERGE(goal%range%min%y, valor_medio%y, octant(2) == 1)
       Calcular_Range%z = MERGE(goal%range%min%z, valor_medio%z, octant(3) == 1)
       
    CASE (1) !Calculate the range for a subcell with either the maximum or the mean value
        !If octant == 1, Calcular_Range = goal%range%max (maximum). If not, Calcular_Range = valor_medio (mean)
       Calcular_Range%x = MERGE(valor_medio%x, goal%range%max%x, octant(1) == 1)
       Calcular_Range%y = MERGE(valor_medio%y, goal%range%max%y, octant(2) == 1)
       Calcular_Range%z = MERGE(valor_medio%z, goal%range%max%z, octant(3) == 1)
       
    END SELECT
    
  END FUNCTION Calcular_Range

  !Recursive subroutine to remove empty subcells in the hierarchical cell structure
  RECURSIVE SUBROUTINE Borrar_empty_leaves(goal)
    TYPE(CELL), POINTER :: goal !Target cell
    INTEGER(INT64) :: i, j, k !Loop indexing variables

    !Check if the subcells exist in the current cell
    IF (ASSOCIATED(goal%subcell(1, 1, 1)%ptr)) THEN
       !Loop through all 8 subcells
       DO i = 1, 2
          DO j = 1, 2
             DO k = 1, 2
                CALL Borrar_empty_leaves(goal%subcell(i, j, k)%ptr) !Recursively call for subcells
                !If a subcell has no particles, deallocate it
                IF (goal%subcell(i, j, k)%ptr%type == 0) THEN
                   DEALLOCATE(goal%subcell(i, j, k)%ptr) !Deallocate empty subcells
                END IF
             END DO
          END DO
       END DO
    END IF
    
  END SUBROUTINE Borrar_empty_leaves

  !Recursive subroutine to delete the entire tree of subcells and their contents
  RECURSIVE SUBROUTINE Borrar_tree(goal)
    TYPE(CELL), POINTER :: goal
    INTEGER(INT64) :: i, j, k !Loop indexing variables

    !Loop through all 8 subcells
    DO i = 1, 2
       DO j = 1, 2
          DO k = 1, 2
             !If the subcell exists, recursively call for deletion of subcells
             IF (ASSOCIATED(goal%subcell(i, j, k)%ptr)) THEN
                CALL Borrar_tree(goal%subcell(i, j, k)%ptr) !Recursively call for subcells
                DEALLOCATE(goal%subcell(i, j, k)%ptr) !Deallocate the subcell
             END IF
          END DO
       END DO
    END DO
    
  END SUBROUTINE Borrar_tree

  !Recursive subroutine to calculate the mass and center of mass for each cell
  RECURSIVE SUBROUTINE Calculate_masses(goal, p)
    TYPE(CELL), POINTER :: goal !Target cell
    TYPE(vector3d) :: c_o_m !Center of mass for the cell
    TYPE(particle3d), INTENT(INOUT) :: p(:) !Particles
    REAL(REAL64) :: mass !Temporary mass variable 
    INTEGER(INT64) :: i, j, k !Loop indexing variables

    !Initialize mass and center of mass
    goal%mass = 0
    goal%c_o_m = vector3d(0, 0, 0)

    SELECT CASE (goal%type)
    CASE (1) !One particle in the cell
       !Set the mass and the center of mass
       goal%mass = p(goal%pos)%m
       goal%c_o_m = point_to_vector(p(goal%pos)%p)

    CASE (2) !Conglomerate in the cell
       !Recursively calculate mass for each subcell
       DO i = 1, 2
          DO j = 1, 2
             DO k = 1, 2
                IF (ASSOCIATED(goal%subcell(i, j, k)%ptr)) THEN
                   CALL Calculate_masses(goal%subcell(i, j, k)%ptr, p) !Recursively call subcells
                   mass = goal%mass
                   goal%mass = goal%mass + goal%subcell(i, j, k)%ptr%mass !Accumulate mass from subcell
                   goal%c_o_m = (mass * goal%c_o_m + &
                                  goal%subcell(i, j, k)%ptr%mass * goal%subcell(i, j, k)%ptr%c_o_m) / goal%mass !Update center of mass based on contributions from subcells
                END IF
             END DO
          END DO
       END DO
    END SELECT
    
  END SUBROUTINE Calculate_masses

  !Subroutine to calculate forces between particles based on their positions and masses
  SUBROUTINE Calculate_forces(head, p, rji, theta, i_start, i_end, rank)
    TYPE(CELL), POINTER :: head !Head of the tree or root cell
    INTEGER(INT64) :: i !Loop indexing variable
    TYPE(particle3d), INTENT(INOUT) :: p(:) !Particles
    TYPE(vector3d), INTENT(INOUT) :: rji !Vector from one particle to another
    REAL(REAL64) :: theta !Parameter that determines the accuracy of the simulation
    INTEGER, INTENT(IN) :: i_start, i_end, rank !First and last index of the range of particles to process and processor's ID
    
    !Calculate forces for each particle
    DO i = i_start, i_end
       CALL Calculate_forces_aux(i, head, p, rji, theta, rank) !Auxiliary subroutine for force calculation 
    END DO
    
  END SUBROUTINE Calculate_forces

  !Recursive subroutine to calculate forces between a specific particle and the tree structure
  RECURSIVE SUBROUTINE Calculate_forces_aux(goal, tree, p, rji, theta, rank)
    INTEGER(INT64) :: goal, i, j, k !Index of the current particle and loop indexing variables
    TYPE(CELL), POINTER :: tree !Current cell or subcell being processed
    TYPE(particle3d), INTENT(INOUT) :: p(:) !Particles
    TYPE(vector3d), INTENT(INOUT) :: rji !Vector from one particle to another
    REAL(REAL64) :: l, D, theta !Length of the side of the cell, distance between particle and center of mass, and parameter that determines the accuracy of the simulation
    REAL(REAL64) :: r2, r3 !Real variables for squared and cubed distances
    INTEGER :: rank !Processor's ID
    
    SELECT CASE (tree%type)
    CASE (1) !One particle in the cell
       IF (goal .NE. tree%pos) THEN
          rji = tree%c_o_m - point_to_vector(p(goal)%p) !Vector from particle to center of mass
          r2 = (NORM(rji))**2 !Square of the vector
          r3 = r2 * SQRT(r2)  !Cube of the vector
          p(goal)%a = p(goal)%a + p(tree%pos)%m * rji / r3 !Acceleration of the particle
       END IF
       
    CASE (2) !Conglomerate in the cell
       l = tree%range%max%x - tree%range%min%x !Calculate the size of the cell
       rji = tree%c_o_m - point_to_vector(p(goal)%p) !Vector from particle to center of mass
       D = NORM(rji) !Distance between particle and center of mass
       
       IF (l/D .LT. theta) THEN !Barnes-Hut approximation: cell treated as a point mass when is sufficiently distant
          r3 = D**3 !Cube of the distance between particle and center of mass
          p(goal)%a = p(goal)%a + tree%mass * rji / r3 !Accelaration of the particle
       ELSE !Cell too close
          DO i = 1, 2
             DO j = 1, 2
                DO k = 1, 2
                   IF (ASSOCIATED(tree%subcell(i, j, k)%ptr)) THEN
                      CALL Calculate_forces_aux(goal, tree%subcell(i, j, k)%ptr, p, rji, theta, rank) !Recursive call into the subcells to calculate forces on individual particles
                   END IF
                END DO
             END DO
          END DO
       END IF
    END SELECT
    
  END SUBROUTINE Calculate_forces_aux

END MODULE barnes

