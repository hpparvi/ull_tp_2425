MODULE barnes_hut
  USE, INTRINSIC :: iso_fortran_env ! for 64-bit reals
  USE geometry
  USE particle
  IMPLICIT NONE

  ! Contains the edges of the cell
  TYPE range
     TYPE(point3d) :: min,max
  END TYPE range

  
  ! Cell pointer
  TYPE CPtr
     TYPE(CELL), POINTER :: ptr
  END TYPE CPtr

  
  ! Cell. Has its own range, whether it contains particles, the position ??
  ! the total mass, the position of the center of mass, and the octree subcells
  TYPE cell
     TYPE(range)      :: range ! edges of the cell
     TYPE(particle3d) :: part  ! The particle inside the cell, if type=1.
     ! part contains the mass of the particle
     INTEGER          :: pos   ! particle ID
     INTEGER          :: type  ! 0 = no particle; 1 = particle; 2 = aggl.
     TYPE(point3d)    :: c_o_m ! center of mass; used when type=2.
     
     TYPE (CPtr), DIMENSION(2,2,2) :: subcell ! Pointers to each of the 8 subcells
     
  END TYPE cell
 

  ! Initialize the head node
  ALLOCATE(head)

  CALL Calculate_ranges(head) ! space occupied by the head node
  head%type = 0 ! no particle (initialization)
  CALL Nullify_Pointers(head) ! remove all pointers
  

  ! Create the initial tree
  DO i = 1,n ! For all particles
     
     ! Locate the cell where this particle should go
     CALL Find_Cell(head, temp_cell, particles(i)%p)
     
     ! Place the particle inside said cell
     CALL Place_Cell(temp_cell, particles(i)%p, i) ! i is the ID of the particle
                                         ! (called n or pos in the subroutine)
     
  END DO

  
  ! Remove subcells with no particles inside
  CALL Delete_empty_leaves(head)

  ! Calculate the masses (recursively)
  CALL Calculate_masses(head)

  ! Get the initial accelerations (recursive)
  a = 0.0
  CALL Calculate_forces(head)


  ! Main loop
  t_out = 0.0
  DO t = 0.0, t_end, dt
     v = v + a * dt/2
     r = r + v * dt

     !! Las posiciones han cambiado, por lo que tenemos que borrar
     !! y reinicializar el ´arbol
     CALL Delete_tree(head)

     CALL Calculate_ranges(head)
     head%type = 0
     CALL Nullify_Pointers(head)

     ! This updates the tree
     DO i = 1,n
        CALL Find_Cell(head,temp_cell, particles(i)%p)
        CALL Place_Cell(temp_cell, particles(i)%p, i)
     END DO
     
     CALL Delete_empty_leaves(head)
     CALL Calculate_masses(head)

     a = 0.0
     CALL Calculate_forces(head)
     v = v + a * dt/2
     t_out = t_out + dt
     IF (t_out >= dt_out) THEN
        DO i = 1,10
           PRINT*, particles(i)%p  ! change this to output into a file
        END DO
        PRINT*, "-----------------------------------"
        PRINT*, ""
        t_out = 0.0
     END IF
  END DO
  ! End of main loop



CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calculates the ranges of the particles in the
!! matrix r in the 3 dimensions and places it in
!! the variable pointed to by goal
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
SUBROUTINE Calculate_Ranges(goal, particles)
  TYPE(CELL),POINTER :: goal
  TYPE(particle3d)   :: particles
  REAL, DIMENSION(3) :: mins,maxs,medios
  REAL :: span
  
  mins = MINVAL(particles%p, DIM=1)
  maxs = MAXVAL(particles%p, DIM=1)
  
  ! Add 10% so that the particles are not
  ! exactly on the edge
  span = MAXVAL(maxs - mins) * 1.1
  medios = (maxs + mins) / 2.0
  goal%range%min = medios - span/2.0
  goal%range%max = medios + span/2.0
  
END SUBROUTINE Calculate_Ranges




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Finds the cell where the particle will be placed.
!! If the cell we are considering has no particle
!! or one particle, this is the cell where we will
!! place the particle.
!! If the one we are considering is an agglomerate,
!! find (with the BELONGS function) which cell of the
!! 8 possible ones it belongs to, and with this
!! subcell we call Find_Cell again.
!!
!! NOTE: When an agglomerate cell is created, all 8
!! cells are created, so we can assume they always
!! exist. The empty cells are deleted at the very
!! end, when the whole tree has been created.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RECURSIVE SUBROUTINE Find_Cell(root,goal,part)
  TYPE(particle3d) :: part
  TYPE(CELL),POINTER :: root,goal,temp
  INTEGER :: i,j,k

  SELECT CASE (root%type)
  ! More than one particle
  CASE (2)
     ! Check every cell of the octant
     out: DO i = 1,2
        DO j = 1,2
           DO k = 1,2
              ! If you find the subcell the particle belongs to,
              ! then Find_Cell and exit the loop
              IF (Belongs(part, root%subcell(i,j,k)%ptr)) THEN
                 
                 CALL Find_Cell(root%subcell(i,j,k)%ptr, temp, part)
                 goal => temp ! Change the pointer
                 
                 EXIT out ! leave
                 
              END IF
           END DO
        END DO 
     END DO out
     
  ! Case of only one particle   
  CASE DEFAULT
     goal => root ! then, the cell we are in needs no more subdivision
  END SELECT
  
END SUBROUTINE Find_Cell



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Is executed after Find_Cell, in the cell
!! that that function returns, so it is always
!! a type 0 cell (no particle) or type 1 (1
!! particle). If it is type 1 it will be
!! necessary to subdivide the cell and place
!! each particle in the corresponding subcell.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RECURSIVE SUBROUTINE Place_Cell(goal,part,n)
  TYPE(CELL),POINTER :: goal,temp
  TYPE(particle3d) :: part
  INTEGER :: n
  
  SELECT CASE (goal%type)
  ! If there was no particle, great!
  CASE (0)
     goal%type = 1   
     goal%part = part
     goal%pos  = n
     
  ! If there was already a particle, create subcells
  CASE (1)
     CALL Create_Subcells(goal)
     CALL Find_Cell(goal,temp,part)
     CALL Place_Cell(temp,part,n)
     
  ! If there is more than one particle, Find_Cell made a mistake
  CASE DEFAULT
     print*,"SHOULD NOT BE HERE. ERROR!"
     
  END SELECT
  
END SUBROUTINE Place_Cell



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This subroutine is called from Place_Cell
!! and is only called when one particle is in
!! the cell, so there must be a subdivision. 
!! It creates 8 subcells emerging from goal
!! and places the particle in the corresponding
!! subcell.
!!
!! To create the subcells it uses the functions
!! CALCULATE_RANGE, BELONGS and NULLIFY_POINTERS
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Create_Subcells(goal)
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
           
           goal%subcell(i,j,k)%ptr%range%min = Calculate_Range (0,goal,octant)
           goal%subcell(i,j,k)%ptr%range%max = Calculate_Range (1,goal,octant)
           
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
  
END SUBROUTINE Create_Subcells


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Nullifies all pointers in the 8 subcells
!! of the "goal" cell
!!
!! Used in main loop and by CREATE_SUBCELLS
!!
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
!! 
!! Returns TRUE if the considered particle
!! is within the range of the "goal" cell
!!
!! Used by FIND_CELL
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION Belongs (part,goal)
  TYPE(particle3d)    :: part
  TYPE(CELL), POINTER :: goal
  LOGICAL :: Belongs

  IF (part%p%x >= goal%range%min(1) .AND.  &
       part%p%x <= goal%range%max(1) .AND. &
       part%p%y >= goal%range%min(2) .AND. &
       part%p%y <= goal%range%max(2) .AND. &
       part%p%z >= goal%range%min(3) .AND. &
       part%p%z <= goal%range%max(3)) THEN

     Belongs = .TRUE.

  ELSE
     Belongs = .FALSE.
  END IF

END FUNCTION Belongs



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Given an octant, calculate ranges based on
!! the ranges of "goal". If "what" = 0, it
!! calculates the minima. If "what" = 1, it
!! calculates the maxima.
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION Calculate_Range (what,goal,octant)
  INTEGER :: what,n
  TYPE(CELL), POINTER :: goal
  INTEGER, DIMENSION(3) :: octant
  REAL, DIMENSION(3) :: Calculate_Range, valor_medio
  valor_medio = (goal%range%min + goal%range%max) / 2.0
  SELECT CASE (what)
  CASE (0)
     WHERE (octant == 1)
        Calculate_Range = goal%range%min
     ELSEWHERE
        Calculate_Range = valor_medio
     ENDWHERE
  CASE (1)
     WHERE (octant == 1)
        Calculate_Range = valor_medio
     ELSEWHERE
        Calculate_Range = goal%range%max
     ENDWHERE
  END SELECT
  
END FUNCTION Calculate_Range



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Called once the tree is completed to delete
!! (deallocate) the empty cells (no particle). 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RECURSIVE SUBROUTINE Delete_empty_leaves(goal)
  TYPE(CELL),POINTER :: goal
  INTEGER :: i,j,k

  ! If the cell contains subdivision
  IF (ASSOCIATED(goal%subcell(1,1,1)%ptr)) THEN
     ! For every subcell
     DO i = 1,2
        DO j = 1,2
           DO k = 1,2
              ! Continue calling until you reach the smallest subcell
              CALL Delete_empty_leaves(goal%subcell(i,j,k)%ptr)

              ! If no particles inside, deallocate
              IF (goal%subcell(i,j,k)%ptr%type == 0) THEN
                 DEALLOCATE (goal%subcell(i,j,k)%ptr)
                 
              END IF
           END DO
        END DO
     END DO
  END IF
  
END SUBROUTINE Delete_empty_leaves


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Deletes the entire tree save for "head"
!!
!! The tree must be re-generated, so we need to
!! delete the previous one to avoid memory leaks.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RECURSIVE SUBROUTINE Delete_tree(goal)
  TYPE(CELL),POINTER :: goal
  INTEGER :: i,j,k
  DO i = 1,2
     DO j = 1,2
        DO k = 1,2
           IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
              CALL Delete_tree(goal%subcell(i,j,k)%ptr)
              DEALLOCATE (goal%subcell(i,j,k)%ptr)
           END IF
        END DO
     END DO
  END DO
  
END SUBROUTINE Delete_tree



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calculates the mass and center-of-mass for
!! all child subcells of "goal".
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RECURSIVE SUBROUTINE Calculate_masses(goal, particles)
  TYPE(CELL),POINTER :: goal
  TYPE(particle3d)   :: particles
  INTEGER :: i,j,k
  REAL :: mass
  TYPE(point3d) :: c_o_m
  
  goal%mass = 0
  goal%c_o_m = 0
  
  SELECT CASE (goal%type)
  ! If there is one particle
  CASE (1)
     goal%mass = particles(goal%pos)%m
     goal%c_o_m = particles(goal%pos)%p

  ! If there is an agglomerate
  CASE (2)
     ! For all subcells
     DO i = 1,2
        DO j = 1,2
           DO k = 1,2
              ! If the cell has subcells, use recursivity
              IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                 
                 CALL Calculate_masses(goal%subcell(i,j,k)%ptr)
                 
                 mass = goal%mass

                 goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass

                 ! c_o_m = m*vec(r); add the considered particle
                 ! and the other subcells at the same level
                 goal%c_o_m = (mass * goal%c_o_m + &
                      goal%subcell(i,j,k)%ptr%mass * &
                      goal%subcell(i,j,k)%ptr%c_o_m) / goal%mass
                 
              END IF
           END DO
        END DO
     END DO
  END SELECT
  
END SUBROUTINE Calculate_masses


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calculate the forces of all the particles
!! against "head".
!!
!! Uses the function "Calculate_forces_aux"
!! which actually does the calculations for
!! each particle
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Calculate_forces(head)
  TYPE(CELL),POINTER :: head
  INTEGER :: i,n
  
  DO i = 1,n
     CALL Calculate_forces_aux(i,head)
  END DO
  
END SUBROUTINE Calculate_forces



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Given a "goal" particle, calculate the forces
!! on it from the "tree" cell.
!! 
!! If "tree" is an agglomerate, see if l/D < theta.
!! If it is, treat the cell as a single particle;
!! otherwise, use all the subcells to calculate the
!! forces.
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RECURSIVE SUBROUTINE Calculate_forces_aux(goal, tree, particles)
  TYPE(CELL),POINTER :: tree
  INTEGER :: i,j,k,goal
  REAL(real64) :: l,D
  TYPE(vector3d) :: rji
  REAL(real64) :: r2, r3

  SELECT CASE (tree%type)
  ! One particle
  CASE (1)
     IF (goal .NE. tree%pos) THEN
        rji = tree%c_o_m - particles(goal)%px
        r2 = SUM(rji**2)
        r3 = r2 * SQRT(r2)
        a(goal) = a(goal) + particles(tree%pos)%m * rji / r3
     END IF

  ! An agglomerate; check how to treat it
  CASE (2)
     !! The range is the same (span) in all 3 dimensions.
     !! Any dimension works to find the side length.
     l = tree%range%max(1) - tree%range%min(1)
     rji = tree%c_o_m - particles(goal)%p
     r2 = SUM(rji**2)
     D = SQRT(r2)

     ! Case l/D < theta: calculate as if it was one particle
     IF (l/D < theta) THEN
        r3 = r2 * D
        a(goal) = a(goal) + tree%mass * rji / r3
        
     ! Case l/D > theta: Calculate the forces using all subcells
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


END MODULE barnes_hut
