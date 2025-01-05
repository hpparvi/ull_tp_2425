module tree_algorithm
  use mpi_f08
  use iso_fortran_env
  use geometry
  use particle
  implicit none
  
  !! Define types used in Barnes-Hut algorithm
  !! to create the tree
  
  ! type range used to save ...
  TYPE RANGE
  REAL(real64), DIMENSION(3) :: min_val,max_val 
  END TYPE RANGE
  
  ! Type cell as a pointer
  TYPE CPtr
    TYPE(CELL), POINTER :: ptr !pointer of the 'target' cell
  END TYPE CPtr
 
 ! Cell include the info of the particle/s in the cell
  TYPE CELL
    TYPE (RANGE) :: range
    TYPE(particle3d) :: part!(particle info in the cell)
    INTEGER :: pos  !! id of the particle
    INTEGER :: type !! 0 = no particle; 1 = particle; 2 = conglomerate
    TYPE(point3d) :: c_o_m ! center of mass of that cell
    TYPE (CPtr), DIMENSION(2,2,2) :: subcell ! subcells in main cell
    ! that contains cells as well 
  END TYPE CELL

  CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_Ranges
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calculate particles ranges using the position
!! contained in particles in three dimensions 
!! and they are saved in the variable aim by goal
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calculate_Ranges(goal, particles)
    TYPE(CELL), POINTER :: goal
    TYPE(particle3d), DIMENSION(:), INTENT(in) :: particles	
    REAL(real64), DIMENSION(3):: mins,maxs,medios 
    REAL (real64) :: span
    ! Find the minimum and maximum values
    ! of the particle positions
    mins(1) = MINVAL(particles%p%x, DIM=1) 
    maxs(1) = MAXVAL(particles%p%x, DIM=1)
    mins(2) = MINVAL(particles%p%y, DIM=1) 
    maxs(2) = MAXVAL(particles%p%y, DIM=1)
    mins(3) = MINVAL(particles%p%z, DIM=1) 
    maxs(3) = MAXVAL(particles%p%z, DIM=1)
  
    ! To calculate the span, it is added a 10% 
    ! to the subtraction of max and min values
    ! that give us the range of the cell, so any
    ! of the particles will be in the edge 
    span = MAXVAL(maxs - mins) * 1.1
    ! Calculate mean of min and max values
    ! in each dimension to use it as center 
    ! of the cell
    medios = (maxs + mins)/2. 
    goal%range%min_val = medios - span/2.  
    goal%range%max_val = medios + span/2.
    
  END SUBROUTINE Calculate_Ranges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Find_Cell
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Find the cell where we place the particle.
!! If the cell that we are considering (root) 
!! hasn't got a particle or just one of them,
!!this is the cell where we place the particle 
!! (goal).

!! If the cell we are considering in a conglomerate, 
!! we look for a which subcell of the 8 of them possible
!! the particle belongs, with the function BELONGS and
!! with this subcell we call Find_Cell again.
!!
!! NOTE: When a cell is created as a type 2 (conglomerate) 
!! it's also created the 8 subcells, thus we can assume that
!! always siempre existen the 8 subcells. The empty cell are
!! removed at the end, when the tree has been creates.

!! This is a recursive subroutine, so it calls itself in order
!! to make easier (and possible) to find the cell in case of a
!! conglomerate of particles, when each cell is divided and 
!! look in every subcell to find where the particle should
!! be placed. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE Find_Cell(root,goal,part)
    type(particle3d), INTENT(in) :: part
    TYPE(CELL),POINTER, INTENT(inout) :: root,goal
    TYPE(CELL), POINTER :: temp
    INTEGER :: i,j,k
    SELECT CASE (root%type)
      ! In case of a conglomerate, then use belongs to 
      ! find the subcell of the particle. 
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
      ! In case that only has one particle, then particle is
      ! in the current cell.
      CASE DEFAULT
        goal => root
    END SELECT
  END SUBROUTINE Find_Cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Belongs
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Return TRUE if the particle "part" is into
!! the range of the cell "goal". 
!!
!! Used by FIND_CELL
!!NOTE: In order to avoid problems, it is imposed that
!! if a particle is in the edge of the cell, then
!! it belongs to the cell where the particle is
!! in the minimum value, instead of the maximum. 
!! It could be the same whether is placed in the maximum
!! or the minimum, but if both limits are included,
!! then the particle belongs to two neighbouring cells,
!! and that's a problem.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION Belongs (part,goal)
    TYPE(particle3d) :: part
    TYPE(CELL), POINTER :: goal
    LOGICAL :: Belongs
    
    IF ((part%p%x >= goal%range%min_val(1)) .AND. &
      (part%p%x < goal%range%max_val(1)) .AND. &
      (part%p%y >= goal%range%min_val(2)) .AND. &
      (part%p%y < goal%range%max_val(2)) .AND. &
      (part%p%z >= goal%range%min_val(3)) .AND. &
      (part%p%z < goal%range%max_val(3))) THEN
        Belongs = .TRUE.
    
    ELSE
      Belongs = .FALSE.
    END IF
  END FUNCTION Belongs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Place_Cell
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! It's executed after Find_Cell, in the cell that
!! the function returned, so it's always a type 0 (without 
!! particle) or type 1 (with one particle). In the case of
!! type 1 cell, the cell must be sudivided and place the 
!! both particles (the one that was already in, and the new
!! one) in their correct subcells.
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE Place_Cell(goal,part,n)
    TYPE(CELL),POINTER, INTENT(inout) :: goal
    TYPE(CELL),POINTER :: temp
    TYPE(particle3d), INTENT(in) :: part
    INTEGER, INTENT(in) :: n
    SELECT CASE (goal%type)
      CASE (0)
        goal%type = 1
        goal%part = part
        goal%pos = n
    
      CASE (1)
        CALL Crear_Subcells(goal)
        CALL Find_Cell(goal,temp,part)
        CALL Place_Cell(temp,part,n)
      ! In case of a type 2 cell, there's still a
      ! conglomerate of particles, so must be a
      ! problem with the previous functions used.
      CASE DEFAULT
        print*,"SHOULD NOT BE HERE. ERROR!"
    END SELECT
  END SUBROUTINE Place_Cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Crear_Subcells
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This function is only call when there is a
!! particle already in the cell, so we need 
!! to subdivide the cell (goal) in 8 subcells
!! and place the particle in goal in their 
!! corresponding subcell. 
!!
!! To create the subcells, it is used:
!! CALCULAR_RANGE, BELONGS and NULLIFY_POINTERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Crear_Subcells(goal)
    TYPE(CELL), POINTER, INTENT(inout) :: goal
    TYPE(particle3d) :: part
    INTEGER :: i,j,k !,n
    INTEGER, DIMENSION(3) :: octant ! aim which octant is
    part = goal%part
    goal%type = 2
    DO i = 1,2
      DO j = 1,2
        DO k = 1,2
          octant = (/i,j,k/) 
          ! creation of the subcells and their range
          ALLOCATE(goal%subcell(i,j,k)%ptr)
          goal%subcell(i,j,k)%ptr%range%min_val = Calcular_Range(0,goal,octant)
          goal%subcell(i,j,k)%ptr%range%max_val = Calcular_Range(1,goal,octant)
          ! assign the particle if it's in the subcell
          IF (Belongs(part,goal%subcell(i,j,k)%ptr)) THEN
            goal%subcell(i,j,k)%ptr%part = part
            goal%subcell(i,j,k)%ptr%type = 1
            goal%subcell(i,j,k)%ptr%pos = goal%pos
          ! if not, keep it empty
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
!! Just NULLIFY pointers of the 8 subcells
!! of the "goal" cell
!!
!! It's used in the main loop of the program and
!! by CREAR_SUBCELLS.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Nullify_Pointers(goal)
    TYPE(CELL), POINTER, INTENT(inout) :: goal
    INTEGER :: i,j,k
    ! check every 'leaves'/octant pointer is null (no direction)
    DO i = 1,2
      DO j = 1,2
        DO k = 1,2
          NULLIFY(goal%subcell(i,j,k)%ptr)
        END DO
      END DO
    END DO
  END SUBROUTINE Nullify_Pointers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calcular_Range
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Given an octant like (1,1,1, 1,1,2 ... 2,2,2),
!! calculate its ranges in base of the ranges of 
!! "goal" cell. If "what" = 0 calculate minimums.
!! If "what" = 1 calculate the maximum.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION Calcular_Range (what,goal,octant)
    INTEGER :: what
    TYPE(CELL), POINTER :: goal
    INTEGER, DIMENSION(3) :: octant
    REAL(real64), DIMENSION(3) :: Calcular_Range, valor_medio
    ! mean value of the ranges to establish where the octant limit
    ! with another octant
    valor_medio = (goal%range%min_val + goal%range%max_val) / 2.0

    SELECT CASE (what)
      CASE (0)
        WHERE (octant == 1)
          Calcular_Range = goal%range%min_val
        ELSEWHERE
          Calcular_Range = valor_medio
        ENDWHERE
    
      CASE (1)
        WHERE (octant == 1)
          Calcular_Range = valor_medio
        ELSEWHERE
          Calcular_Range = goal%range%max_val
        ENDWHERE
    END SELECT
  END FUNCTION Calcular_Range

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Borrar_empty_leaves
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! When the tree has been completed, this 
!! subroutine remove the empty leaves of the 
!! cell(i.e. the octants without particle).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE Borrar_empty_leaves(goal)
    TYPE(CELL),POINTER, INTENT(inout) :: goal
    INTEGER :: i,j,k
    ! choose any of the octant of the subcell
    ! that at first is pointed at a particle 
    IF (ASSOCIATED(goal%subcell(1,1,1)%ptr)) THEN
      DO i = 1,2
        DO j = 1,2
          DO k = 1,2
            ! the subcells without particle are removed
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
!! Remove the complete tree, except "head".
!!
!! The tree change every time so we have to 
!! delete the last one to avoid "memory leaks".
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE Borrar_tree(goal)
    TYPE(CELL), POINTER, INTENT(inout) :: goal
    INTEGER :: i,j,k
    DO i = 1,2
      DO j = 1,2
        DO k = 1,2
          ! if the cell is associated to subcells
          ! then now the whole tree is removed
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
!! Calculate for all the cells contained in
!! "goal" their mass and their center-of-mass.
!!
!! As this subroutine is also recursive, 
!! it is calculated from the largest cell, and
!! subdivided in each iteration, looking for
!! the subcells that contains only one particle.
!! The larger cell (contains a conglomerate),
!! has a mass that corresponds to the addition
!! of the masses of every particle in it, and
!! the c.o.m of the particles distribution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE Calculate_masses(goal, particles)
    TYPE(CELL),POINTER, INTENT(inout) :: goal
    TYPE(particle3d), DIMENSION(:), INTENT(in) :: particles
    INTEGER :: i,j,k
    REAL(real64) :: mass
    type(point3d) :: c_o_m
    ! First both of them have null values
    goal%part%m = 0.
    goal%c_o_m = point3d(0.,0.,0.)
    
    SELECT CASE (goal%type)
      
      CASE (1) ! if only one particle, cell mass and c.o.m. are
      ! the particle values of mass and position 
      goal%part%m = particles(goal%pos)%m ! mass of the particle in the cell 
      goal%c_o_m= particles(goal%pos)%p ! position of the particle
      
      CASE (2) ! in case of a conglomerate, look in every subcell to see if 
      ! it has one (or various) particle(s).
        DO i = 1,2
          DO j = 1,2
            DO k = 1,2
              IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                CALL Calculate_masses(goal%subcell(i,j,k)%ptr, particles)
                 mass = goal%part%m 
                 goal%part%m = goal%part%m + goal%subcell(i,j,k)%ptr%part%m
	         
	         goal%c_o_m%x = (mass * goal%c_o_m%x + &
	         goal%subcell(i,j,k)%ptr%part%m * goal%subcell(i,j,k)%ptr%c_o_m%x) / goal%part%m
	         goal%c_o_m%y = (mass * goal%c_o_m%y + &
	         goal%subcell(i,j,k)%ptr%part%m * goal%subcell(i,j,k)%ptr%c_o_m%y) / goal%part%m
	         goal%c_o_m%z = (mass * goal%c_o_m%z + &
	         goal%subcell(i,j,k)%ptr%part%m * goal%subcell(i,j,k)%ptr%c_o_m%z) / goal%part%m
	      END IF
	    END DO
          END DO
        END DO
    END SELECT
    !print *, 'The mass of the cell is', goal%part%m
    !print *, 'The com of the cell is ', goal%c_o_m
  END SUBROUTINE Calculate_masses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_forces
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Calculate forces of all of the particles 
!! exerted on "head". 
!! Use the function Calculate_forces_aux, 
!! which is the one that actually calculate
!! the force for every particle.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calculate_forces(head,particles,n,acc)
    TYPE(CELL), POINTER, INTENT(in) :: head
    TYPE(particle3d), DIMENSION(:), INTENT(in) :: particles
    TYPE(vector3d), DIMENSION(:), INTENT(inout) :: acc
    INTEGER:: n 
    INTEGER :: i
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!! Parallel programming with OpenMP !!!!      
!!!!! This subroutine parallelise do loop and
!!!!! that calculate (for each particle) the 
!!!!! interaction between particles. 
!!!!! WARNING: omp parallel must be opened in the main code
!!!!! in order this works (it could be a better optimization
!!!!! of the code if a conditional open the  omp parallel if
!!!!! it's not previously opened(?).)      
    !print *, "Before parallel"
    !!$omp parallel private(nt, tid, i) shared(head,particles,acc)
    
    !!$ nt = omp_get_num_threads()
    !!$ tid = omp_get_thread_num()
    !!$omp do
    DO i = 1,n 
      CALL Calculate_forces_aux(i,head,particles,acc)
      !print '(xA,i3,xA,i2,xA,i3)', "Do loop, i =", i, &
          !& "thread ", tid, " of ", nt
    END DO
    !!$omp end do
    !!$omp end parallel
    !print *, "After parallel"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
  END SUBROUTINE Calculate_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate_forces_aux
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Given a particle "goal" (which is the index
!!of the target particle), calculate the forces 
!! exerted on it in the cell "tree". 
!! If "tree" is a cell that contains only one particle, 
!! then the attraction is between two particles, and
!! that's simple :)
!!
!! But if "tree" is a conglomerate cell, first we see if
!! l/r < theta, i.e. if the range of the cell (l)
!! divided by the distance of the particle "goal" to
!! the center_of_mass of the cell "tree" (r) is less than theta.
!! In case it's true, we treat the cell like a single 
!! particle. In case that it's not less than theta,
!! then we have to consider all the subcells of the 
!! tree and for each of them call Calculate_forces_aux in a
!! recursive way.
!! Note: The original condition we set (theta = 1.) is that  
!! the distance between particles must be greater than the 
!! range of the cell, which in some cases may be insufficient 
!! accuracy, as was noted by some classmates, and it is necessary 
!! for the calculation to reduce the value of theta.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE Calculate_forces_aux(goal,tree,particles,acc)
    TYPE(CELL),POINTER, INTENT(in) :: tree
    TYPE(particle3d), DIMENSION(:), INTENT(in)  :: particles
    TYPE(vector3d), DIMENSION(:), INTENT(inout) :: acc
    INTEGER :: i,j,k,goal
    REAL(real64) :: l,r,r3, theta = 0.01
    TYPE(vector3d) :: rji
    
    SELECT CASE (tree%type)
      CASE (1) ! only one particle in the tree
        IF (goal .NE. tree%pos) THEN
          rji = tree%c_o_m - particles(goal)%p 
          r = distance(tree%c_o_m, particles(goal)%p)
          acc(goal) = acc(goal) + particles(tree%pos)%m * rji / (r**3)
        END IF
    
      CASE (2)
      !! The range has the same span in the three dimensions 
      !! so we can consider any dimension to calculate the 
      !! side of the cell (in this case the x-dimension)
        l = tree%range%max_val(1) - tree%range%min_val(1) 
        rji = tree%c_o_m - particles(goal)%p
        r = distance(tree%c_o_m, particles(goal)%p)
        IF (l/r < theta) THEN
        r3 = r**3
        acc(goal) = acc(goal) + tree%part%m * rji / r3
      
        ELSE
          DO i = 1,2
            DO j = 1,2
	      DO k = 1,2
	        IF (ASSOCIATED(tree%subcell(i,j,k)%ptr)) THEN
	          CALL Calculate_forces_aux(goal,tree%subcell(i,j,k)%ptr,particles,acc)
	        END IF
	      END DO
	    END DO
          END DO
        END IF
    END SELECT
    ! note: this code could have issues in the integration if the particles are close enough (r aprox 0),
    ! so it could be necessary to add a softening length in some cases (as in ex1)
  END SUBROUTINE Calculate_forces_aux

end module tree_algorithm
