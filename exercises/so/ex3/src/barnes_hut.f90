MODULE barnes_hut
  use geometry
  use particles

  use mpi_f08

  !$use_omp_lib
  
  IMPLICIT NONE

  ! To set the range of each cell in the tree
  TYPE RANGE
     REAL, DIMENSION(3) :: min,max
  END TYPE RANGE

  ! Cell pointer
  TYPE CPtr
     TYPE(CELL), POINTER :: ptr
  END TYPE CPtr

  ! Cell (or node) of the tree
  TYPE CELL
     TYPE (RANGE) :: range
     TYPE(particle) :: part                     ! particles
     INTEGER :: pos                             ! index of the particle contained in it
     INTEGER :: type                            ! 0 = no particle; 1 = particle; 2 = group of particles
     REAL :: mass                               ! total mass of the cell
     type(point3d) :: c_o_m                     ! position of the center of mass
     TYPE (CPtr), DIMENSION(2,2,2) :: subcell   ! pointer to the subcells of the current cell
  END TYPE CELL

CONTAINS
!-------------------------------------------------------------
  ! Calculations of the ranges of the particles in the
  ! array of particles p_arr
  ! in 3d and saves it in goal variable
  
  SUBROUTINE Calculate_Ranges(goal, p_arr)
    TYPE(CELL), POINTER :: goal
    TYPE(PARTICLE), DIMENSION(:), INTENT(IN) :: p_arr
    REAL, DIMENSION(3) :: mins, maxs, medios
    REAL :: span

    mins = [MINVAL([p_arr(:)%p%x]), MINVAL([p_arr(:)%p%y]), MINVAL([p_arr(:)%p%z])]
    maxs = [MAXVAL([p_arr(:)%p%x]), MAXVAL([p_arr(:)%p%y]), MAXVAL([p_arr(:)%p%z])]
    ! Adding a 10% to the span to avoid particles on borders
    span = MAXVAL(maxs - mins) * 1.1
    medios = (maxs + mins) / 2.0
    goal%range%min = medios - span/2.0
    goal%range%max = medios + span/2.0
  END SUBROUTINE Calculate_Ranges

!-------------------------------------------------------------
  ! Nullifyies the pointers of a given cell (goal)
  ! Used in the main loop (Block 2) and in the
  ! subroutine CREAR_SUBCELLS below

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

!------------------------------------------------------------- 
  !! Finds the cell where particle p_i is to be placed.
  ! If the cell is empty or has just 1 particle, it places it there
  ! (type =0 or type =1).
  ! If it is already populated by a group of particles (type=2),
  ! we locate the specific subcell were it BELONGS and
  ! apply recursively the find_cell subrourtine.
  !
  ! Note: when a cell has a group of particles (type=2), the
  ! 8 subcells are created so we can assume they always exist.
  ! The unused empty cells are erased when the whole tree is created.

  RECURSIVE SUBROUTINE Find_Cell(root, goal, p_i)
    type(particle), intent(in) :: p_i
    TYPE(CELL), POINTER :: root, goal, temp
    INTEGER :: i, j, k
    ! Depending on the situation
    SELECT CASE (root%type)
    CASE (2) ! if root has a group of particles we run through its subcells
       out: DO i = 1,2
          DO j = 1,2
             DO k = 1,2
                IF (Belongs(p_i,root%subcell(i,j,k)%ptr)) THEN
                   ! if it is inside, call to the routine again with the
                   ! subcell as the new root
                   CALL Find_Cell(root%subcell(i,j,k)%ptr,temp,p_i)
                   goal => temp ! save the value and exit loop
                   EXIT out
                END IF
             END DO
          END DO
       END DO out 
    CASE DEFAULT
       goal => root
    END SELECT
  END SUBROUTINE Find_Cell

!------------------------------------------------------------- 
  ! Belongs is TRUE if the particle "part" is inside the range
  ! of the cell "goal", and FALSE if it is not.
  !
  ! Used by FIND_CELL
  
  FUNCTION Belongs (p_i,goal)
    type(particle), intent(in) :: p_i
    TYPE(CELL), POINTER :: goal
    LOGICAL :: Belongs
    IF (p_i%p%x > goal%range%min(1) .AND. &
         p_i%p%x <= goal%range%max(1) .AND. &
         p_i%p%y > goal%range%min(2) .AND. &
         p_i%p%y <= goal%range%max(2) .AND. &
         p_i%p%z > goal%range%min(3) .AND. &
         p_i%p%z <= goal%range%max(3)) THEN
       Belongs = .TRUE.
    ELSE
       Belongs = .FALSE.
    END IF
  END FUNCTION Belongs

!-------------------------------------------------------------  
  ! Place_Cell is run after Find_cell, to store a particle in
  ! a given cell (always empty or one particle, type=0,1).
  ! If type=1, it divides the cell and places each particle in
  ! the corresponding subcell.

  RECURSIVE SUBROUTINE Place_Cell(goal, p_i, n)
    TYPE(CELL),POINTER :: goal, temp
    type(particle), intent(in) :: p_i
    INTEGER :: n
    
    SELECT CASE (goal%type)
    CASE (0) !if the cell is empty, directly assign it
       goal%type = 1
       goal%part = p_i
       goal%pos = n
    CASE (1) !if it has one particle, create subcells, find the corresponding
       CALL Create_Subcells(goal)
       CALL Find_Cell(goal,temp,p_i)
       CALL Place_Cell(temp,p_i,n)
    CASE DEFAULT
       print*,"SHOULD NOT BE HERE. ERROR!"
    END SELECT
  END SUBROUTINE Place_Cell

!-------------------------------------------------------------  
  ! Create_Subcells is called in Place_Cell when there is already
  ! a particle in the goal cell and subcells are required.
  ! Creates 8 subcells below goal and places its particle
  ! in the corresponding subcell

  SUBROUTINE Create_Subcells(goal)
    TYPE(CELL), POINTER :: goal
    type(particle) :: part
    INTEGER :: i, j, k
    INTEGER, DIMENSION(3) :: octant ! id of the octant
    part = goal%part
    goal%type=2
    DO i = 1,2
       DO j = 1,2
          DO k = 1,2
             octant = (/i,j,k/)
             ALLOCATE(goal%subcell(i,j,k)%ptr)
             goal%subcell(i,j,k)%ptr%range%min = Calculate_Range_minmax(0,goal,octant)
             goal%subcell(i,j,k)%ptr%range%max = Calculate_Range_minmax(1,goal,octant)
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

!-------------------------------------------------------------  
  ! Calculate_Range_minmax calculates the ranges of a given
  ! octant id (i,j,k), respective to the ranges of goal.
  ! - what = 0 calculates mins
  ! - what = 1 calculates maxs
  ! Used in Create_Subcells

  FUNCTION Calculate_Range_minmax (what,goal,octant)
    INTEGER :: what
    TYPE(CELL), POINTER :: goal
    INTEGER, DIMENSION(3) :: octant
    REAL, DIMENSION(3) :: Calculate_Range_minmax, mean_value
    mean_value = (goal%range%min + goal%range%max) / 2.0
    SELECT CASE (what)
    CASE (0)
       WHERE (octant == 1)
          Calculate_Range_minmax = goal%range%min
       ELSEWHERE
          Calculate_Range_minmax = mean_value
       ENDWHERE
    CASE (1)
       WHERE (octant == 1)
          Calculate_Range_minmax = mean_value
       ELSEWHERE
          Calculate_Range_minmax = goal%range%max
       ENDWHERE
    END SELECT
  END FUNCTION Calculate_Range_minmax

!-------------------------------------------------------------  
  ! Delete_empty_leaves is called when the tree is completed
  ! and it deallocates all empty cells

  RECURSIVE SUBROUTINE Delete_empty_leaves(goal)
    TYPE(CELL),POINTER :: goal
    INTEGER :: i, j, k
    IF (ASSOCIATED(goal%subcell(1,1,1)%ptr)) THEN
       DO i = 1,2
          DO j = 1,2
             DO k = 1,2
                CALL Delete_empty_leaves(goal%subcell(i,j,k)%ptr)
                IF (goal%subcell(i,j,k)%ptr%type == 0) THEN
                   DEALLOCATE (goal%subcell(i,j,k)%ptr)
                END IF
             END DO
          END DO
       END DO
    END IF
  END SUBROUTINE Delete_empty_leaves

!-------------------------------------------------------------  
  ! Calculate_masses calculates the mass and
  ! center-of-mass of all cells hanging from goal

  RECURSIVE SUBROUTINE Calculate_masses(goal, p_arr)
    TYPE(CELL),POINTER :: goal
    type(particle), dimension(:), intent(in) :: p_arr
    INTEGER :: i,j,k
    REAL :: mass

    goal%mass = 0
    goal%c_o_m = point3d(0,0,0)
    SELECT CASE (goal%type)
    CASE (1)
       goal%mass = p_arr(goal%pos)%m
       goal%c_o_m = p_arr(goal%pos)%p
    CASE (2)
       DO i = 1,2
          DO j = 1,2
             DO k = 1,2
                IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                   CALL Calculate_masses(goal%subcell(i,j,k)%ptr, p_arr)
                   mass = goal%mass
                   goal%mass = goal%mass + goal%subcell(i,j,k)%ptr%mass
                   goal%c_o_m = (mass * goal%c_o_m + &
                        goal%subcell(i,j,k)%ptr%mass * goal%subcell(i,j,k)%ptr%c_o_m) / goal%mass
                END IF
             END DO
          END DO
       END DO
    END SELECT
  END SUBROUTINE Calculate_masses

!-------------------------------------------------------------  
  ! Calculate_forces calculates the gravitatory force felt by
  ! the cell head. Real calculations are performed by
  ! the recursive subroutine Calculate_forces_aux()

  SUBROUTINE Calculate_forces(head, p_arr, n, theta, epsilon, n_send, master_rank, rank)
    TYPE(CELL),POINTER :: head
    type(particle), dimension(:) :: p_arr
    INTEGER :: i, n
    real, intent(in) :: theta, epsilon
    integer :: i_send_min, i_send_max
    integer, intent(in) :: master_rank, rank
    integer, dimension(:), intent(in) :: n_send

    ! first particle by each node: rank (from 0 to size) times the amount of
    if (rank.eq.master_rank) then
       i_send_min = 1
    else
       i_send_min = SUM(n_send(:rank)) +1
    end if
    i_send_max = SUM(n_send(:rank+1))
    DO i = i_send_min, i_send_max
       CALL Calculate_forces_aux(i,head, p_arr, theta, epsilon)
    END DO
    
  END SUBROUTINE Calculate_forces

!-------------------------------------------------------------  
  ! Calculate_forces_aux calculates the gravitational forces
  ! of the whole tree under the cell 'tree' on the particle 'goal'.
  !
  ! If there is only one particle in the cell 'tree', it is an easy case.
  !
  ! If 'tree' contains a group of particles:
  ! - if the side of the cell (l) compared to the distance (D) between its
  !   center of mass and the goal particle is small (l/D<theta),
  !   it will be considered as a single particle itself.
  ! - if the cell is close to the goal particle (l/D>theta), all subcells
  !   under the cell have to be considered, calling the function recursively.

  RECURSIVE SUBROUTINE Calculate_forces_aux(goal, tree, p_arr, theta, epsilon)
    TYPE(CELL), POINTER :: tree
    INTEGER :: i, j, k, goal
    type(particle), dimension(:) :: p_arr
    REAL :: l, D, r
    real, intent(in) :: theta, epsilon
    type(vector3d) :: rji_v
    
    SELECT CASE (tree%type)
    CASE (1)
       IF (goal .NE. tree%pos) THEN ! not considering the goal particle for calculations!
          rji_v = vecpp(p_arr(goal)%p, tree%c_o_m)
          r = norm(rji_v)
          p_arr(goal)%a = p_arr(goal)%a + p_arr(tree%pos)%m * rji_v /((r**2 + epsilon**2) * r)
       END IF
    CASE (2)
       ! The span of the range is the same in every dimension, so dim1 is chosen for example
       l = tree%range%max(1) - tree%range%min(1)
       rji_v = vecpp(p_arr(goal)%p, tree%c_o_m)
       D = norm(rji_v)
       ! if far enough to approximate by center of mass
       IF (l/D < theta) THEN
          !r3 = D**3
          p_arr(goal)%a = p_arr(goal)%a + tree%mass * rji_v / ((D**2 + epsilon**2) * D)
       ! if not, go through all subcells recursively
       ELSE
          DO i = 1,2
             DO j = 1,2
                DO k = 1,2
                   ! empty cells have been deleted so we can avoid them
                   IF (ASSOCIATED(tree%subcell(i,j,k)%ptr)) THEN
                      CALL Calculate_forces_aux(goal,tree%subcell(i,j,k)%ptr, p_arr, theta, epsilon)
                   END IF
                END DO
             END DO
          END DO
       END IF
    END SELECT
  END SUBROUTINE Calculate_forces_aux

!-------------------------------------------------------------  
  ! Deletes the whole tree except the 'head' cell.
  ! In this way the tree can be re-generated from
  ! scratch in each iteration.

  RECURSIVE SUBROUTINE Delete_Tree(goal)
    TYPE(CELL),POINTER :: goal
    INTEGER :: i,j,k
    DO i = 1,2
       DO j = 1,2
          DO k = 1,2
             IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
                CALL Delete_Tree(goal%subcell(i,j,k)%ptr)
                DEALLOCATE (goal%subcell(i,j,k)%ptr)
             END IF
          END DO
       END DO
    END DO
  END SUBROUTINE Delete_Tree

!-------------------------------------------------------------  
  ! Make_Tree calls all previous functions to make the tree from scratch
  
  SUBROUTINE Make_Tree(root, p_arr, n, theta, epsilon, t_calcs, n_send, master_rank, rank)
    type(CELL), pointer :: root, temp_cell
    type(particle), dimension(:) :: p_arr
    real, dimension(7), intent(inout) :: t_calcs
    integer :: t_dummy, t_dummy_2   !, t_rate
    INTEGER :: n, i
    real, intent(in) :: theta, epsilon
    integer, dimension(:), intent(in) :: n_send
    integer, intent(in) :: master_rank, rank
    
    ! calculates its min and max values
    !call system_clock(t_dummy, t_rate)
    if (rank.eq.master_rank) then
       t_dummy = MPI_Wtime()
    end if
    
    CALL Calculate_Ranges(root, p_arr)
    !call system_clock(t_dummy_2)
    if (rank.eq.master_rank) then
       t_dummy_2 = MPI_Wtime()
       t_calcs(1) = t_calcs(1) + (t_dummy_2 - t_dummy)
    end if
    
    root%type = 0 ! empty cell
    ! remove the subcell pointers
    CALL Nullify_Pointers(root)
    !call system_clock(t_dummy)
    if (rank.eq.master_rank) then
       t_dummy = MPI_Wtime()
       t_calcs(2) = t_calcs(2) + (t_dummy- t_dummy_2)
    end if
    
    ! Creation of the tree, particle by particle
    DO i = 1,n
       ! find the corresponding cell and place particle there
       CALL Find_Cell(root,temp_cell,p_arr(i))
       CALL Place_Cell(temp_cell,p_arr(i),i)
    END DO
    ! call system_clock(t_dummy_2)
    if (rank.eq.master_rank) then
       t_dummy_2 = MPI_Wtime()
       t_calcs(3) = t_calcs(3) + (t_dummy_2-t_dummy)
    end if
    
    ! Clean the tree 
    CALL Delete_empty_leaves(root)
    !call system_clock(t_dummy)
    if (rank.eq.master_rank) then
       t_dummy = MPI_Wtime()
       t_calcs(4) = t_calcs(4) + (t_dummy-t_dummy_2)
    end if

    ! Calculate masses in each cell
    CALL Calculate_masses(root, p_arr)
    !call system_clock(t_dummy_2)
    if (rank.eq.master_rank) then
       t_dummy_2 = MPI_Wtime()
       t_calcs(5) = t_calcs(5) + (t_dummy_2-t_dummy)
    end if
    
    ! Initialization of accelerations
    do i = 1,n
     p_arr(i)%a = vector3d(0,0,0)
    end do
    !call system_clock(t_dummy)
    if (rank.eq.master_rank) then
       t_dummy = MPI_Wtime()
       t_calcs(6) = t_calcs(6) + (t_dummy-t_dummy_2)
    end if

    ! Calculation of forces
    ! distribute particle information among ranks
    !!! ME QUEDO POR AQUÍ. FALTA MANDAR LA INFO DE QUÉ PARTÍCULAS USAR A CADA NODO...
    
    CALL Calculate_forces(root, p_arr, n, theta, epsilon, n_send, master_rank, rank)
    !call system_clock(t_dummy_2)
    if (rank.eq.master_rank) then
       t_dummy_2 = MPI_Wtime()
       t_calcs(7) = t_calcs(7) + (t_dummy_2-t_dummy)
    end if
    ! wait for all processes
    CALL MPI_Barrier(MPI_COMM_WORLD)

    ! all processes should gather the updated particle information in each other process
    !CALL MPI_Allgatherv(MPI_IN_PLACE, n_send(rank + 1), MPI_particle, &
    !                    partics, n_send, displs, MPI_particle, MPI_COMM_WORLD, ierr)

    
  END SUBROUTINE Make_Tree

End MODULE barnes_hut
