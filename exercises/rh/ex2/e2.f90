program e2 
  use iso_fortran_env
  use geometry
  use particle
  IMPLICIT NONE
  
  INTEGER :: i,j,k,n,rc
  REAL (real64):: dt, t_end, t, dt_out, t_out
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: particles !p,v,m position, velocity and mass of each particles
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: acc !acceleration
  CHARACTER(len=*), PARAMETER :: filename = 'initial_conditions.dat', outname = 'results.dat' ! i.c. input/output files names
  
  TYPE RANGE
    REAL, DIMENSION(3) :: min_val,max_val ! it should be a vector, shouldn't be?
  END TYPE RANGE
  
  ! new type cell as a pointer
  TYPE CPtr
    TYPE(CELL), POINTER :: ptr !pointer of the 'target' cell
  END TYPE CPtr
 
  TYPE CELL
    TYPE (RANGE) :: range
    TYPE(particle3d) :: part!(particle info in the cell)
    INTEGER :: pos  !! id of the particle
    INTEGER :: type !! 0 = no particle; 1 = particle; 2 = conglomerado
    TYPE(point3d) :: c_o_m ! center of mass of that cell
    TYPE (CPtr), DIMENSION(2,2,2) :: subcell ! subcells in main cell
    ! that contains cells as well 
  END TYPE CELL

  TYPE (CELL), POINTER :: head, temp_cell ! create cell (as pointer)
  
  
  !! Lectura de datos
  !!!!!!!!!!!!!!!!!!!
  ! open the input file
  OPEN (file = filename, action = 'read', status = 'old', unit = 3, iostat = rc)
  IF (rc/=0) WRITE (*,*) 'Cannot open file ' , filename  
  ! include read initial_conditions file
  READ (3, *) dt
  READ (3, *) dt_out
  READ (3, *) t_end
  READ (3, *) n
  ALLOCATE(particles(n))
  ALLOCATE(acc(n))
  DO i = 1, n
    READ (3, *) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z, &
          & particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
  END DO
  
  
  !! Inicialización head node
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(head)
  CALL Calculate_ranges(head) 
  head%type = 0 ! no particle
  CALL Nullify_Pointers(head) ! null all the pointers first just in case
  
  !! Creación del árbol inicial
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print*, n
  DO i = 1,n
    CALL Find_Cell(head,temp_cell,particles(i)) 
    CALL Place_Cell(temp_cell,particles(i),i)

  END DO

  CALL Borrar_empty_leaves(head)
  CALL Calculate_masses(head)
  
  
  !! Calcular aceleraciones iniciales
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    acc = vector3d(0.0,0.0,0.0)
    CALL Calculate_forces(head)
    
  ! open the output file 
  OPEN (file = outname, action = 'write', status = 'replace', unit = 4, iostat = rc) 
    IF (rc/=0) WRITE (*,*) 'Cannot open file ' , outname
  
  ! the first record of the output file is the initial positions of the particles
  WRITE(4, *) particles%p
  
  !! Bucle principal
  !!!!!!!!!!!!!!!!!!
    t_out = 0.0
    DO  WHILE (t <= t_end)
      particles%v = particles%v + acc * dt/2
      particles%p = particles%p + particles%v * dt

!! Las posiciones han cambiado, por lo que tenemos que borrar
!! y reinicializar el árbol
      CALL Borrar_tree(head) !remove previous tree
      CALL Calculate_ranges(head) !calculate head range again
      head%type = 0 !return to 0 head type
      CALL Nullify_Pointers(head) 
      
      DO i = 1,n
        CALL Find_Cell(head,temp_cell,particles(i))
        CALL Place_Cell(temp_cell,particles(i),i)
      END DO
      
      CALL Borrar_empty_leaves(head)
      CALL Calculate_masses(head)

      acc = vector3d(0.0,0.0,0.0)
      CALL Calculate_forces(head)
      
      particles%v = particles%v + acc * (dt/2.)
      
      t_out = t_out + dt
      IF (t_out >= dt_out) THEN
        WRITE(4, *) particles%p ! positions in one row (one particle position after another)
          DO i = 1,n
          PRINT*, particles(i)%p
         END DO
        PRINT*, "-----------------------------------"
        PRINT*, ""
        t_out = 0.0
      END IF
      t = t + dt
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
    TYPE(CELL), POINTER :: goal
    REAL(real64), DIMENSION(3):: mins,maxs,medios !should be vectors? 
    REAL (real64) :: span

    mins(1) = MINVAL(particles%p%x) 
    maxs(1) = MAXVAL(particles%p%x)
    mins(2) = MINVAL(particles%p%y) 
    maxs(2) = MAXVAL(particles%p%y)
    mins(3) = MINVAL(particles%p%z) 
    maxs(3) = MAXVAL(particles%p%z)
    
! Al calcular span le sumo un 10% para que las
! particulas no caigan justo en el borde

    span = MAXVAL(maxs - mins) * 1.1
    medios = (maxs + mins)/2. 
    goal%range%min_val = medios - span/2.  
    goal%range%max_val = medios + span/2. 
  END SUBROUTINE Calculate_Ranges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Find_Cell
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Encuentra la celda donde colocaremos la particula.
!! Si la celda que estamos considerando (root) no tiene
!! particula o tiene una particula, es esta celda donde
!! colocaremos la particula (goal).

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
    type(particle3d) :: part
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
    TYPE(particle3d) :: part
    INTEGER :: n
    print *, goal%type
    print *, goal%pos
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
    TYPE(particle3d) :: part
    INTEGER :: i,j,k,n
    INTEGER, DIMENSION(3) :: octant !stays as an array because just aim which octant is
    part = goal%part
    goal%type=2
    DO i = 1,2
      DO j = 1,2
        DO k = 1,2
          octant = (/i,j,k/) !???
          ALLOCATE(goal%subcell(i,j,k)%ptr)
          goal%subcell(i,j,k)%ptr%range%min_val = Calcular_Range(0,goal,octant)
          goal%subcell(i,j,k)%ptr%range%max_val = Calcular_Range(1,goal,octant)
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
    TYPE(particle3d) :: part
    TYPE(CELL), POINTER :: goal
    LOGICAL :: Belongs
    IF ((part%p%x >= goal%range%min_val(1)) .AND. &
      (part%p%x < goal%range%max_val(1)) .AND. &
      (part%p%y >= goal%range%min_val(2)) .AND. &! change bc if not if a particle is in the middle u have issues to assign a subcell
      (part%p%y < goal%range%max_val(2)) .AND. &
      (part%p%z >= goal%range%min_val(3)) .AND. &
      (part%p%z < goal%range%max_val(3))) THEN
      Belongs = .TRUE.
      print *, 'found it'
    
    ELSE
    print *, 'not here bro'
      Belongs = .FALSE.
    END IF
  END FUNCTION Belongs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calcular_Range
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Dado un octante "octant" (1,1,1, 1,1,2 ... 2,2,2),
!! calcula sus rangos en base a los rangos de
!! "goal". Si "what" = 0 calcula los minimos. Si what=1
!! calcula los maximos.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION Calcular_Range (what,goal,octant)
    INTEGER :: what
    TYPE(CELL), POINTER :: goal
    INTEGER, DIMENSION(3) :: octant
    INTEGER, DIMENSION(3) :: Calcular_Range, valor_medio
    valor_medio = (goal%range%min_val + goal%range%max_val) / 2.

    SELECT CASE (what)
    CASE (0)
      WHERE (octant == 1)
        Calcular_Range = goal%range%min_val
      ELSEWHERE
        Calcular_Range = valor_medio
      END WHERE
    
    CASE (1)
      WHERE (octant == 1)
        Calcular_Range = valor_medio
      ELSEWHERE
        Calcular_Range = goal%range%max_val
      END WHERE
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
    REAL(real64) :: mass
    type(point3d) :: c_o_m
    goal%part%m = 0
    goal%c_o_m = point3d(0.,0.,0.)
    SELECT CASE (goal%type)
    CASE (1) !only one particle
    goal%part%m = particles(goal%pos)%m !is the mass of the particle in the cell 
    goal%c_o_m= particles(goal%pos)%p !is position of the particle
    CASE (2) !in case of a conglomerate
    DO i = 1,2
      DO j = 1,2
        DO k = 1,2
          IF (ASSOCIATED(goal%subcell(i,j,k)%ptr)) THEN
            CALL Calculate_masses(goal%subcell(i,j,k)%ptr)
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
  END SUBROUTINE Calculate_masses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    INTEGER :: i !j,k,start,end
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
    REAL(real64) :: l,D,r
    TYPE(vector3d) :: rji
    
    SELECT CASE (tree%type)
    CASE (1) ! only one particle in the tree
    IF (goal .NE. tree%pos) THEN
      rji = tree%c_o_m - particles(goal)%p !distance between the com and particle
      r = distance(tree%c_o_m, particles(goal)%p)
    !r2 = SUM(rji**2) !
    !r3 = r2 * SQRT(r2)
      acc(goal) = acc(goal) + particles(tree%pos)%m * rji / r**3
    END IF
    
    CASE (2)
!! El rango tiene el mismo span en las 3 dimensiones
!! por lo que podemos considerar una dimension cualquiera
!! para calcular el lado de la celda (en este caso la
!! dimension 1)
    l = tree%range%max_val(1) - tree%range%min_val(1) 
    rji = tree%c_o_m - particles(goal)%p
    r = distance(tree%c_o_m, particles(goal)%p)
    ! r2 = SUM(rji**2)
    !D = r SQRT(r2)
    IF (l/r < 1.) THEN
!! Si conglomerado, tenemos que ver si se cumple l/D < theta
      !r3 = r**2 * D 
      acc(goal) = acc(goal) + tree%part%m * rji / r**3
      
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


end program e2
