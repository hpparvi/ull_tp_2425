program e2 
  !$ Use omp_lib 
  use, intrinsic ::  iso_fortran_env
  use geometry
  use particle
  use tree_algorithm
  IMPLICIT NONE
  
  INTEGER :: i,j,k,n,rc
  REAL (real64):: dt, t_end, t, dt_out, t_out
  REAL(real64), PARAMETER :: theta = 1
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: particles !p,v,m position, velocity and mass of each particles
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: acc !acceleration
  CHARACTER(len=*), PARAMETER :: filename = 'init_files/example.dat', outname = 'output.dat' ! i.c. input/output files names
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
  CLOSE(UNIT=3)
  
  !! Inicialización head node
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(head)
  CALL Calculate_ranges(head, particles) 
  head%type = 0 ! no particle
  CALL Nullify_Pointers(head) ! null all the pointers first just in case
  
  !! Creación del árbol inicial
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1,n
    CALL Find_Cell(head,temp_cell,particles(i)) 
    CALL Place_Cell(temp_cell,particles(i),i)

  END DO

  CALL Borrar_empty_leaves(head)
  CALL Calculate_masses(head, particles)
  
  
  !! Calcular aceleraciones iniciales
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    acc = vector3d(0.0,0.0,0.0)
    CALL Calculate_forces(head,particles,acc)
    
  ! open the output file 
  OPEN (file = outname, action = 'write', status = 'replace', unit = 4, iostat = rc) 
    IF (rc/=0) WRITE (*,*) 'Cannot open file ' , outname
  
  ! the first record of the output file is the initial positions of the particles
  WRITE(4, *) t, particles%p
  
  !! Bucle principal
  !!!!!!!!!!!!!!!!!!
    t_out = 0.0
    DO  WHILE (t <= t_end)
    
      particles%v = particles%v + acc * (dt/2.)
      particles%p = particles%p + particles%v * dt

!! Las posiciones han cambiado, por lo que tenemos que borrar
!! y reinicializar el árbol
      CALL Borrar_tree(head) !remove previous tree
      CALL Calculate_ranges(head, particles) !calculate head range again
      head%type = 0 !return to 0 head type
      CALL Nullify_Pointers(head) 
      
      DO i = 1,n
        CALL Find_Cell(head,temp_cell,particles(i))
        CALL Place_Cell(temp_cell,particles(i),i)
      END DO
      
      CALL Borrar_empty_leaves(head)
      CALL Calculate_masses(head, particles)

      acc = vector3d(0.0,0.0,0.0)
      CALL Calculate_forces(head,particles,acc)
      
      !$OMP WORKSHARE
      particles%v = particles%v + acc * (dt/2.)
      !$OMP END WORKSHARE []
      t_out = t_out + dt
      IF (t_out >= dt_out) THEN
        WRITE(4, *) t, particles%p ! positions in one row (one particle position after another)
        t_out = 0.0
      END IF
      t = t + dt
    END DO 

  CLOSE(UNIT=4)

end program e2
