PROGRAM ex2

  !$ USE omp_lib
  USE geometry
  USE particle
  USE bh
  IMPLICIT NONE

  
  integer :: start_count, end_count, count_rate, n, i
  REAL ::  dt, t_end, t, dt_out, t_out, execution_time, r2, r3
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: par
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a
  REAL(8), PARAMETER :: theta = 1.0d0  ! Usando REAL(8) para precisión doble


  
  !! Lectura de datos
  !!!!!!!!!!!!!!!!!!!
  print*, "Enter value of the timestep, dt: "
  read*, dt
  print*, "Enter value of the output time, dt_out: "
  read*, dt_out
  print*, "Enter value of the final time, t_end: "
  read*, t_end
  print*, "Enter value of the number of particles, n: "
  read*, n
  
  ALLOCATE(par(n))
	ALLOCATE(a(n))
  
  ! Open the file for reading
  open(unit=10, file="IC.txt", status='old')
  ! Read and process particle data
  do i = 1, n
   read(10, *) par(i)%m, par(i)%p, par(i)%v
  end do
  close(unit=10)
  
  ! Measure the starting time
  call system_clock(count=start_count, count_rate=count_rate)

  
  !! Inicializaci´on head node
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(head)
  
  CALL Calculate_ranges(head, par)
  head%type = 0
  CALL Nullify_Pointers(head)
  
  !! Creaci´on del ´arbol inicial
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1,n
    CALL Find_Cell(head,temp_cell,par(i))
    CALL Place_Cell(temp_cell, par(i), i)
  END DO

  
  CALL Borrar_empty_leaves(head)
  CALL Calculate_masses(head, par)
  
  !! Calcular aceleraciones iniciales
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  a=vector3d(0.0,0.0,0.0)
  CALL Calculate_forces(head,par,n,theta) ! Parallelized in bh module (loops Calculate_forces_aux)
  
  !! Fichero para guardar resultados
  open(unit=11, file="output.dat", status='replace')
  
  t_out = 0.0

  !! Bucle principal
  !!!!!!!!!!!!!!!!!!

  DO WHILE (t <= t_end)
		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
    !$OMP DO
    DO i=1,n
      par(i)%v = par(i)%v + mulrv(REAL(dt/2,KIND=8), a(i))
   	  par(i)%p = par(i)%p + mulrv(REAL(dt,KIND=8), par(i)%v) 
   	END DO
    !$OMP END DO
    !$OMP END PARALLEL

    !! Las posiciones han cambiado, por lo que tenemos que borrar
    !! y reinicializar el ´arbol
    
    CALL Borrar_tree(head)  
    CALL Calculate_ranges(head, par)
    head%type = 0
    CALL Nullify_Pointers(head)

    DO i = 1,n
      CALL Find_Cell(head,temp_cell,par(i))
      CALL Place_Cell(temp_cell,par(i),i)
    END DO
    
    CALL Borrar_empty_leaves(head)
    CALL Calculate_masses(head, par)
    
    a=vector3d(0.0,0.0,0.0)
      
    CALL Calculate_forces(head,par,n,theta) ! Parallelized in bh module (loops Calculate_forces_aux)
    
 		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
		!$OMP DO
		DO i=1,n
 	     par(i)%v = par(i)%v + mulrv(REAL(dt/2,KIND=8), a(i))
		END DO
		!$OMP END DO
		!$OMP END PARALLEL
    
    t_out = t_out + dt
    
    IF (t_out >= dt_out) THEN
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
			!$OMP DO
      DO i = 1,n
        WRITE(11,*) i, par(i)%p
      END DO
 		  !$OMP END DO
		  !$OMP END PARALLEL
        
      t_out = 0.0
    END IF
    t = t + dt
  END DO

  close(unit=11)
  
  ! Measure the ending time
  call system_clock(count=end_count)

  ! Calculate the execution time
  execution_time = real(end_count - start_count) / real(count_rate)

  print*, "Total execution time:", execution_time, "seconds"
  
END PROGRAM ex2PROGRAM ex2

  !$ USE omp_lib
  USE geometry
  USE particle
  USE bh
  IMPLICIT NONE

  
  integer :: start_count, end_count, count_rate, n, i
  REAL ::  dt, t_end, t, dt_out, t_out, execution_time, r2, r3
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: par
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: a

  
  !! Lectura de datos
  !!!!!!!!!!!!!!!!!!!
  print*, "Enter value of the timestep, dt: "
  read*, dt
  print*, "Enter value of the output time, dt_out: "
  read*, dt_out
  print*, "Enter value of the final time, t_end: "
  read*, t_end
  print*, "Enter value of the number of particles, n: "
  read*, n
  
  ALLOCATE(par(n))
	ALLOCATE(a(n))
  
  ! Open the file for reading
  open(unit=10, file="IC.txt", status='old')
  ! Read and process particle data
  do i = 1, n
   read(10, *) par(i)%m, par(i)%p, par(i)%v
  end do
  close(unit=10)
  
  ! Measure the starting time
  call system_clock(count=start_count, count_rate=count_rate)

  
  !! Inicializaci´on head node
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(head)
  
  CALL Calculate_ranges(head, par)
  head%type = 0
  CALL Nullify_Pointers(head)
  
  !! Creaci´on del ´arbol inicial
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1,n
    CALL Find_Cell(head,temp_cell,par(i))
    CALL Place_Cell(temp_cell, par(i), i)
  END DO

  
  CALL Borrar_empty_leaves(head)
  CALL Calculate_masses(head, par)
  
  !! Calcular aceleraciones iniciales
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  a=vector3d(0.0,0.0,0.0)
  CALL Calculate_forces(head,par,a,n,theta) ! Parallelized in bh module (loops Calculate_forces_aux)
  
  !! Fichero para guardar resultados
  open(unit=11, file="output.dat", status='replace')
  
  t_out = 0.0

  !! Bucle principal
  !!!!!!!!!!!!!!!!!!

  DO WHILE (t <= t_end)
		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
    !$OMP DO
    DO i=1,n
      par(i)%v = par(i)%v + mulrv(REAL(dt/2,KIND=8), a(i))
   	  par(i)%p = par(i)%p + mulrv(REAL(dt,KIND=8), par(i)%v) 
   	END DO
    !$OMP END DO
    !$OMP END PARALLEL

    !! Las posiciones han cambiado, por lo que tenemos que borrar
    !! y reinicializar el ´arbol
    
    CALL Borrar_tree(head)  
    CALL Calculate_ranges(head, par)
    head%type = 0
    CALL Nullify_Pointers(head)

    DO i = 1,n
      CALL Find_Cell(head,temp_cell,par(i))
      CALL Place_Cell(temp_cell,par(i),i)
    END DO
    
    CALL Borrar_empty_leaves(head)
    CALL Calculate_masses(head, par)
    
    a=vector3d(0.0,0.0,0.0)
      
    CALL Calculate_forces(head,par,a,n,theta) ! Parallelized in bh module (loops Calculate_forces_aux)
    
 		!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
		!$OMP DO
		DO i=1,n
 	     par(i)%v = par(i)%v + mulrv(REAL(dt/2,KIND=8), a(i))
		END DO
		!$OMP END DO
		!$OMP END PARALLEL
    
    t_out = t_out + dt
    
    IF (t_out >= dt_out) THEN
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
			!$OMP DO
      DO i = 1,n
        WRITE(11,*) i, par(i)%p
      END DO
 		  !$OMP END DO
		  !$OMP END PARALLEL
        
      t_out = 0.0
    END IF
    t = t + dt
  END DO

  close(unit=11)
  
  ! Measure the ending time
  call system_clock(count=end_count)

  ! Calculate the execution time
  execution_time = real(end_count - start_count) / real(count_rate)

  print*, "Total execution time:", execution_time, "seconds"
  
END PROGRAM ex2
