PROGRAM ex2

  USE geometry
  USE ex2_deriv_types
  USE barnes_hut
  USE ISO_FORTRAN_ENV
  !$ use omp_lib
  IMPLICIT NONE      
  
  REAL(REAL64) :: dt, t_end, t, dt_out, t_out    
  character(len=25) :: input, output, orbit
  

  ! To test time of execution
  INTEGER :: count_begin, count_end, count_rate
  
  
!! Data reading
!!!!!!!!!!!!!!!!!!!

  input = "initial_conditions.dat"
  output = "output.dat"
  
  OPEN (1,file = input, status = 'old', action = 'read')
  

  READ(1,*) dt
  READ(1,*) dt_out
  READ(1,*) t_end
  READ(1,*) n
  ALLOCATE(pt(n))
  ALLOCATE(a(n))
  DO i = 1, n
    READ(1,*) pt(i)%m, pt(i)%p, pt(i)%v
  END DO

  CLOSE (UNIT=1) 
  
  ! Start measuring time here
  CALL SYSTEM_CLOCK(count_begin, count_rate)
  !PRINT*, "The count rate is ", count_rate
  
  
!! Initialization of head node
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ALLOCATE(head)
  CALL Calculate_ranges(head)
  head%type = 0
  CALL Nullify_Pointers(head)
  
  
!! Creation of the initial tree
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO i = 1,n
    CALL Find_Cell(head,temp_cell,pt(i)%p)
    CALL Place_Cell(temp_cell,pt(i)%p,i)
  END DO
  CALL Borrar_empty_leaves(head)
  CALL Calculate_masses(head)


!! Calculate initial acceleration
!! Parallelization will be applied, and when it is parallelized, 
!! a message will be show with the number of cores.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  a = vector3d(0,0,0)
  !$omp parallel private(i) shared(head,temp_cell,pt,a)  
  !$ IF (omp_get_thread_num() == 0) THEN
  !$   PRINT'(A,I2,A)', 'Parallelized process with ', omp_get_num_threads(), ' cores' 
  !$ END IF
  CALL Calculate_forces(head)
  !$omp end parallel

!! Main loop
!!!!!!!!!!!!!!!!!!
  
  OPEN(3, file = output, status = 'replace', action = 'write')
  t = 0.0
  t_out = 0.0
  DO 
    t = t + dt
    t_out = t_out + dt
    pt%v = pt%v + a * (dt/2)
    pt%p = pt%p + pt%v * dt
    
    !! The positions have changed, so the tree needs to be reset.

    CALL Borrar_tree(head)
    CALL Calculate_ranges(head)
    head%type = 0
    CALL Nullify_Pointers(head)
    DO i = 1,n
      CALL Find_Cell(head,temp_cell,pt(i)%p)
      CALL Place_Cell(temp_cell,pt(i)%p,i)
    END DO
    CALL Borrar_empty_leaves(head)
    CALL Calculate_masses(head)
    a = vector3d(0,0,0)
    
    !! The calculation of forces is parallelized.
    
    !$omp parallel private(i) shared(head,temp_cell,pt,a)
    CALL Calculate_forces(head)
    !$omp workshare
    pt%v = pt%v + a * (dt/2)
    !$omp end workshare
    !$omp end parallel
    
    IF (t_out >= dt_out) THEN
      write (3, '(3F12.4)', advance='no') t
      DO i = 1,n
        IF(i<n) THEN
          write (3, '(6F16.6)', advance='no') pt(i)%p
        ELSE
          write (3, '(6F16.6)') pt(i)%p
        END IF
      END DO
      t_out = 0.0 ! reset of output parameter 
    END IF
    ! Exit conditions
    IF (t>= t_end) THEN
      EXIT
    END IF
  END DO
  
  CLOSE (3)
  !! End measuring time here
  CALL SYSTEM_CLOCK(count_end)

  PRINT'(A,F7.3,A)', "The elapsed time is ", &
       REAL(count_end-count_begin) / REAL(count_rate), &
       " seconds"  
 
END PROGRAM ex2