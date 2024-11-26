program e2
  ! import required libraries/modules 
  !$ Use omp_lib 
  use, intrinsic ::  iso_fortran_env
  use geometry
  use particle
  use tree_algorithm
  IMPLICIT NONE

  INTEGER :: i,j,k,n,rc
  REAL (real64) :: dt, t_end, t, dt_out, t_out
  INTEGER :: start_time, end_time, count_rate
  REAL :: elapsed_time
  REAL(real64), PARAMETER :: theta = 1
  
  TYPE(particle3d), DIMENSION(:), ALLOCATABLE :: particles 
  TYPE(vector3d), DIMENSION(:), ALLOCATABLE :: acc !acceleration
  CHARACTER(len=*), PARAMETER :: filename = 'particle_position.dat', outname = 'output.dat' ! i.c. input/output files names
  TYPE (CELL), POINTER :: head, temp_cell ! create cell (as pointer)
  
  ! open the input file
  OPEN (file = filename, action = 'read', status = 'old', unit = 3, iostat = rc)
  IF (rc/=0) WRITE (*,*) 'Cannot open file ' , filename  
  ! save the initial conditions
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
  
  !! Initialise head node
  ALLOCATE(head)
  CALL Calculate_ranges(head, particles) 
  head%type = 0 ! no particle
  CALL Nullify_Pointers(head) ! null all the pointers first just in case
  
  ! Start the timer
  call system_clock(start_time, count_rate)

  ! Create initial tree
  DO i = 1,n
    CALL Find_Cell(head,temp_cell,particles(i)) 
    CALL Place_Cell(temp_cell,particles(i),i)

  END DO
  
  CALL Borrar_empty_leaves(head)
  CALL Calculate_masses(head, particles)

  ! Calculate initial accelerations
    acc = vector3d(0.0,0.0,0.0)
    CALL Calculate_forces(head,particles,acc)
    
  ! open the output file 
  OPEN (file = outname, action = 'write', status = 'replace', unit = 4, iostat = rc) 
    IF (rc/=0) WRITE (*,*) 'Cannot open file ' , outname
  
  ! the first record of the output file is the initial positions of the particles
  WRITE(4, *) t, particles%p
  
  !! Main loop 
  !!!!!!!!!!!!!!!!!!
    t_out = 0.0
    DO  WHILE (t <= t_end)
    
      !$omp parallel
      ! update velocities and positions of the particles
      !$OMP WORKSHARE
      particles%v = particles%v + acc * (dt/2.)
      particles%p = particles%p + particles%v * dt
      !$OMP END WORKSHARE
      
      !The positions have changed, so we have to remove and initialise the tree again
      !$OMP SINGLE
      CALL Borrar_tree(head) ! remove previous tree
      CALL Calculate_ranges(head, particles) ! calculate head range again
      head%type = 0 
      CALL Nullify_Pointers(head) 
      
      DO i = 1,n
        CALL Find_Cell(head,temp_cell,particles(i))
        CALL Place_Cell(temp_cell,particles(i),i)
      END DO
      
      CALL Borrar_empty_leaves(head)
      CALL Calculate_masses(head, particles)
      
      acc = vector3d(0.0,0.0,0.0)
      !$OMP END SINGLE 
      CALL Calculate_forces(head,particles,acc)
      
      !$OMP WORKSHARE
      particles%v = particles%v + acc * (dt/2.)
      !$OMP END WORKSHARE
      !$omp end parallel
      
      t_out = t_out + dt
      IF (t_out >= dt_out) THEN
        WRITE(4, *) t, particles%p ! time and positions in one row (one particle position after another)
        t_out = 0.0
      END IF
      t = t + dt
    END DO 

  CLOSE(UNIT=4)
  
    ! End the timer
  call system_clock(end_time)

  ! Calculate elapsed time
  elapsed_time = REAL(end_time-start_time) / REAL(count_rate)
  
  ! Output the elapsed time
  PRINT *, "The elapsed time is ", &
       elapsed_time, " seconds"

end program e2
