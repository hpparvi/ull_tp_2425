PROGRAM ex2
  USE, INTRINSIC :: iso_fortran_env
  use omp_lib
  USE geometry
  USE particle
  USE tree
  IMPLICIT NONE

  TYPE(particle3d), ALLOCATABLE :: parts(:)
  CHARACTER(len=*), PARAMETER :: filename = 'initial_conditions.dat', outname = 'output.dat'
  INTEGER :: n=0
  INTEGER :: i, j, stat
  REAL(real64) :: dt, dt_out, t_end, theta = 1., t_out=0., t=0.
  INTEGER :: start_time, end_time, rate
  REAL :: elapsed_time
  TYPE(vector3d), ALLOCATABLE :: aa(:)

  TYPE(cell), POINTER :: head, temp_cell

  
  !Start timer
  call system_clock(count_rate = rate)
  call system_clock(count = start_time)

  OPEN (file = filename, action = 'read', status = 'old', unit = 3, iostat = stat) !opens input file
  IF (stat/=0) WRITE (*,*) 'Cannot open file ' , filename

  DO i = 1, 3, 1 ! skips the first 3 lines as those are not particles
     READ(3, *)
  END DO

  DO ! this loop counts the amount of particles
     READ(3, *, iostat = stat)
     IF (stat/=0) EXIT
     n = n + 1
  END DO

  REWIND(3) ! goes back to the beginning of the file

  READ (3, *) dt
  READ (3, *) dt_out
  READ (3, *) t_end

  ALLOCATE(parts(n)) !allocates the particles now that it knows how many there are
  ALLOCATE(aa(n))

  DO i = 1, n !reads initial conditions for all particles
     READ (3, *) parts(i)
  END DO
  CLOSE(3)

  !! Initializing the head node
  ALLOCATE(head)
  CALL calculate_ranges(parts, head)
  head%TYPE = 0
  CALL nullify_pointers(head)

  !! Creating the starting tree
  DO i = 1, n
     CALL find_cell(head, temp_cell, parts(i))
     CALL place_cell(temp_cell, parts(i), i)
  END DO

  CALL delete_empty_leaves(head)
  CALL calculate_masses(head)

  !! Compute initial accelerations
  aa = vector3d(0.,0.,0.)
  CALL calculate_forces(head, aa, parts, theta)

  !! Open output file
  OPEN (file = outname, action = 'write', status = 'replace', unit =&
       & 4, iostat = stat)
  IF (stat/=0) WRITE (*,*) 'Cannot open file ', outname

  !! Main loop
  t_out = 0.0

  DO WHILE (t .LE. t_end)

     !$omp parallel shared(parts, aa, dt)
     !$omp workshare
     parts%v = parts%v + aa * (dt/2.)
     parts%p = parts%p + parts%v * dt
     !$omp end workshare
     !$omp single

     !! Redo the tree
     call delete_tree(head)

     call calculate_ranges(parts, head)
     head%type = 0
     call nullify_pointers(head)

     DO i = 1, n
        CALL find_cell(head, temp_cell, parts(i))
        CALL place_cell(temp_cell, parts(i), i)
     END DO
     
     
     CALL delete_empty_leaves(head)
     CALL calculate_masses(head)

     aa = vector3d(0.,0.,0.)
     call calculate_forces(head, aa, parts, theta)

     !$omp end single
     !$omp workshare
     parts%v = parts%v + aa * (dt/2.)
     !$omp end workshare
     !$omp end parallel
     t_out = t_out + dt
     IF (t_out .GE. dt_out) THEN
        WRITE(4, fmt='(F11.3)', advance='no') t
        DO i = 1, n
           WRITE(4, fmt='(4ES13.4)', advance='no') parts(i)%p
        END DO
        WRITE(4,*) ''
        t_out = 0.0
     END IF
     t = t + dt

  END DO

  CLOSE(4)

  call system_clock(count = end_time)

  elapsed_time = real(end_time-start_time)/real(rate)

  print *, 'Elapsed time: ', elapsed_time, ' seconds'

END PROGRAM ex2
