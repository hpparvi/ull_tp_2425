PROGRAM tree           !This program simulates the N-body problem using the Barnes-Hut algorithm.
  
  use iso_fortran_env  !This module ensures all variables are defined as 64-bit.
  use omp_lib          !This module is used for parallelization with OpenMP
  use geometry         !This module defines 3D vector and point operations for vector3d and point3d types.
  use definitions      !This module defines key data structures for N-body simulations.
  use barneshut        !This module implements the Barnes-Hut algorithm.
  IMPLICIT NONE

  ! Variables for time measurement.
  integer(int64) :: start, finish, rate
  real :: elapsed_time

  ! Simulation parameters and arrays.
  integer(int64) :: i
  type(vector3d), allocatable :: a(:)           !Acceleration of each particle.
  type(particle3d), allocatable :: p(:)         !Array of particles.
  TYPE (CELL), POINTER :: head, temp_cell       !Head of the octree and a temporary cell pointer.
  type(vector3d) :: rji                         !Vector representing the distance between particles.


  character(len=30) :: data    !Input file containing initial particle data.
  character(len=30) :: orbits  !Output file for storing particle positions over time.

  data = 'data_input.dat'      !Assign input file name.
  orbits = 'data_output.dat'   !Assign output file name.

  ! Read the simulation parameters and initial conditions from the input file.
  call system_clock(count_rate=rate)
  call system_clock(count=start)

  OPEN(10, file=data, status='old', action='read')

  ! Read the simulation time step, output time step, total time, and number of particles.
  READ(10, *) dt
  PRINT *, "The read value of the time step (dt) is ", dt

  READ(10, *) dt_out
  PRINT *, "The read value of the output time step (dt_out) is", dt_out

  READ(10, *) t_end
  PRINT *, "The read value of the total simulation time (t_end) is ", t_end

  READ(10, *) n
  PRINT *, "The read value of the number of bodies (n) is ", n

  ALLOCATE(a(n))  !Allocate memory for acceleration vectors.
  ALLOCATE(p(n))  !Allocate memory for particle structures.

  ! Read positions, velocities, and masses of each particle.
  DO i = 1, n
     READ(10, *) p(i)%m, p(i)%p%x, p(i)%p%y, p(i)%p%z, &            
                                      p(i)%v%x, p(i)%v%y, p(i)%v%z  
     PRINT *, "Particle", i, ":"
     PRINT *, "  Position:", p(i)%p
     PRINT *, "  Velocity:", p(i)%v
     PRINT *, "  Mass:    ", p(i)%m
  END DO

  CLOSE(10)

  ! Initialize the Head Node of the Octree
  ALLOCATE(head)                 
  CALL Calculate_ranges(head, p) !Calculate the spatial range of the octree.
  head%type = 0                  !Initialize the head cell as empty.
  CALL Nullify_Pointers(head)    !Set all subcell pointers to null.

  ! Construct the Initial Octree
  DO i = 1, n
     CALL Find_Cell(head, temp_cell, p(i))  !Find the cell where the particle belongs.
     CALL Place_Cell(temp_cell, p(i), i)    !Place the particle in the appropriate cell.
  END DO
  CALL Borrar_empty_leaves(head)            !Remove empty leaf cells.
  CALL Calculate_masses(head, p)            !Calculate mass and center of mass for each cell.

  ! Calculate Initial Accelerations
  DO i = 1, n
     a(i) = vector3d(0.0, 0.0, 0.0)         
  END DO

  OPEN(11, file=orbits, status='old', action='write')

  !$OMP PARALLEL PRIVATE(i, rji) SHARED(head, p, a)  !Using OMP to parallelize the force calculation.
  CALL Calculate_forces(head, n, p, rji, a)          !Calculate gravitational forces.
  !$OMP END PARALLEL

  ! Main Simulation Loop
  t_out = 0.0
  t = 0.0
  DO WHILE (t <= t_end)
     !$OMP PARALLEL PRIVATE(i, rji) SHARED(head, p, a)
     DO i = 1, n
        p(i)%v = p(i)%v + a(i) * (dt / 2)
     END DO

     DO i = 1, n
        p(i)%p = p(i)%p + p(i)%v * dt
     END DO

     !$OMP SINGLE
     ! Rebuild the octree after positions have changed.
     CALL Borrar_tree(head)
     CALL Calculate_ranges(head, p)
     head%type = 0
     CALL Nullify_Pointers(head)
     DO i = 1, n
        CALL Find_Cell(head, temp_cell, p(i))
        CALL Place_Cell(temp_cell, p(i), i)
     END DO
     CALL Borrar_empty_leaves(head)
     CALL Calculate_masses(head, p)
     
     DO i = 1, n
        a(i) = vector3d(0.0, 0.0, 0.0)
     END DO
     !$OMP END SINGLE

     CALL Calculate_forces(head, n, p, rji, a)  !Recalculate forces after tree update.

     !$OMP DO
     DO i = 1, n
        p(i)%v = p(i)%v + a(i) * (dt / 2)
     END DO
     !$OMP END DO
     !$OMP END PARALLEL

     ! Write output particle positions.
     t_out = t_out + dt
     IF (t_out >= dt_out) THEN
        WRITE(11, '(E12.2)', ADVANCE='NO') t
        DO i = 1, n
           WRITE(11, '(3X, E12.6, 3X, E12.6, 3X, E12.6)', ADVANCE='NO') & 
               p(i)%p%x, p(i)%p%y, p(i)%p%z
        END DO
        WRITE(11, *)
        t_out = 0.0
     END IF

     t = t + dt
  END DO

  CLOSE(11)

  ! Calculate and print the elapsed time of the simulation.
  CALL system_clock(count=finish)
  elapsed_time = real(finish - start, kind=real64) / real(rate, kind=real64)
  IF (elapsed_time >= 60) THEN
     PRINT *, "Elapsed time:", floor(elapsed_time / 60), "minutes", &
              floor(mod(elapsed_time, 60)), "seconds"
  ELSE
     PRINT *, "Elapsed time:", elapsed_time, "seconds"
  END IF

END PROGRAM tree
