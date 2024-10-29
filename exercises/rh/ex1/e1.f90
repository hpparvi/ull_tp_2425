program e1
  ! Import the modules that we are going to use
  !use iso_fortran_env
  use geometry
  use particle
  implicit none

  INTEGER :: i, j ! loop indexes
  INTEGER :: rc ! variable to read line of input file (with the i.c.)
  INTEGER :: n = 0 ! number of particles
  REAL(kind = 8) :: dt, t_end, t, dt_out, t_out ! time variables
  type(particle3d), dimension(:), allocatable :: particles ! particle type array (contains particles info)
  type(vector3d), dimension(:), allocatable :: acc ! acceleration array (for each particle) 
  CHARACTER(len=*), PARAMETER :: filename = 'initial_conditions.dat', outname = 'results.dat' ! i.c. input/output files names

  ! open the input file
  OPEN (file = filename, action = 'read', status = 'old', unit = 3, iostat = rc)
  IF (rc/=0) WRITE (*,*) 'Cannot open file ' , filename

  ! read the info of the file:
  READ (3, *) dt ! time step for the integration
  READ (3, *) t_end ! time when the integration stops
  READ (3, *) dt_out ! time step that controls if a value is saved in the output file
  
  ! read the number of particles and allocate the array dimension
  READ (3, *) n 
  ALLOCATE(particles(n))

  ! read the info of all the particles
  ! each row contains the info of one particles
  ! the data are in columns: first mass, then position and velocity
  ! (m x y z vx vy vz)
  DO i = 1, n
     READ (3, *) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z,&
          &particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
  END DO   
  CLOSE(3)
  
  ! using the variable n we allocate also the acceleration array
  ALLOCATE(acc(n))
  
  ! initially the acceleration array is null
  acc = vector3d(0.0,0.0,0.0)
  ! calculate the acceleration using the accelerations subroutine
  CALL accelerations(particles,acc)
  
  t_out = 0.0 ! time counter to make a new record in the output file
  t = 0.0 ! time counter of the integration
  
  ! open the output file 
  OPEN (file = outname, action = 'write', status = 'replace', unit = 4, iostat = rc) 
    IF (rc/=0) WRITE (*,*) 'Cannot open file ' , outname
  
  ! the first record of the output file is the initial position of the particle
  WRITE(4, *) particles%p
  
  ! loop for the integration, stops when t reach the ending time
  DO WHILE (t <= t_end) 
    t = t + dt ! add a time step to the time counter
    
    ! calculate the movement of the particles using leapfrog algorithm:  
    particles%v = particles%v + acc * (dt/2.0)  
    particles%p = particles%p + particles%v * dt 
    
    acc = vector3d(0.0,0.0,0.0) 
    CALL accelerations(particles, acc) 
    
    particles%v = particles%v + acc * (dt/2.0)
  
    t_out = t_out + dt ! add to the output counter
    
    ! if dt_out is greater than dt_out, we save the position of the particles
    IF (t_out >= dt_out) THEN
      WRITE(4, *) particles%p ! positions in one row (one particle position after another)
      t_out = 0.0 ! put the counter to zero 
    END IF
    
  END DO
  
  ! close the file
  CLOSE(4) 

  
  CONTAINS
  
  ! subroutine to calculate the accelerations
  SUBROUTINE accelerations(particles,a)
    TYPE(particle3d), DIMENSION(:), INTENT(in) :: particles 
    TYPE(vector3d), DIMENSION(:), INTENT(inout) :: a
    TYPE(vector3d) :: rji ! vector that goes from one particle to another 
    REAL(kind = 8) :: r ! distance between particles
    INTEGER :: i,j ! loop indexes
    
    ! loop to calculate the effect of each particle in the acceleration of themselves	
    DO i = 1,n
      DO j = i+1,n
        rji = particles(j)%p - particles(i)%p
    	r =  distance(particles(j)%p, particles(i)%p)
    	a(i) = a(i) + particles(j)%m * rji / r**3 ! calculate attraction made by other particles 
    	a(j) = a(j) - particles(i)%m * rji / r**3 ! calculate attraction of that particle to the rest of them
      END DO
    END DO

  END SUBROUTINE accelerations
  
end program e1 
