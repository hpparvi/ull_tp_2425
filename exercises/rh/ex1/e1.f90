program e1

  use geometry
  use particle
  implicit none

  INTEGER :: i, j ! loop indexes
  INTEGER :: rc 
  INTEGER :: n = 0 ! number of particles
  REAL :: dt, t_end, t, dt_out, t_out = real64 ! time variables
  REAL :: r = real64
  type(particle3d), dimension(:), allocatable :: particles ! particle type array (contains particles info)
  type(vector3d), dimension(:), allocatable :: a ! acceleration array (for each particle) 
  type(vector3d) :: rji
  CHARACTER(len=*), PARAMETER :: filename = 'initial_conditions.dat', outname = 'results.dat'


  OPEN (file = filename, action = 'read', status = 'old', unit = 3, iostat = rc)
  IF (rc/=0) WRITE (*,*) 'Cannot open file ' , filename

  
  READ (3, *) dt
  READ (3, *) t_end
  READ (3, *) dt_out
  READ (3, *) n

  ALLOCATE(particles(n))

  DO i = 1, n
     READ (3, *) particles(i)%m, particles(i)%p%x, particles(i)%p%y, particles(i)%p%z,&
          &particles(i)%v%x, particles(i)%v%y, particles(i)%v%z
  END DO   
  CLOSE(3)

  ALLOCATE(a(n))

  a = vector3d(0.0,0.0,0.0)
  CALL accelerations(particles)
  
  t_out = 0.0
  t = 0.0
  
  OPEN (file = outname, action = 'write', status = 'replace', unit = 4, iostat = rc) 
    IF (rc/=0) WRITE (*,*) 'Cannot open file ' , outname
  
  WRITE(4, *) particles%p
  DO WHILE (t <= t_end) 
    t = t + dt
    
    particles%v = particles%v + a * dt/2. 
    particles%p = particles%p + particles%v * dt 
    
    a = vector3d(0.0,0.0,0.0) 
    CALL accelerations(particles) 
    
    particles%v = particles%v + a * dt/2.
  
    t_out = t_out + dt
    
    IF (t_out >= dt_out) THEN
      WRITE(4, *) particles%p
      t_out = 0.0
    END IF
    
  END DO

  CLOSE(4) 

  
  CONTAINS
  
  SUBROUTINE accelerations(particles)
    TYPE(particle3d), DIMENSION(:), INTENT(in) :: particles ! idk if this is necessary neither :-(
    INTEGER :: i,j ! idk if this is necessary	
    DO i = 1,n
      DO j = i+1,n
        rji = particles(j)%p - particles(i)%p
    	r =  distance(particles(j)%p, particles(i)%p)
    	a(i) = a(i) + particles(j)%m * rji / r**3
    	a(j) = a(j) - particles(i)%m * rji / r**3
      END DO
    END DO

  END SUBROUTINE accelerations
  
end program e1 
