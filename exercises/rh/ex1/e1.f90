program e1

  use geometry
  use particle
  implicit none

INTEGER :: i
INTEGER :: n ! number of particles
type(particle3d), dimension(:), allocatable :: particles ! particle type array (contains the particles info)
real, dimension(:), allocatable :: m

READ *, n ! number of particles
ALLOCATE(particles(n)) 
ALLOCATE(m(n)) 

! probably this part is easier using a file with the particle info, but oky
DO i = 1, n
  PRINT *, 'Insert the info of the particle: '
  READ *, particles(i)%m, particles(i)%p, particles(i)%v !intro values for each particle of mass, position and velocity  
END DO

m = particles(:)%m
print *, m

end program e1 
