program ex1
  use geometry
  use particle
  implicit none
  integer, parameter :: n = 3
  real, parameter :: dt = 0.001, t_end = 100, dt_out = 0.1, t_out = 100
  real :: t
  type(particle3d) :: star1

  

  open (file = 'file.txt', action = 'read', status = 'old', unit = 3)

  read (3,*) star1%m, star1%p%xx, star1%p%yy, star1%p%zz, star1%v%xx, star1%v%yy, star1%v%zz

  close(3)
    
end program ex1
