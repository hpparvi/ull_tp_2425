program ex1
  use geometry
  use particle
  implicit none
  real :: dt, t_end, dt_out, t_out
  real :: t
  type(particle3d), allocatable :: particles(:)
  integer :: rc, stat, i, j
  integer :: n = 0
  character(len=*), parameter :: filename = 'data.txt'
  type(vector3d), dimension(:), allocatable :: aa
  type(vector3d) :: rji
  real :: r2

  

  open (file = filename, action = 'read', status = 'old', unit = 3, iostat = rc)
  if (rc/=0) write (*,*) 'Cannot open file ' , filename

  do i = 1, 3, 1 ! skips the first 3 lines as those are not particles
     read(3, *)
  end do
  
  do ! this loop counts the amount of particles
     read(3, *, iostat = stat)
     if (stat/=0) exit
     n = n + 1
  end do
  
  rewind(3)
  
  read (3, *) dt
  read (3, *) dt_out
  read (3, *) t_end

  allocate(particles(n))

  do i = 1, n
     read (3, *) particles(i)%m, particles(i)%p%xx, particles(i)%p%yy, particles(i)%p%zz,&
          &particles(i)%v%xx, particles(i)%v%yy, particles(i)%v%zz
  end do   
  close(3)

  allocate(aa(n))
  aa = vector3d(0.,0.,0.)
  
  
  do i = 1, n
     do j = i+1, n
        rji =  normalize(particles(i)%p - particles(j)%p)
        r2 = distance(particles(i)%p, particles(j)%p)
        aa(i) = aa(i)
     end do
  end do
  
    
end program ex1
