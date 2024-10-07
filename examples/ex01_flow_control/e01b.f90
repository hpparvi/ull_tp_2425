program an0102
  implicit none
  integer :: i, n, f = 1

  print *, "Give n:"
  read *, n

  do i = 1, n
     f = f * i
     print *, f
  end do
  
end program an0102
