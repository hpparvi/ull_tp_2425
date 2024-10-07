program an0101
  implicit none
  integer :: i, n

  print *, "How many times should we print the text?"
  read *, n

  if (n==1) then
     print *, "Meh"
  end if
  
  do i = 1, n
     print *, "Hello world!"
  end do
  
end program an0101
