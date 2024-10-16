program test
  use slist
  implicit none

  type(sortedlist) :: list

  call init_cell(list%head, "river")
  call extend(list, "spark")
  call extend(list, "melody")
  call extend(list, "whisper")
  call extend(list, "canyon")
  call extend(list, "drift")
  call extend(list, "lantern")
  call extend(list, "echo")
  call extend(list, "quartz")
  call extend(list, "breeze")
  call print_forward(list%head)
  print *
  call print_reverse(list%head)

end program test
