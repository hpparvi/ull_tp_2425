module slist
  implicit none
  
  type :: cell
     type(cell), pointer :: next => null()
     character(:), allocatable :: d
  end type cell
  
  type :: sortedlist
     integer :: length = 0
     type(cell), pointer :: head => null()
  end type sortedlist

contains
  
   subroutine init_cell(n, d)
    type(cell), pointer, intent(out) :: n
    character(len=*) :: d
    allocate(n)
    n%d = d
  end subroutine init_cell

  subroutine extend(list, s)
    type(sortedlist) :: list
    type(cell), pointer :: cur, new
    character(len=*) :: s

    call init_cell(new, s)
    
    cur => list%head
    do while (associated(cur%next))
       if (new%d(1:1) > cur%next%d(1:1)) then
          cur => cur%next
       else
          exit
       end if
    end do

    if (associated(cur, list%head)) then
       if (new%d(1:1) < cur%d(1:1)) then
          new%next => cur
          list%head => new
       else
          new%next => cur%next
          cur%next => new  
       end if
    else
       new%next => cur%next
       cur%next => new
    end if
  end subroutine extend

  type(cell) function pop(list)
    type(sortedlist) :: list
    type(cell), pointer :: c, p => null()

    c => list%head
    do while (associated(c%next))
       p => c
       c => c%next
    end do

    pop = c

    if (associated(p)) then
       nullify(p%next)
    end if
  end function pop
  
  recursive subroutine print_forward(n)
    type(cell), pointer, intent(in) :: n
    print *, n%d
    if (associated(n%next)) call print_forward(n%next)
  end subroutine print_forward

  recursive subroutine print_reverse(n)
    type(cell), pointer, intent(in) :: n
    if (associated(n%next)) call print_reverse(n%next)
    print *, n%d
  end subroutine print_reverse
  
end module slist
