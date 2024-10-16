! This module is to define the way that linked lists are created.

module slist
  implicit none

  ! If the list is the chain, these are the links. They contain the info
  ! and a reference (pointer) to the next link
  type :: cell
     ! The pointer is initialized as null, so it is defined but not
     ! associated to anything (might not be anything else yet)
     type(cell), pointer :: next => null()
     ! This is allocatable because we will introduce the data later
     character(:), allocatable :: d
  end type cell

  ! This is the list, and it contains the information for the first link
  ! and the length (number of links)
  type :: sortedlist
     integer :: length = 0
     type(cell), pointer :: head => null()
  end type sortedlist

contains

  ! Initialize the cell: allocate space for the string of unknown length
  subroutine init_cell(n, d)
    ! The pointer is supposed to leave the subroutine, so we indicate it.
    ! This references the cell we are creating. 
    type(cell), pointer, intent(out) :: n
    character(len=*) :: d
    allocate(n) ! Allocate space for this cell
    n%d = d     ! Assign the new string to the data "attr." of the cell 
  end subroutine init_cell
  
  ! This subroutine includes one new cell (string) to the list
  subroutine extend(list, s)
    type(sortedlist) :: list  ! Receive the list
    type(cell), pointer :: cur, new ! Two new pointers, to new cell and
    ! previous one

    ! New data
    character(len=*) :: s
    
    ! Create new cell
    call init_cell(new, s)

    ! point to the first element and cycle through to place new one
    cur => list%head
    do while (associated(cur%next)) ! this is a boolean
       ! Note that this only checks the first letter, not the full
       ! alphabetical order. Also, Fortran does not like (1), needs
       ! to be (1:1) for strings.
       if (new%d(1:1) > cur%next%d(1:1)) then
          cur => cur%next   ! Continue doing this until it is not ass.
       else    ! When it finally isn't, end the loop
          exit
       end if
    end do

    ! If it's associated with the head, update the head because it
    ! is a new first element
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

    ! Update list length
    list%length = list%length + 1
    
  end subroutine extend



  
  ! this function outputs a list element, therefore the "type(cell)"
  ! MIGHT BE BUGGED
  type(cell) function pop(list) 
    type(sortedlist) :: list
    type(cell), pointer :: c, p => null()

    c => list%head
    do while (associated(c%next)) ! when it reaches the end (ghost)
       ! it will not be associated anymore
       p => c
       c => c%next
    end do

    ! Here you say what the returned cell will be (same name as func.)
    pop = c

    if (associated(p)) then
       nullify(p%next) ! eliminate the memory of the popped cell
    end if

    ! Update list length
    list%length = list%length - 1
  end function pop



  ! This is recursive to reference the next "next" in a chain
  recursive subroutine print_forward(n)
    type(cell), pointer, intent(in) :: n
    print *, n%d
    if (associated(n%next)) call print_forward(n%next)
  end subroutine print_forward

  ! This is the EXACT same but reversing the two lines inside
  ! (save for the n declaration). 
  recursive subroutine print_reverse(n)
    type(cell), pointer, intent(in) :: n
    if (associated(n%next)) call print_reverse(n%next)
    print *, n%d
  end subroutine print_reverse
  
end module slist
