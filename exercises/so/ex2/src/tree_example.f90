module tree
  use geometry
  use particles
  
  implicit none

  type :: npointer
     type(node), pointer :: node => null()
  end type npointer

  type :: node
     !type(particle), pointer          :: part                  ! particle
     real                             :: m_tot = 0.0           ! total mass
     type(point3d)                    :: c_o_m                 ! Center of mass
     integer                          :: level = 0, cid = -1   ! Depth in the tree, id of each child
     type(node), pointer              :: parent => null()
     type(npointer), dimension(2,2,2) :: child                 ! All 8 children nodes
  end type node
  
end module tree


program treetest
  use tree
  implicit none
  type(node), pointer :: root => null()

  call init_node(root, null(), 3.0, 0)
  call set_node(root, [1], 4.0)
  call set_node(root, [1,1], 2.0)
  call set_node(root, [1,2], 4.0)
  call set_node(root, [1,3], 3.0)
  call set_node(root, [2,2,4], 3.0)
  call set_node(root, [3], 2.0)
  call set_node(root, [3,1], 1.0)
  
  print *, "Print forward"
  call print_tree_forward(root)
  print *
  print *, "Print reverse"
  call print_tree_reverse(root)
  print *
  print *, "Cumulative sum"
  call csum(root)
  
contains


  subroutine init_node(n, p, v, cid)
    type(node), pointer :: p, n
    real, intent(in) :: v
    integer :: i, cid
    
    if (.not. associated(n)) then
       allocate(n)
       n%v = v
       n%cid = cid
       n%parent => p
       if (associated(p)) then
          n%level = p%level + 1
       end if
    end if
  end subroutine init_node

  
  subroutine add_node(n, cid, v)
    type(node), pointer :: n
    integer :: cid
    real :: v
    call init_node(n%child(cid)%node, n, v, cid)
  end subroutine add_node

  
  subroutine set_node(n, ids, v)
    type(node), pointer :: n, t
    integer, dimension(:) :: ids
    real :: v
    integer :: nlevels, i

    nlevels = size(ids)
    t => n
    do i = 1, nlevels
       if (.not. associated(t%child(ids(i))%node)) then
          call add_node(t, ids(i), 0.0)
       end if
       t => t%child(ids(i))%node
    end do
    t%v = v
  end subroutine set_node
  
  
  subroutine print_node(n)
    type(node), pointer :: n
    call print_index(n)
    print '(a,i2,2x,f4.2,2x,f6.2)', char(9), n%level, n%v, n%vsum
  end subroutine print_node

  
  recursive subroutine print_index(n)
    type(node), pointer :: n, t
    if (associated(n%parent)) then
       call print_index(n%parent)
       write(*, fmt="(a)", advance="no") "."
    end if
    write(*, fmt="(i0)", advance="no") n%cid
  end subroutine print_index

  
 recursive subroutine print_tree_forward(n)
    type(node), pointer :: n, t
    integer :: i

    call print_node(n)
    do i = 1, 4
       if (associated(n%child(i)%node)) then
          call print_tree_forward(n%child(i)%node)
       end if
    end do
  end subroutine print_tree_forward

  
  recursive subroutine print_tree_reverse(n)
    type(node), pointer :: n
    integer :: i
    
    do i = 1, 4
       if (associated(n%child(i)%node)) then
          call print_tree_reverse(n%child(i)%node)
       end if
    end do
    call print_node(n)
  end subroutine print_tree_reverse

  
  recursive subroutine csum(n)
    type(node), pointer :: n
    integer :: i

    n%vsum = n%v
    do i = 1, 4
       if (associated(n%child(i)%node)) then
          call csum(n%child(i)%node)
          n%vsum = n%vsum + n%child(i)%node%vsum
       end if
    end do
    call print_node(n)
  end subroutine csum
end program treetest
