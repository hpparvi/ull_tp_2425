program e1v2
  implicit none
  integer :: ix, iy
  integer :: sax, say, sbx, sby
  real :: A(3,4), B(4,3)
  real, allocatable :: C(:,:)

  say = size(A, 1)
  sax = size(A, 2)
  sby = size(B, 1)
  sbx = size(B, 2)
  
  A = reshape([3, 2, 4, 1,  &
       &       2, 4, 2, 2,  &
       &       1, 2, 3, 7], &
       & [say, sax], order=[2, 1])

  B = reshape([3, 2, 4,  &
       &       2, 4, 2,  &
       &       1, 2, 3,  &
       &       0, 2, 1], &
       & [sby, sbx], order=[2, 1])

  allocate(C(say, sbx))
  
  do iy = 1, say
     do ix = 1, sbx
        C(iy, ix) = sum(A(iy, :) *  B(:, ix))
     end do
  end do

  do iy = 1, size(C, 1)
     print '(*(F6.1,2X))', C(iy, :)
  end do
  
end program e1v2
  
