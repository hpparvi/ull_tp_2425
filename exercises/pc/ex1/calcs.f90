module calcs
  use particle !only using the module particle because it's using the module geometry as well
  implicit none

  integer :: i, j
  real :: r2, r3

contains

  subroutine position(n, dt, p)
    integer, intent(in) :: n
    real, intent(in) :: dt
    type(particle3d), intent(inout) :: p(:)
    do i = 1, n
       p(i)%p = p(i)%p + p(i)%v *  dt
    end do
  end subroutine position

  subroutine velocity(n, dt, p)
    integer, intent(in) :: n
    real, intent(in) :: dt 
    type(particle3d), intent(inout) :: p(:)
    do i = 1, n
       p(i)%v = p(i)%v + p(i)%a * dt * 0.5
    end do
  end subroutine velocity

  subroutine set_acceleration(n, p)
    integer, intent(in) :: n
    type(particle3d), intent(inout) :: p(:)
    do i = 1, n
       p(i)%a = vector3d(0.0, 0.0, 0.0) !Setting the initial particles' accelerations to zero
    end do
  end subroutine set_acceleration

  subroutine acceleration(n, dt, p, rji)
    integer, intent(in) :: n
    real, intent(in) :: dt
    type(particle3d), intent(inout) :: p(:)
    type(vector3d), intent(inout) :: rji
    do i = 1, n
       do j = i + 1, n
          rji = p(j)%p - p(i)%p
          r2 = (norm(rji))**2
          r3 = r2 * sqrt(r2)
          p(i)%a = p(i)%a + p(j)%m * rji/r3
          p(j)%a = p(j)%a - p(i)%m * rji/r3
       end do
    end do
  end subroutine acceleration

end module calcs

