module calcs
  use geometry !Importing the geometry module
  use particle !Importing the particle module
  use barnes   !Importing the barnes module where the Barnes-Hut algorithm is defined
  use iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  implicit none

contains

  !Subroutine to update the position of particles over time
  subroutine position(dt, p)
    real(real64), intent(in) :: dt  !Time step
    type(particle3d), intent(inout) :: p(:) !Particles
    INTEGER(INT64) :: i !Loop indexing variable
    do i = 1, size(p)
       p(i)%p = p(i)%p + p(i)%v * dt !Update position using velocity
    end do
  end subroutine position

  !Subroutine to update the velocity of particles over time
  subroutine velocity(dt, p)
    real(real64), intent(in) :: dt  !Time step
    type(particle3d), intent(inout) :: p(:) !Particles
    INTEGER(INT64) :: i !Loop indexing variable
    do i = 1, size(p)
       p(i)%v = p(i)%v + p(i)%a * (dt * 0.5) !Update velocity using acceleration
    end do
  end subroutine velocity

  !Subroutine to set acceleration of all particles to zero
  subroutine set_acceleration(p)
    type(particle3d), intent(inout) :: p(:) !Particles
    INTEGER(INT64) :: i !Loop indexing variable
    do i = 1, size(p)
       p(i)%a = vector3d(0.0, 0.0, 0.0) !Set initial acceleration of each particle to zero
    end do
  end subroutine set_acceleration

end module calcs

