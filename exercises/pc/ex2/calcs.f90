module calcs
  use geometry !Importing the geometry module
  use particle !Importing the particle module
  use barnes   !Importing the barnes module where the Barnes-Hut algorithm is defined
  use iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  !$use omp_lib  !Importing the omp library to use OpenMP
  implicit none

contains

  !Subroutine to update the position of particles over time
  subroutine position(n, dt, p)
    integer(int64), intent(in) :: n !Number of particles
    real(real64), intent(in) :: dt  !Time step
    type(particle3d), intent(inout) :: p(:) !Particles
    INTEGER(INT64) :: i !Loop indexing variable
    !$omp do
    do i = 1, n
       p(i)%p = p(i)%p + p(i)%v * dt !Update position using velocity
    end do
    !$omp end do
  end subroutine position

  !Subroutine to update the velocity of particles over time
  subroutine velocity(n, dt, p)
    integer(int64), intent(in) :: n !Number of particles
    real(real64), intent(in) :: dt  !Time step
    type(particle3d), intent(inout) :: p(:) !Particles
    INTEGER(INT64) :: i !Loop indexing variable
    !$omp do
    do i = 1, n
       p(i)%v = p(i)%v + p(i)%a * (dt * 0.5) !Update velocity using acceleration
    end do
    !$omp end do
  end subroutine velocity

  !Subroutine to set acceleration of all particles to zero
  subroutine set_acceleration(n, p)
    integer(int64), intent(in) :: n !Number of particles
    type(particle3d), intent(inout) :: p(:) !Particles
    INTEGER(INT64) :: i !Loop indexing variable
    !$omp do
    do i = 1, n
       p(i)%a = vector3d(0.0, 0.0, 0.0) !Set initial acceleration of each particle to zero
    end do
    !$omp end do
  end subroutine set_acceleration

end module calcs

