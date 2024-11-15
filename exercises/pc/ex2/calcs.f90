module calcs
  use geometry !Importing the geometry module
  use particle !Importing the particle module
  use barnes   !Importing the barnes module where the Barnes-Hut algorithm is defined
  use iso_fortran_env !Importing iso_fortran_env to specify the number of bits of the variables
  implicit none

contains

  !Subroutine to update the position of particles over time
  subroutine position(n, dt, p)
    integer(int64), intent(in) :: n !Number of particles
    real(real64), intent(in) :: dt  !Time step
    type(particle3d), intent(inout) :: p(:) !Particles
    do i = 1, n
       p(i)%p = p(i)%p + p(i)%v * dt !Update position using velocity
    end do
  end subroutine position

  !Subroutine to update the velocity of particles over time
  subroutine velocity(n, dt, p)
    integer(int64), intent(in) :: n !Number of particles
    real(real64), intent(in) :: dt  !Time step
    type(particle3d), intent(inout) :: p(:) !Particles
    do i = 1, n
       p(i)%v = p(i)%v + p(i)%a * (dt * 0.5) !Update velocity using acceleration
    end do
  end subroutine velocity

  !Subroutine to set acceleration of all particles to zero
  subroutine set_acceleration(n, p)
    integer(int64), intent(in) :: n !Number of particles
    type(particle3d), intent(inout) :: p(:) !Particles
    do i = 1, n
       p(i)%a = vector3d(0.0, 0.0, 0.0) !Set initial acceleration of each particle to zero
    end do
  end subroutine set_acceleration

  !Subroutine to calculate gravitational acceleration on each particle due to others
  subroutine acceleration(n, dt, p, rji)
    integer(int64), intent(in) :: n !Number of particles
    real(real64), intent(in) :: dt  !Time step
    type(particle3d), intent(inout) :: p(:) !Particles
    type(vector3d), intent(inout) :: rji !Vector from one particle to another

    !Calculate acceleration between particles
    do i = 1, n
       do j = i + 1, n
          rji = p(j)%p - p(i)%p !Determine vector between particle j and particle i
          r2 = (norm(rji))**2 !Square of the vector
          r3 = r2 * sqrt(r2)  !Cube of the vector
          p(i)%a = p(i)%a + p(j)%m * rji / r3 !Update acceleration on particle i
          p(j)%a = p(j)%a - p(i)%m * rji / r3 !Update acceleration on particle j
       end do
    end do
  end subroutine acceleration

end module calcs

