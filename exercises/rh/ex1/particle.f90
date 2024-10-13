
module particles

  use geometry
  
  type particle
    type(point3d) :: p ! particle position
    type(vector3d) :: v ! velocity of the particle
    real :: m ! mass of the particle
  end type particle
  
end module particles
