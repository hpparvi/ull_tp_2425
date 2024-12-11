module definitions     !This module defines key data structures for N-body simulations.

  use iso_fortran_env  !This module ensures all variables are defined as 64-bit.
  use omp_lib          !This module is used for parallelization with OpenMP
  use geometry         !This module defines 3D vector and point operations for vector3d and point3d types.
  implicit none

  type :: RANGE
     real(real64), dimension(3) :: min, max  !Minimum and maximum bounds in 3D.
  end type RANGE

  
  type :: CPtr
     type(CELL), pointer :: ptr  !Pointer to a CELL object.
  end type CPtr


  type :: particle3d
     type(point3d) :: p          !Position of the particle.
     type(vector3d) :: v         !Velocity of the particle.
     real(real64) :: m           !Mass of the particle.
  end type particle3d


  type :: CELL
     type(RANGE) :: range        !The spatial range of the cell.
     type(particle3d) :: p       !The particle contained in the cell, if any.
     integer(int64) :: pos       !Position index of the cell.
     integer(int64) :: type      !Type of cell: 0 = no particle, 1 = contains a particle, 2 = contains a cluster of particles.
     real(real64) :: mass        !Total mass of the cell (particle or cluster).
     type(point3d) :: c_o_m      !Center of mass of the cell.
     type(CPtr), dimension(2,2,2) :: subcell  ! Pointers to the eight subcells (octree structure).
  end type CELL

end module definitions
