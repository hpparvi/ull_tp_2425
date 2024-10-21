program nbody_program
   use geometry
   use particle
   implicit none

   type ( particle3d ) :: b1

   b1%v = vector3d(1, 2, 3)
   b1%p = point3d(1, 2, 3)
   b1%m = 10

   print *, b1
   
end program nbody_program

