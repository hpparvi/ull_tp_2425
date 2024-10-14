module geometry
  implicit none

  type point3d
     real :: xx, yy, zz
  end type point3d
  
  type vector3d
     real :: xx, yy, zz
  end type vector3d

  interface operator(+)
     module procedure sumvv, sumvp, sumpv
  end interface operator(+)

  interface operator(-)
     module procedure subvv, subvp, subpv
  end interface operator(-)

  interface operator(*)
     module procedure mulrv, mulvr
  end interface operator(*)

  interface operator(/)
     module procedure divvr
  end interface operator(/)
  

contains
  
  pure type(vector3d) function sumvv(vector1, vector2) !sums v1 and v2
    type(vector3d), intent(in) :: vector1, vector2

    sumvv%xx = vector1%xx + vector2%xx
    sumvv%yy = vector1%yy + vector2%yy
    sumvv%zz = vector1%zz + vector2%zz
  end function sumvv

  pure type(vector3d) function subvv(vector1, vector2) !subtracts v2 from v1
    type(vector3d), intent(in) :: vector1, vector2

    subvv%xx = vector1%xx - vector2%xx
    subvv%yy = vector1%yy - vector2%yy
    subvv%zz = vector1%zz - vector2%zz
  end function subvv

  pure type(point3d) function sumpv(point1, vector1) !sums p1 and v1, returns point
    type(vector3d), intent(in) :: vector1
    type(point3d), intent(in) :: point1
    
    sumpv%xx = point1%xx + vector1%xx
    sumpv%yy = point1%yy + vector1%yy
    sumpv%zz = point1%zz + vector1%zz
  end function sumpv

  pure type(point3d) function sumvp(vector1, point1) !sums p1 and v1, returns point
    type(vector3d), intent(in) :: vector1
    type(point3d), intent(in) :: point1
    
    sumvp = sumpv(point1, vector1)
  end function sumvp

  
  pure type(point3d) function subvp(vector1, point1) !subs p1 from v1, returns point
    type(vector3d), intent(in) :: vector1
    type(point3d), intent(in) :: point1
    
    subvp%xx = vector1%xx - point1%xx
    subvp%yy = vector1%yy - point1%yy
    subvp%zz = vector1%zz - point1%zz
  end function subvp

  pure type(point3d) function subpv(point1, vector1) !subs v1 from p1, returns point
    type(vector3d), intent(in) :: vector1
    type(point3d), intent(in)  :: point1
    
    subpv%xx = point1%xx - vector1%xx
    subpv%yy = point1%yy - vector1%yy
    subpv%zz = point1%zz - vector1%zz
  end function subpv
    

  pure type(vector3d) function mulrv(scalar1, vector1) !does a*v1
    type(vector3d), intent(in) :: vector1
    real, intent(in)           :: scalar1

    mulrv%xx = vector1%xx * scalar1
    mulrv%yy = vector1%yy * scalar1
    mulrv%zz = vector1%zz * scalar1    
  end function mulrv

  pure type(vector3d) function mulvr(vector1, scalar1) !v1*a
    type(vector3d), intent(in) :: vector1
    real, intent(in) :: scalar1

    mulvr = mulrv(scalar1, vector1)
  end function mulvr

  pure type(vector3d) function divvr(vector1, scalar1) !v1/a
    type(vector3d), intent(in) :: vector1
    real, intent(in)           :: scalar1

    divvr%xx = vector1%xx / scalar1
    divvr%yy = vector1%yy / scalar1
    divvr%zz = vector1%zz / scalar1
  end function divvr

  pure real function dotprod(vector1, vector2) !(v1, v2)
    type(vector3d), intent(in) :: vector1, vector2

    dotprod = vector1%xx * vector2%xx + &
         &vector1%yy * vector2%yy + &
         &vector1%zz * vector2%zz
  end function dotprod

  pure real function distance(point1, point2)
    type(point3d), intent(in) :: point1, point2

    distance = sqrt((point1%xx - point2%xx)**2 + &
         &(point1%yy - point2%yy)**2 + &
         &(point1%zz - point2%zz)**2)
  end function distance

  pure real function vecnorm(vector1) !norm(v1) for later use
    type(vector3d), intent(in) :: vector1

    vecnorm = sqrt(vector1%xx**2 + vector1%yy**2 + vector1%zz**2)
  end function vecnorm

  pure real function angle(vector1, vector2) !inner angle between v1, v2, in radians.
    type(vector3d), intent(in) :: vector1, vector2

    angle = acos(dotprod(vector1, vector2)/(vecnorm(vector1)*vecnorm(vector2)))
  end function angle

  pure type(vector3d) function normalize(vector1) !v1/norm(v1)
    type(vector3d), intent(in) :: vector1

    normalize = vector1/vecnorm(vector1)
  end function normalize

  pure type(vector3d) function cross_product(vector1, vector2)
    type(vector3d), intent(in) :: vector1, vector2

    cross_product%xx = vector1%xx * vector2%zz - vector1%zz * vector2%yy
    cross_product%yy = vector1%zz * vector2%xx - vector1%xx * vector2%zz
    cross_product%zz = vector1%xx * vector2%yy - vector1%yy * vector2%xx
  end function cross_product

  pure type(vector3d) function orthv(vector1, vector2)
    type(vector3d), intent(in) :: vector1, vector2
    
    orthv = cross_product(vector1, vector2)
  end function orthv
  
end module geometry
