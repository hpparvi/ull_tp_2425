PROGRAM leapfrog

  use geometry
  use particles

  IMPLICIT NONE

  INTEGER :: i,j,k, i_t, ios
  INTEGER :: n, n_t
  INTEGER :: num_args       ! to count the number of given arguments
  character(len=100) :: filename     ! to import the stars data
  character(len=256) :: dummy_line ! to ignore first line of file when reading data
  character(len=100)  :: output_file = "output.txt"
  real :: dt, t_end, t, dt_out, t_out
  real :: rs, r2, r3
  
  type(particle), dimension(:), allocatable :: partics
  type(vector3d)                            :: rji_v
  real                                      :: r3_2
  
  real, DIMENSION(:), ALLOCATABLE :: m
  real, DIMENSION(:,:), ALLOCATABLE :: r,v,a
  real, DIMENSION(3) :: rji

  ! Get the number of command-line arguments
    num_args = command_argument_count()

    if (num_args < 1) then
        print *, "Error: No file name provided. Usage: ./program <filename>"
        stop
    end if

    ! Get the input file name from the command-line arguments
    call get_command_argument(1, filename)

    ! Remove any trailing spaces from the filename
    filename = trim(adjustl(filename))

    ! Open the file
    open(unit=10, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) then
        print *, "Error: Could not open file ", filename
        stop
    end if

    ! Count the number of stars
    n = 0
    do
        read(10, '(A)', iostat=ios) ! reading each line as a whole string
        if (ios /= 0) exit
        n = n + 1
     end do
     n = n-1 ! to remove the first line with time parameters
    rewind(10) !so that pointer is at the beginnin for reading and importing data later

    ! Allocate arrays based on the number of stars
! Particle method for calculations
 allocate(partics(n))
! No particle method (original one) used for checking. Results must match
 ALLOCATE(m(n))
 ALLOCATE(r(n,3))
 ALLOCATE(v(n,3))
 ALLOCATE(a(n,3))
    
 ! Read the data into arrays (no particle method)
 read(10, '(A)', iostat=ios) dummy_line
 if (ios /= 0) then
        print *, "Error: Could not read the first line"
        stop
     end if
     
    do i = 1, n
        read(10, *, iostat=ios) m(i), r(i,1), r(i,2), r(i,3), v(i,1), v(i,2), v(i,3)
        if (ios /= 0) then
            print *, "Error reading file, stopped at star:", i
            exit
        end if
     end do
     print *, "Values saved to stars (no particle method)"
rewind(10)
     
! Repeat reading for particle variables
read(10, *, iostat = ios) dt, dt_out, t_end
 if (ios /= 0) then
        print *, "Error: Could not read the first line"
        stop
     end if
     
    do i = 1, n
        read(10, *, iostat=ios) partics(i)%m, partics(i)%p%x, partics(i)%p%y, partics(i)%p%z, partics(i)%v
        if (ios /= 0) then
            print *, "Error reading file, stopped at star:", i
            exit
        end if
     end do
     print *, "Values saved to particle type variables"

    ! Close the file
    close(10)

  n_t = t_end/dt

! initialize accelerations
  do i = 1,n
     partics(i)%a = vector3d(0,0,0)
  end do
! initialize values
  DO i = 1,n
     DO j = i+1,n
        ! Particles method
        rji_v = vecpp(partics(i)%p, partics(j)%p)
        r3_2 = norm(rji_v)**3
        partics(i)%a = partics(i)%a + (partics(j)%m * rji_v / r3_2)
        partics(j)%a = partics(j)%a - (partics(i)%m * rji_v / r3_2)
        ! Original method
        rji = r(j,:) - r(i,:)
        r2 = SUM(rji**2)
        r3 = r2 * SQRT(r2)
        a(i,:) = a(i,:) + m(j) * rji / r3
        a(j,:) = a(j,:) - m(i) * rji / r3
     END DO
  END DO
  
  t_out = 0.0
  
  ! Open the output file
  open(unit=20, file=output_file, status="replace", action="write", iostat=ios)
     if (ios /= 0) then
        print *, "Error: Could not open file ", filename
        stop
     end if
     
! start all calculations
  DO i_t = 0, n_t
     t = i_t * dt
     ! Particles method
     partics(:)%v = partics(:)%v + (partics(:)%a * (dt/2))
     partics(:)%p = partics(:)%p + (partics(:)%v * dt)
     partics(:)%a = vector3d(0,0,0)
     ! Original method
     v = v + a * dt/2
     r = r + v * dt
     a = 0.0
     
     DO i = 1,n
        DO j = i+1,n
           ! Particles method
           rji_v = vecpp(partics(i)%p, partics(j)%p)
           r3_2 = norm(rji_v)**3
           partics(i)%a = partics(i)%a + (partics(j)%m * rji_v / r3_2)
           partics(j)%a = partics(j)%a - (partics(i)%m * rji_v / r3_2)
           ! Original method
           rji = r(j,:) - r(i,:)
           r2 = SUM(rji**2)
           r3 = r2 * SQRT(r2)
           a(i,:) = a(i,:) + m(j) * rji / r3
           a(j,:) = a(j,:) - m(i) * rji / r3
        END DO
        ! Store all values of positions:
     write(20, '(I5, F12.2, 3F20.15)') i_t, t, partics(i)%p
     write(20, '(I5, F12.2, 3F20.15)') i_t, t, r(i,:)
     END DO
     !Particles method
     partics(:)%v = partics(:)%v + (partics(:)%a * (dt/2))
     ! Original method
     v = v + a * dt/2
     
     t_out = t_out + dt
     IF (t_out >= dt_out) THEN
        print*, "Time = ", t, "", "Iteration = ", i_t
        DO i = 1,n
           print*, "Particle ", i
           print*, partics(i)%p
           print*, r(i,:)
        END DO
        t_out = 0.0
     END IF
     
  END DO

  close(20)
  
END PROGRAM leapfrog
