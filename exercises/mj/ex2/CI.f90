program particulas
    implicit none
    integer :: i, n
    integer :: values(1:8), k
    integer, dimension(:), allocatable :: seed
    real :: mass, rx, ry, rz


    ! Open the file for writing
    open(unit=10, file='IC.txt', status='replace')

    ! Ask the user for the number of bodies
    print*, "Number of bodies?"
    read*, n
    mass = 1.0 / n

    ! Generate and write particle positions to the file
    do i = 1, n
        call random_number(rx)
        do
            call random_number(ry)
            if ((rx**2 + ry**2) <= 1) exit
        end do
        do
            call random_number(rz)
            if ((rx**2 + ry**2 + rz**2) <= 1) exit
        end do
        write(10, '(F6.3, 3F11.8, 3I2)') mass, rx, ry, rz, 0, 0, 0
    end do

    ! Close the file
    close(unit=10)

    print*, "Results saved to IC"
end program particulas
