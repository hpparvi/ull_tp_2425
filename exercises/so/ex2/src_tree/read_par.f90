module parameter_reader
    implicit none
    ! Declare parameter variables
    real :: dt, dt_out, t_end, epsilon, theta, radius
    character(len=100) :: input_file, output_file
    logical :: create_bodies
    integer :: N_bodies

contains

    subroutine read_parameters(filename)
        implicit none
        character(len=*), intent(in) :: filename
        character(len=256) :: line
        character(len=30) :: key, key_short
        character(len=20) :: dummy_char, dummy_two, dummy_one
        integer :: ios

        open(unit=10, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error: Could not open file ", filename
            stop
        end if

        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! End of file
            line = adjustl(line)  ! Remove leading spaces

            if (line(1:1) == "#" .or. trim(line) == "") cycle  ! Skip comments/empty lines

            ! Parse the line
            read(line, '(A)', iostat=ios) dummy_char
            read(dummy_char, *) key_short
            key = trim(key_short)
            line = trim(adjustl(line))
            select case (trim(adjustl(key)))
                case ("dt")
                    read(line, *) dummy_one, dummy_two, dt
                case ("dt_out")
                    read(line, *) dummy_one, dummy_two, dt_out
                case ("t_end")
                    read(line, *) dummy_one, dummy_two, t_end
                case ("input_file")
                   read(line, *) dummy_one, dummy_two, input_file
                case ("output_file")
                   read(line, *) dummy_one, dummy_two, output_file
                   print*, output_file
                case ("create_bodies")
                    read(line, *) dummy_one, dummy_two, create_bodies
                case ("N_bodies")
                   read(line, *) dummy_one, dummy_two, N_bodies
                case ("radius")
                   read(line, *) dummy_one, dummy_two, radius
                case ("epsilon")
                   read(line, *) dummy_one, dummy_two, epsilon
                case ("theta")
                   read(line, *) dummy_one, dummy_two, theta
                case default
                    print *, "Warning: Unrecognized parameter:", key
                 end select
        end do

        close(10)
    end subroutine read_parameters
end module parameter_reader
