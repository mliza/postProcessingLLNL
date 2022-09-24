subroutine  scalar_reader(n_max, box_path, temp_path, var_in, time_step)
    implicit none 
    integer :: i_max, j_max, k_max, n_blocks
    double precision :: mach_number, ang_of_attack, reynolds_number, iteration
    integer, intent(in) :: n_max
    double precision, dimension(n_max) :: scalar_variable 
    character(len=*), intent(in) :: box_path, temp_path, var_in, time_step
    character(len=80) file_in, file_out 

    ! Input output grid 
    file_in  = box_path//'/'//var_in//'.'//time_step//'.q'
    file_out = temp_path//'/'//time_step//'_'//var_in//'.dat'

    print *, 'Loading: ', file_in  ! Print statement  
    open(unit=7, file=file_in, & 
         form='unformatted', access ='stream')
    read(unit=7) n_blocks  
    read(unit=7) i_max, j_max, k_max 
    read(unit=7) mach_number, ang_of_attack, reynolds_number, iteration
    read(unit=7) scalar_variable  
    close(unit=7) 

    print *, 'Saving: ', file_out ! Print statement  
    open(unit=8, file=file_out, &
         form='unformatted', action='write', status='replace')
    write(unit=8) scalar_variable 
    close(unit=8)

end subroutine  scalar_reader 
