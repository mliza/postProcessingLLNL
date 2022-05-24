subroutine vector_reader(n_max, box_path, temp_path, var_in, time_step)
    implicit none 
    integer :: i_max, j_max, k_max, n_blocks 
    double precision :: mach_number, ang_of_attack, reynolds_number, iteration
    integer, intent(in) :: n_max 
    double precision, dimension(n_max) :: x_variable, y_variable, z_variable 
    character(len=*), intent(in) :: box_path, temp_path, var_in, time_step
    character(len=50) file_in, file_out_x, file_out_y, file_out_z 

    ! Input output grid 
    !file_in    = '../../plate_data/data_12/smallBOX_BIN/U.0848001.q' 
    !file_out_x = '../../plate_data/data_12/temp_data/Ux.dat'
    !file_out_y = '../../plate_data/data_12/temp_data/Uy.dat'
    !file_out_z = '../../plate_data/data_12/temp_data/Uz.dat'
    file_in    = box_path//'/'//var_in//'.'//time_step//'.q'
    file_out_x = temp_path//'/'//var_in//'x.dat'
    file_out_y = temp_path//'/'//var_in//'y.dat'
    file_out_z = temp_path//'/'//var_in//'z.dat'

    ! Using cases 
    print *, 'Loading: ', file_in  ! Print statement  
    open(unit=7, file=file_in, & 
         form='unformatted', access ='stream')

    ! For a 3D vector variable 
    read(unit=7) n_blocks  
    read(unit=7) i_max, j_max, k_max 
    read(unit=7) mach_number, ang_of_attack, reynolds_number, iteration
    read(unit=7) x_variable
    read(unit=7) y_variable
    read(unit=7) z_variable
    close(unit=7) 
    ! Writing the vector in 3 different files 
    ! X-vector  
    print *, 'Saving: ', file_out_x  ! Print statement  
    open(unit=10, file=file_out_x, &
         form='unformatted', action='write', status='replace')
    write(unit=10) x_variable 
    close(unit=10)
     ! Y-vector  
    print *, 'Saving: ', file_out_y  ! Print statement  
    open(unit=11, file=file_out_y, &
         form='unformatted', action='write', status='replace')
    write(unit=11) y_variable 
    close(unit=11)
    ! Z-vector  
    print *, 'Saving: ', file_out_z  ! Print statement  
    open(unit=12, file=file_out_z, &
         form='unformatted', action='write', status='replace')
    write(unit=12) z_variable 
    close(unit=12)

end subroutine vector_reader 
