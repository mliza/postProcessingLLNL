! Reads the Grid 
subroutine grid_reader(n_max, box_path, temp_path, var_in)
    implicit none 
    integer :: i_max, j_max, k_max, n_blocks
    double precision :: mach_number, angle_of_attack, reynolds_number, time_step
    integer, intent(in):: n_max
    double precision, dimension(n_max) :: x_variable, y_variable, z_variable 
    character(len=*), intent(in) :: box_path, temp_path, var_in
    character(len=53) grid_in, grid_out_x, grid_out_y, grid_out_z 

    ! Input and output files  
    grid_in     = box_path//'/'//var_in//'.xyz'
    grid_out_x  = temp_path//'/X.dat'
    grid_out_y  = temp_path//'/Y.dat'
    grid_out_z  = temp_path//'/Z.dat'
    !grid_in    = '../../plate_data/data_12/smallBOX_BIN/T.xyz'
    !grid_out_x = '../../plate_data/data_12/temp_data/X.dat'
    !grid_out_y = '../../plate_data/data_12/temp_data/Y.dat'
    !grid_out_z = '../../plate_data/data_12/temp_data/Z.dat'

    ! Opening and reading file 
    print *, 'Loading: ', grid_in  ! Print statement  
    open(unit=7, file=grid_in, & 
         form='unformatted', access ='stream')
    read(unit=7) n_blocks  
    read(unit=7) i_max, j_max, k_max 
    read(unit=7) x_variable
    read(unit=7) y_variable
    read(unit=7) z_variable

    ! X-vector  
    print *, 'Saving: ', grid_out_x  ! Print statement  
    open(unit=10, file=grid_out_x, &
         form='unformatted', action='write', status='replace')
    write(unit=10) x_variable 
    close(unit=10)

    ! Y-vector  
    print *, 'Saving: ', grid_out_y  ! Print statement  
    open(unit=11, file=grid_out_y, &
         form='unformatted', action='write', status='replace')
    write(unit=11) y_variable 
    close(unit=11)

    ! Z-vector  
    print *, 'Saving: ', grid_out_z  ! Print statement  
    open(unit=12, file=grid_out_z, &
         form='unformatted', action='write', status='replace')
    write(unit=12) z_variable 
    close(unit=12)

end subroutine grid_reader 
