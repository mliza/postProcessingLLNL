program scalar_reader
    implicit none 
    integer :: i_max, j_max, k_max, n_blocks
    double precision :: mach_number, angle_of_attack, reynolds_number, time_step
    integer, parameter :: n_max = 30 !change me, depends on nx, ny, nz   
    double precision, dimension(n_max) :: scalar_variable 
    character(len=50) file_in, file_out 

    ! Input output grid 
    file_in    = '../../plate_data/data_12/smallBOX_BIN/T.0848001.q'
    file_out   = '../../plate_data/data_12/temp_data/T.dat'

    ! Loading data 
    print *, 'Loading: ', file_in  
    open(unit=7, file=file_in, & 
         form='unformatted', access ='stream')
    read(unit=7) n_blocks  
    read(unit=7) i_max, j_max, k_max 
    read(unit=7) mach_number, angle_of_attack, reynolds_number, time_step
    read(unit=7) scalar_variable  
    close(unit=7) 
    ! Saving data 
    print *, 'Saving: ', file_out 
    open(unit=8, file=file_out, &
         form='unformatted', action='write', status='replace')
    write(unit=8) scalar_variable 
    close(unit=8)

end program scalar_reader 
