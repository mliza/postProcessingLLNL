program vector_reader
    implicit none 
    integer :: i_max, j_max, k_max, n_blocks 
    double precision :: mach_number, angle_of_attack, reynolds_number, time_step
    integer, parameter :: n_max = 77166738 
    double precision, dimension(n_max) :: x_variable, y_variable, z_variable 
    character(len=20) file_in, file_out_x, file_out_y, file_out_z 

    ! Input output grid 
    !file_in    = '../data/U.0838100.q' 
    file_in    = '../data/U.0838100.q' 
    file_out_x = '../dataOut/Ux.dat'
    file_out_y = '../dataOut/Uy.dat'
    file_out_z = '../dataOut/Uz.dat'

    ! Using cases 
    print *, 'Loading: ', file_in  ! Print statement  
    open(unit=7, file=file_in, & 
         form='unformatted', access ='stream')

    ! For a 3D vector variable 
    read(unit=7) n_blocks  
    read(unit=7) i_max, j_max, k_max 
    read(unit=7) mach_number, angle_of_attack, reynolds_number, time_step
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

end program vector_reader 

