program scalar_reader
    implicit none 
    integer :: i_max, j_max, k_max, n_blocks
    double precision :: mach_number, angle_of_attack, reynolds_number, time_step
    integer, parameter :: n_max = 77166738 
    double precision, dimension(n_max) :: scalar_variable 
    character(len=30) file_in, file_out 

    ! Input output grid 
    file_in    = '../data/GRADRHOMAG.0838100.q' 
    file_out   = '../dataOut/GRADRHOMAG.dat' 

    print *, 'Loading: ', file_in  ! Print statement  
    open(unit=7, file=file_in, & 
         form='unformatted', access ='stream')
    read(unit=7) n_blocks  
    read(unit=7) i_max, j_max, k_max 
    read(unit=7) mach_number, angle_of_attack, reynolds_number, time_step
    read(unit=7) scalar_variable  
    close(unit=7) 
    print *, 'Saving: ', file_out ! Print statement  
    open(unit=8, file=file_out, &
         form='unformatted', action='write', status='replace')
    write(unit=8) scalar_variable 
    close(unit=8)

end program scalar_reader 
