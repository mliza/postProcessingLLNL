program mapping
    implicit none
    integer :: nx, ny, nz, i 
    integer, parameter :: n_max = 77166738
    integer, dimension(n_max) :: x_mapping, y_mapping, z_mapping 
    character(len=30) mapping_vector 

    ! Define sizes of nx, ny, nz
    nx = 1359 
    ny = 89 
    nz = 638 
    ! Output file 
    mapping_vector = '../data_out/mappingVector.dat' 

    ! Write mapping for loop 
    do i = 1, (nx * ny * nz) 
        x_mapping(i) = mod(i-1, nx) 
        y_mapping(i) = int(mod(i-1, nx * ny) / nx)
        z_mapping(i) = int((i-1) / (nx * ny)) 
    end do 

    ! Saving mapping 
    print *, 'Saving :', mapping_vector 
    open(unit=10, file=mapping_vector, & 
        form='unformatted', action='write', status='replace')
    write(unit=10) x_mapping
    write(unit=10) y_mapping
    write(unit=10) z_mapping
    close(unit=10) 

end program mapping
