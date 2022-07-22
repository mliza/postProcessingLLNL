program mapping
    implicit none
    integer :: nx, ny, nz, i 
    integer, parameter :: n_max = 78036970 !change me, depends on nx, ny, nz
    integer, dimension(n_max, 4) :: mapping_matrix   
    character(len=53) mapping_vector 

    ! Change me, depends on nx, ny, nz 
    nx = 1439
    ny = 85
    nz = 638

    ! Output file 
    mapping_vector = '../../plate_data/data_10/temp_data/mappingVector.dat' 

    ! Write mapping for loop 
    do i = 1, n_max
        mapping_matrix(i,1) = mod(i-1, nx) 
        mapping_matrix(i,2) = int(mod(i-1, nx * ny) / nx)
        mapping_matrix(i,3) = int((i-1) / (nx * ny)) 
        mapping_matrix(i,4) = i - 1
    end do

    ! Printing testing 
    !do i = 1, n_max
    !    print *, mapping_matrix(i,:)
    !end do 
        

    ! Saving mapping 
    print *, 'Saving :', mapping_vector 
    open(unit=10, file=mapping_vector, & 
        form='unformatted', action='write', status='replace')
        do i = 1, n_max
            write(unit=10) mapping_matrix(i,:) 
        end do 
    close(unit=10) 

end program mapping
