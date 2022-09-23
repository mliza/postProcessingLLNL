subroutine mapping(nx, ny, nz, mapping_path) 
    implicit none
    integer, intent(in) :: nx, ny, nz
    character(len=*), intent(in) :: mapping_path 
    integer :: i, n_max 
    integer, dimension(nx * ny * nz, 4) :: mapping_matrix   
    character(len=55) mapping_vector 
    n_max = nx * ny * nz

    ! Output file 
    mapping_vector = mapping_path//'/mappingVector.dat' 

    ! Write mapping for loop 
    do i = 1, n_max
        mapping_matrix(i,1) = mod(i-1, nx) 
        mapping_matrix(i,2) = int(mod(i-1, nx * ny) / nx)
        mapping_matrix(i,3) = int((i-1) / (nx * ny)) 
        mapping_matrix(i,4) = i - 1
    end do

    ! Saving mapping 
    print *, 'Saving :', mapping_vector 
    open(unit=10, file=mapping_vector, & 
        form='unformatted', action='write', status='replace')
        do i = 1, n_max
            write(unit=10) mapping_matrix(i,:) 
        end do 
    close(unit=10) 

end subroutine mapping
