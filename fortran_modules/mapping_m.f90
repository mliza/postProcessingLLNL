
! Creates a temp file with a plot3D mapping 
subroutine mapping(nx, ny, nz, mapping_path)
    implicit none
    integer :: i
    integer, intent(in) :: nx, ny, nz
    integer, dimension(nx * ny * nz) :: x_mapping, y_mapping, z_mapping 
    character(len=*), intent(in) :: mapping_path
    character(len=53) mapping_vector

    ! Output file 
    mapping_vector = mapping_path//'/mappingVector.dat'

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

end subroutine mapping

