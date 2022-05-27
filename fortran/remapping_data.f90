program remapping_data
    implicit none
    integer :: nx, ny, nz, i, j, k, n 
    integer, parameter :: n_max = 30 !change me, depends on nx, ny, nz
    integer, dimension(n_max, 4) :: mapping_matrix   
    double precision, dimension(n_max) :: data_array 
    double precision, dimension(nx, ny, nz) :: mdim_array 
    
    character(len=53) mapping_in, data_in, data_out 

    ! Change me, depends on nx, ny, nz 
    nx = 5
    ny = 3 
    nz = 2

    ! Output file 
    mapping_in = '../../plate_data/data_12/temp_data/mappingVector.dat' 
    data_in    = '../../plate_data/data_12/temp_data/X.dat' 
    data_out   = '../../plate_data/data_12/temp_data/X_3D.dat' 

    ! Loading mapping matrix 
    print *, 'Loading :', mapping_in 
    open(unit=10, file=mapping_in, & 
        form='unformatted', access='sequential') 
        do n = 1, n_max
            read(unit=10) mapping_matrix(n,:) 
        end do 
    close(unit=10) 

    ! Loading 1D array 
    print *, 'Loading :', data_in  
    open(unit=20, file=data_in, form='unformatted', access='sequential')
    read(unit=20) data_array 
    close(unit=20) 

    ! Storing as a multidimensional array 
    do n = 1, n_max
        i = mapping_matrix(n,1)
        j = mapping_matrix(n,2)
        k = mapping_matrix(n,3)
        mdim_array(i,j,k) = data_array(n) 
    end do
    










end program remapping_data 
