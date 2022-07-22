#!/opt/homebrew/bin/python3
'''
    Date:   05/24/2022
    Author: Martin E. Liza
    File:   fortran_modules.py 
    Def:    Wraps Fortran code and uses to process box data on plot3D format. 

    Author		    Date		Revision
    ----------------------------------------------------
    Martin E. Liza	05/24/2022	Initial version.
'''
import os 
import sys
sys.path.append('../fortran_modules') 
# Fortran subroutines 
import f_mapping 
import f_gridReader 
import f_scalarReader
import f_vectorReader

# Path definitions and user inputs  
path_in   = '../../plate_data/data_11'
path_temp = os.path.join(path_in, 'temp_data')
path_box  = os.path.join(path_in, 'BOX') 
nx        = 1439
ny        = 85
nz        = 638 
n_max     = nx * ny * nz 
scalar    = [ 'GRADRHOMAG', 'P', 
              'RHO', 'RHOE', 'T', 
              'GRADV_11', 'GRADV_12', 'GRADV_13', 
              'GRADV_21', 'GRADV_22', 'GRADV_23',
              'GRADV_31', 'GRADV_32', 'GRADV_33' ]
time_stp  = '0900000'

# Calls fortran subroutines  
f_mapping.mapping(nx, ny, nz, path_temp) 
f_gridReader.grid_reader(n_max, path_box, path_temp, 'U') 
f_vectorReader.vector_reader(n_max, path_box, path_temp, 'U', time_stp)
for i in scalar: 
    f_scalarReader.scalar_reader(n_max, path_box, path_temp, i, time_stp)
