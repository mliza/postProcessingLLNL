#!/usr/bin/python3
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
path_in   = '../../plate_data/data_12'
path_temp = os.path.join(path_in, 'temp_data')
path_box  = os.path.join(path_in, 'smallBOX_BIN') 
nx        = 5
ny        = 3
nz        = 2 
n_max     = nx * ny * nz 

# Calls fortran subroutines  
f_mapping.mapping(nx, ny, nz, path_temp) 
f_gridReader.grid_reader(n_max, path_box, path_temp, 'T') 
f_scalarReader.scalar_reader(n_max, path_box, path_temp, 'T', '0848001')
f_vectorReader.vector_reader(n_max, path_box, path_temp, 'U', '0848001')
