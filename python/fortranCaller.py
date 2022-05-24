#!/usr/bin/python3
import sys
import os 
sys.path.append('../fortran_modules') 
import f_mapping 
import f_gridReader 
import f_scalarReader
import f_vectorReader


path_in   = '../../plate_data/data_12'
path_temp = os.path.join(path_in, 'temp_data')
box_data  = os.path.join(path_in, 'smallBOX_BIN') 
grid_in   = os.path.join(box_data, 'T.xyz')
mapping_out = os.path.join(path_temp, 'mappingVector.dat')
nx    = 5
ny    = 3
nz    = 2 
n_max = nx * ny * nz 

# Call fortran routines 
f_mapping.mapping(nx, ny, nz, path_temp) 
f_gridReader.grid_reader(n_max, box_data, path_temp, 'T') 
f_scalarReader.scalar_reader(n_max, box_data, path_temp, 'T', '0848001')
f_vectorReader.vector_reader(n_max, box_data, path_temp, 'U', '0848001')

