#!/opt/homebrew/bin/python3.9
'''
    Date:   05/24/2022
    Author: Martin E. Liza
    File:   main_box.py 
    Def:    Main file to process box data. 

    Author		    Date		Revision
    ----------------------------------------------------
    Martin E. Liza	07/22/2022	Initial version.
'''
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
import pandas as pd 
import IPython 
import sys 
import os 
scripts_path   = os.environ.get('SCRIPTS')
python_scripts = os.path.join(scripts_path, 'Python')
fortran_path   = os.path.abspath('../fortran_modules')
sys.path.append(python_scripts)
sys.path.append(fortran_path)
# Helper class 
import helper_class as helper 
import aerodynamics_class as aero
import box_class as box  
import box_plots 
# Fortran subroutines 
import f_mapping
import f_gridReader
import f_scalarReader 
import f_vectorReader

# User Inputs 
data_path    = '../../plate_data/data_15'
pickle_path  = os.path.join(data_path, 'pickle_n')
temp_path    = os.path.join(data_path, 'temp_data_n')  
box_path     = os.path.join(data_path, 'BOX')  
nx           = 1439
ny           = 89
nz           = 638 
time_step    = '0930000' #This will changing for each time step 
fortran_flag = False  
mapping_flag = False  
writing_flag = False   
working_flag = True  
add_dat_flag = False 
fluct_flag   = False 
rms_flag     = False 
scalar_in    = [ 'T', 'RHO', 'P',
                 'RHOE', 'GRADRHOMAG', 
                 'GRADV_11', 'GRADV_12', 'GRADV_13',
                 'GRADV_21', 'GRADV_22', 'GRADV_23',
                 'GRADV_31', 'GRADV_32', 'GRADV_33' ]

# Loading my classes 
helper = helper.Helper()
aero   = aero.Aero()
box    = box.Box(nx=nx, ny=ny, nz=nz)

# Using Fortran subroutines 
# Convert fortran binaries into python readable and creates the mapping file 
n_max = nx * ny * nz
f_mapping.mapping(nx, ny, nz, temp_path)
f_gridReader.grid_reader(n_max, box_path, temp_path, 'U') 
IPython.embed(colors = 'Linux')
#for i in time_step: 
i = time_step 

f_vectorReader.vector_reader(n_max, box_path, temp_path, 'U', i) 

# For HPC 
files_in   = os.listdir(box_path)
steps_in = [idx for idx in files_in if idx.startswith('U')]
steps_in.remove('U.xyz')
steps_in.remove('U.README')
time_steps = [ ] 
for i in steps_in: time_steps.append(i.split('.')[1])







# Use fortran subroutines 
if fortran_flag:
    n_max = nx * ny * nz
    f_mapping.mapping(nx, ny, nz, temp_path)
    f_gridReader.grid_reader(n_max, box_path, temp_path, 'U') 
    f_vectorReader.vector_reader(n_max, box_path, temp_path, 'U', time_step) 
    for i in scalar_in:
        f_scalarReader.scalar_reader(n_max, box_path, 
                                     temp_path, i, time_step)

# Opens mapping fortran output and saves it as a pickle file
if mapping_flag:
    mapping_file_in = os.path.join(temp_path, 'mappingVector.dat') 
    mapping         = box.mapping_reader(mapping_data_in=mapping_file_in, 
                                         pickle_path=pickle_path) 

# Stores data as a 1D and 3D pickle dictionary
if writing_flag:
    var_in = scalar_in + ['X', 'Y', 'Z', 'Ux', 'Uy', 'Uz']
    var_in.sort() 
    # Dictionaries 
    dict_1D = { } 
    dict_3D = { }
    mapping = helper.pickle_manager(pickle_name_file='mapping', 
                                    pickle_path=pickle_path)
    for i in var_in: 
        dict_1D[i] = helper.fortran_data_loader(variable_in=f'{i}.dat',  
                                                abs_path_in=temp_path) 
        # Splits 1D array into a 3D array
        dict_3D[i] = box.split_plot3D(array_1D=dict_1D[i], mapping=mapping)
    # Saving 1D and 3D arrays 
    helper.pickle_manager(pickle_name_file='dict_1D', pickle_path=pickle_path,
                          data_to_save=dict_1D)
    helper.pickle_manager(pickle_name_file='dict_3D', pickle_path=pickle_path,
                          data_to_save=dict_3D)

# Testing and playing around scripts 
if working_flag:
    # Loading dictionaries 
    data_in3D = helper.pickle_manager(pickle_name_file='dict_3D', 
                                      pickle_path=pickle_path)
    mapping   = helper.pickle_manager(pickle_name_file='mapping', 
                                      pickle_path=pickle_path)

    # Only loads after data is being proceed 
    if not fluct_flag and not add_dat_flag: 
        fluct_3D = helper.pickle_manager(pickle_name_file='fluctuations_dict_3D', 
                                      pickle_path=pickle_path)
        if not rms_flag:
            rms_3D   = helper.pickle_manager(pickle_name_file='rms_dict_3D', 
                                      pickle_path=pickle_path)

    # Adding data to the dictionaries 
    if add_dat_flag:
        grad_1D        = box.gradient_fields(data_in1D) 
        grad_1D['MU']  = aero.sutherland_law(data_in1D['T'])
        grad_1D['SoS'] = aero.speed_of_sound(data_in1D['T'])  
        grad_1D['M']   = grad_1D['UMAG'] / grad_1D['SoS']  
        dict_temp1D    = { }
        dict_temp3D    = { }
        for i in grad_1D.keys():
            dict_temp1D[i] = grad_1D[i] 
            dict_temp3D[i] = box.split_plot3D(array_1D=grad_1D[i], 
                                            mapping=mapping)
        # Add new data to dictionaries 
        for i, key in enumerate(dict_temp1D.keys()):
            if i == 0: 
                helper.pickle_dict_add(var_in_data=dict_temp1D[key], 
                                       var_in_str=key,
                                       pickle_path=pickle_path, 
                                       pickle_dict_in='dict_1D',
                                       pickle_dict_out='new_dict_1D')
                helper.pickle_dict_add(var_in_data=dict_temp3D[key], 
                                       var_in_str=key,
                                       pickle_path=pickle_path, 
                                       pickle_dict_in='dict_3D',
                                       pickle_dict_out='new_dict_3D')
            else:
                helper.pickle_dict_add(var_in_data=dict_temp1D[key], 
                                       var_in_str=key,
                                       pickle_path=pickle_path, 
                                       pickle_dict_in='new_dict_1D',
                                       pickle_dict_out='new_dict_1D')
                helper.pickle_dict_add(var_in_data=dict_temp3D[key], 
                                       var_in_str=key,
                                       pickle_path=pickle_path, 
                                       pickle_dict_in='new_dict_3D',
                                       pickle_dict_out='new_dict_3D')

    # Adding fluctuations dictionary  
    if fluct_flag:
        # Key values 
        keys = list(data_in3D.keys()) 
        keys.remove('X')
        keys.remove('Y')
        keys.remove('Z') 
        fluctuations_3D = { } 
        # Calculate Reynolds decomposition 
        for i in keys:
            fluctuations_3D[i] = box.reynolds_decomposition(data_in3D[i]) 

        # Turbulent Kinetic Energy 
        TKE = 0.5 * ( fluctuations_3D['Ux']**2 + 
                      fluctuations_3D['Uy']**2 + 
                      fluctuations_3D['Uz']**2 ) 
        fluctuations_3D['TKE'] = TKE
        fluctuations   = ( 
                        box.mean_fields(fluctuations_3D['Ux']**2)['mean_xy'] + 
                        box.mean_fields(fluctuations_3D['Uy']**2)['mean_xy'] + 
                        box.mean_fields(fluctuations_3D['Uz']**2)['mean_xy'] ) 
        fluctuations_3D['Mt'] = (np.sqrt(fluctuations) / 
                                box.mean_fields(data_in3D['SoS'])['mean_xy'])

        # Saving 3D arrays 
        helper.pickle_manager(pickle_name_file='fluctuations_dict_3D', 
                              pickle_path=pickle_path,
                              data_to_save=fluctuations_3D)

    # Adding fluctuation dictionary 
    if rms_flag:
        # Key values 
        keys   = list(fluct_3D.keys()) 
        keys.remove('Mt') 
        rms_3D = { } 
        # Calculate Reynolds decomposition 
        for i in keys:
            rms_3D[i] = np.sqrt(box.mean_fields(fluct_3D[i]**2)['mean_xy']) 
        
        # Saving 3D arrays 
        helper.pickle_manager(pickle_name_file='rms_dict_3D', 
                              pickle_path=pickle_path,
                              data_to_save=rms_3D)

    # List keys 
    fluct_keys = list(fluct_3D.keys()) 
    data_keys  = list(data_in3D.keys()) 
