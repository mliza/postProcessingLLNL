#!/opt/homebrew/bin/python3.9
'''
    Date:   09/27/2022
    Author: Martin E. Liza
    File:   main_box.py 
    Def:    Main file to process box data. 

    Author		    Date		Revision
    ----------------------------------------------------
    Martin E. Liza	07/22/2022	Initial version.
'''
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
# Fortran subroutines 
import f_mapping
import f_gridReader
import f_scalarReader 
import f_vectorReader

# User Inputs 
data_path         = '/p/lustre1/liza1/dns_margot'
pickle_path       = os.path.join(data_path, 'pickle')
temp_path         = os.path.join(data_path, 'temp_data')  
fluct_pickle_path = os.path.join(data_path, 'fluct_pickle')
rms_pickle_path   = os.path.join(data_path, 'rms_pickle')
box_path          = os.path.join(data_path, 'BOX')  
nx                = 1439
ny                = 89
nz                = 638 
mapping_flag      = False  
grid_flag         = False
time_flag         = False  
a                 = 85
b                 = 84

# Loading my classes 
helper = helper.Helper()
aero   = aero.Aero()
box    = box.Box(nx=nx, ny=ny, nz=nz)

# Using Fortran subroutines 
# Convert fortran binaries into python readable and creates the mapping file 
n_max = nx * ny * nz
# Run the mapping flag only once 
if mapping_flag:
    f_mapping.mapping(nx, ny, nz, temp_path)
# Run the grid flag only once 
if grid_flag:
    f_gridReader.grid_reader(n_max, box_path, temp_path, 'U') 

# List all files in box folder 
files_in = os.listdir(box_path)

# List time steps and creates a vector 
if time_flag:
    steps_lst = [idx for idx in files_in if idx.startswith('U')]
    steps_lst.remove('U.xyz')
    steps_lst.remove('U.README')
    time_steps = [ ] 
    for i in steps_lst: time_steps.append(i.split('.')[1])
    time_steps.sort() 
    # Saving a time vector  
    helper.pickle_manager(pickle_name_file='time_steps', 
                              pickle_path=pickle_path,
                              data_to_save=time_steps)
else:
    time_steps = helper.pickle_manager(pickle_name_file='time_steps', 
                                       pickle_path=pickle_path)

# For running multiple times 
time_steps = time_steps[a:b]  # cutting processing time 

# List scalar fields and create as a vector  
scalar_fields = [ ]
scalar_lst    = [idx for idx in files_in if time_steps[0] in idx]
for i in scalar_lst: scalar_fields.append(i.split('.')[0])
scalar_fields.remove('U')
scalar_fields.append('Ux')
scalar_fields.append('Uy')
scalar_fields.append('Uz')
scalar_fields.sort() 

# Loading and saving mapping (run this once) 
if mapping_flag:
    print('Loading fortran mapping and saving as a pickle file')
    mapping_file_in = os.path.join(temp_path, 'mappingVector.dat') 
    mapping         = box.mapping_reader(mapping_data_in=mapping_file_in, 
                                         pickle_path=pickle_path) 
else:
    print('Loading pickle mapping') 
    mapping = helper.pickle_manager(pickle_name_file='mapping', 
                                    pickle_path=pickle_path)

# Create grid dictionary with the positions (run this once)
if grid_flag:
    grid_in   = ['X', 'Y', 'Z']
    grid_dict = { }
    print('Processing grid data')
    for i in grid_in:
        print(f'Processing {i}') 
        tmp_array1D  = helper.fortran_data_loader(variable_in=f'{i}.dat',  
                                                  abs_path_in=temp_path) 
        grid_dict[i] = box.split_plot3D(array_1D=tmp_array1D, mapping=mapping)
        os.remove(os.path.join(temp_path, f'{i}.dat')) 
    # Saving grid as a pickle file 
    helper.pickle_manager(pickle_name_file='grid_3D', 
                          pickle_path=pickle_path,
                          data_to_save=grid_dict)

# Iterates throught time steps and variables for pickle, fluctuations and rms 
for i in time_steps:
    f_vectorReader.vector_reader(n_max, box_path, temp_path, 'U', i) 
    dict_3D = { }
    for j in scalar_fields:
        # Creates temp files 
        if j[0] != 'U': 
            f_scalarReader.scalar_reader(n_max, box_path, 
                                     temp_path, j, i)

        # Loads 1D array, splits data in 3D and saves it
        tmp_array1D = helper.fortran_data_loader(variable_in=f'{i}_{j}.dat',  
                                                 abs_path_in=temp_path) 
        print(f'Splitting raw data: {i}_{j}')
        dict_3D[j]  = box.split_plot3D(array_1D=tmp_array1D, mapping=mapping)
        os.remove(os.path.join(temp_path, f'{i}_{j}.dat')) 
    
    print(f'Processing gradient data: {i}')
    grad_3D        = box.gradient_fields(dict_3D) 
    dict_3D['MU']  = aero.sutherland_law(dict_3D['T'])
    grad_3D['SoS'] = aero.speed_of_sound(dict_3D['T'])
    grad_3D['M']   = grad_3D['UMAG'] / grad_3D['SoS']  
    dict_3D.update(grad_3D)

    # Save dictionary as a 3D field and add gradient quantities to it 
    helper.pickle_manager(pickle_name_file=f'{i}_dict3D', 
                          pickle_path=pickle_path,
                          data_to_save=dict_3D)

    # Fluctuation data and rms data  
    var_keys = list(dict_3D.keys())
    fluct_3D = { } 
    rms_2D   = { }
    print(f'Processing fluctuation and rms data: {i}')
    for k in var_keys:
        print(f'Fluctuating data: {i}_{k}')
        fluct_3D[k] = box.reynolds_decomposition(dict_3D[k])

        print(f'RMS data: {i}_{k}')
        rms_2D[k] = box.mean_fields(fluct_3D[k]**2)['mean_xy'] 

    # Turbulent Kinetic enrgy and Mt  in rms_2D
    rms_2D['TKE'] = 0.5 * (rms_2D['Ux'] + rms_2D['Uy'] + rms_2D['Uz']) 
    rms_2D['Mt']  = np.sqrt(2 * rms_2D['TKE']) / rms_2D['SoS']  

    # Saving fluctuation dictionaries 
    helper.pickle_manager(pickle_name_file=f'{i}_fluct3D', 
                          pickle_path=fluct_pickle_path, 
                          data_to_save=fluct_3D)
    # Saving RMS dictionaries 
    helper.pickle_manager(pickle_name_file=f'{i}_rms2D', 
                          pickle_path=rms_pickle_path,
                          data_to_save=rms_2D)
