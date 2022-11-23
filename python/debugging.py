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
import numpy as np 
import matplotlib.pyplot as plt 
scripts_path   = os.environ.get('SCRIPTS')
python_scripts = os.path.join(scripts_path, 'Python')
sys.path.append(python_scripts)

# Loading my classes 
import helper_class as helper 
import box_class as box  

# User Inputs 
data_path         = '/p/lustre1/liza1/dns_results'
pickle_path       = os.path.join(data_path, 'box_pickle')
pickle_results    = os.path.join(data_path, 'sub_pickle')  
nx                = 1439
ny                = 89
nz                = 638 

# Loading my classes 
helper = helper.Helper()
box    = box.Box(nx=nx, ny=ny, nz=nz)
spatial_avg_flag  = False
matrix_name_out   = 'temporal_matrix'
assemble_name_out = 'temporal_average'

# Loading time steps and grid_3D pickle files   
time_steps = helper.pickle_manager(pickle_name_file='time_steps', 
                                   pickle_path=pickle_path)
time_len   = len(time_steps) 

mean_matrix = { }
for count, val in enumerate(time_steps):
    print(count) 
    # Loading dict_3D, fluctuation_3D and rms_2D
    field_3D = helper.pickle_manager(pickle_name_file=f'{val}_dict3D', 
                             pickle_path=pickle_path)
    # Spatial average
    if spatial_avg_flag:
        mean_field_3D  = { }
        for k in field_3D.keys():
            mean_field_3D[k] = box.mean_fields(field_3D[k])['mean_xy']
            # Creates empty matrix at each time [time, nx, ny] 
            if count == 0:
                mean_matrix[k] = np.empty([time_len, nx, ny])
            # Populates spatial average matrix  
            mean_matrix[k][count] = mean_field_3D[k]

    # Doesn't do a spatial average 
    if not spatial_avg_flag:
        for k in field_3D.keys():
            # Creates empty matrix at each time [time, nx, ny, nz] 
            if count == 0:
                mean_matrix[k] = np.empty([time_len, nx, ny, nz])
            # Populates spatial average matrix  
            mean_matrix[k][count] = field_3D[k]
        
# Save spatial average data  
'''
helper.pickle_manager(pickle_name_file=f'{matrix_name_out}',
                      pickle_path=pickle_results,
                      data_to_save=mean_matrix)
'''

# Calculates esemble average  
ensemble_avg = { }
for k in field_3D.keys():
    # Calcualtes time average 
    ensemble_avg[k] = box.time_average(mean_matrix[k])  


# Save ensemble average
helper.pickle_manager(pickle_name_file=f'{assemble_name_out}',
                      pickle_path=pickle_results,
                      data_to_save=ensemble_avg)
