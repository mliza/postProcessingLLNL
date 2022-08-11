#!/opt/homebrew/bin/python3
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
import matplotlib.pyplot as plt 
import IPython 
import sys 
import os 
scripts_path   = os.environ.get('SCRIPTS')
python_scripts = os.path.join(scripts_path, 'Python')
sys.path.append(python_scripts)
sys.path.append('../fortran_modules')
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
data_path    = '../../plate_data/data_15'
pickle_path  = os.path.join(data_path, 'pickle')
temp_path    = os.path.join(data_path, 'temp_data')  
box_path     = os.path.join(data_path, 'BOX')  
saving_path  = '/Users/martin/Desktop/results'
saving_path  = os.path.join(data_path, 'results') 
nx           = 1439
ny           = 89
nz           = 638 
time_step    = '0930000'
U_init       = 3000     #[m/s] 
T_init       = 216.66   #[K] 
RHO_init     = 0.18874  #[kg/m3] 
P_init       = 11737    #[Pa] 
fortran_flag = False  
mapping_flag = False  
writing_flag = False   
working_flag = True  
add_dat_flag = True  
fluct_flag   = False      
scalar_in    = [ 'T', 'RHO', 'P',
                 'RHOE', 'GRADRHOMAG', 
                 'GRADV_11', 'GRADV_12', 'GRADV_13',
                 'GRADV_21', 'GRADV_22', 'GRADV_23',
                 'GRADV_31', 'GRADV_32', 'GRADV_33' ]

# Loading my classes 
helper = helper.Helper()
aero   = aero.Aero()
box    = box.Box(nx=nx, ny=ny, nz=nz)

# Calculate freestream conditions 
sos_init     = aero.speed_of_sound(T_init) 
mu_init      = aero.sutherland_law(T_init) 
mach_init    = U_init / sos_init 
normal_dict  = aero.normal_shock_relations(mach_init) 
oblique_dict = aero.oblique_shock_relations(mach_init, shock_angle_deg=45)  

# Downstream properties, assumes a normal shock wave  
T_2   = T_init * oblique_dict['T_ratio'] #[K]
sos_2 = aero.speed_of_sound(T_2)         #[m/s]
U_2   = oblique_dict['mach_2'] * sos_2   #[m/s]

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
    data_in1D = helper.pickle_manager(pickle_name_file='dict_1D', 
                                      pickle_path=pickle_path)
    data_in3D = helper.pickle_manager(pickle_name_file='dict_3D', 
                                      pickle_path=pickle_path)
    mapping   = helper.pickle_manager(pickle_name_file='mapping', 
                                      pickle_path=pickle_path)
    # Only loads after data is being proceed 
    if not fluct_flag and not add_dat_flag:
        fluct_1D  = helper.pickle_manager(pickle_name_file='fluctuations_dict_1D', 
                                      pickle_path=pickle_path)
        fluct_3D  = helper.pickle_manager(pickle_name_file='fluctuations_dict_3D', 
                                      pickle_path=pickle_path)

    # Calculating mean fields 
    x_mean = box.mean_fields(data_in3D['X'])['mean_x']
    y_mean = box.mean_fields(data_in3D['Y'])['mean_y']
    z_mean = box.mean_fields(data_in3D['Z'])['mean_z'] 
    mean_position_dict = {'mean_x' : x_mean, 
                          'mean_y' : y_mean, 
                          'mean_z' : z_mean} 
    # Calculate boundary layer thickness planes 
    temperature_plane_dict = box.boundary_layer_thickness(data_in3D['T'], 
                                                          data_in3D['Y'],
                                                    freestream_value=T_2)
    velocity_plane_dict    = box.boundary_layer_thickness(data_in3D['Ux'], 
                                                          data_in3D['Y'],
                                                    freestream_value=U_2)
    # Plot boundary layer planes 
    box.plot_boundary_layers(velocity_plane_dict, temperature_plane_dict,
                             mean_position_dict, velocity_freestream=U_2,
                             temperature_freestream=T_2,
                             saving_path=saving_path) 
    #box.wall_shear_stress(velocity_plane_dict, mean_position_dict) 


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

    # Fluctuations 
    if fluct_flag:
        # Key values 
        keys = list(data_in1D.keys()) 
        keys.remove('X')
        keys.remove('Y')
        keys.remove('Z') 
        keys.remove('SoS')
        keys.remove('M')
        fluctuations_1D = { } 
        fluctuations_3D = { } 
        # Calculate Reynolds decomposition 
        for i in keys:
            temp_1D            = box.reynolds_decomposition(data_in1D[i]) 
            fluctuations_1D[i] = temp_1D
            fluctuations_3D[i] = box.split_plot3D(array_1D=temp_1D, 
                                                 mapping=mapping)
        # Turbulent Kinetic Energy 
        turb_k1D = (0.5 * (fluctuations_1D['Ux']**2 + 
                           fluctuations_1D['Uy']**2 + 
                           fluctuations_1D['Uz']**2)) 
        mach_t1D = np.sqrt(2 * turb_k1D) / np.mean(data_in1D['SoS'])  

        fluctuations_1D['K']  = turb_k1D
        fluctuations_3D['K']  = box.split_plot3D(array_1D=turb_k1D,
                                                 mapping=mapping)
        # Turbulent mach number 
        fluctuations_1D['Mt'] = mach_t1D
        fluctuations_3D['Mt'] = box.split_plot3D(array_1D=mach_t1D,
                                                 mapping=mapping)
        # Reynolds stress structure parameters
        rssp_1D                 = box.reynolds_stres_structure_parameters(
                                  fluctuation_1D)
        fluctuations_1D['RSSP'] = rssp_1D
        fluctuations_3D['RSSP'] = box.split_plot3D(array_1D=rssp_1D,
                                                 mapping=mapping)

        # Saving 1D and 3D arrays 
        helper.pickle_manager(pickle_name_file='fluctuations_dict_1D', 
                              pickle_path=pickle_path,
                              data_to_save=fluctuations_1D)
        helper.pickle_manager(pickle_name_file='fluctuations_dict_3D', 
                              pickle_path=pickle_path,
                              data_to_save=fluctuations_3D)

    '''
    # Generate plots  
    box.plot_lineXY(data_in3D, 'Ux', 'Y', x_dim=700, z_dim=300, 
                    saving_path=saving_path) 
    box.plot_lineXY(data_in3D, 'Uy', 'Y', x_dim=700, z_dim=300, 
                    saving_path=saving_path) 
    box.plot_lineXY(data_in3D, 'Z', 'Uz', x_dim=700, y_dim=40, 
                    saving_path=saving_path) 
    box.plot_contour(data_in3D, grid_x='X', grid_y='Z', field='RHO', 
                     slice_cut=40, slice_direction='Y', 
                     levels=500, saving_path=saving_path) 
    box.plot_contour(data_in3D, grid_x='X', grid_y='Z', field='DIL', 
                     slice_cut=40, slice_direction='Y', 
                     levels=500, saving_path=saving_path) 
    box.plot_contour(data_in3D, grid_x='X', grid_y='Z', field='VORTMAG', 
                     slice_cut=40, slice_direction='Y', 
                     levels=500, saving_path=saving_path) 
    box.plot_contour(data_in3D, grid_x='X', grid_y='Z', field='T', 
                     slice_cut=40, slice_direction='Y', 
                     levels=500, saving_path=saving_path) 
    box.plot_contour(data_in3D, grid_x='X', grid_y='Z', field='Ux', 
                     slice_cut=40, slice_direction='Y', 
                     levels=500, saving_path=saving_path) 
    '''
