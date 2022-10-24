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
import numpy as np 
import matplotlib.pyplot as plt 
scripts_path   = os.environ.get('SCRIPTS')
python_scripts = os.path.join(scripts_path, 'Python')
sys.path.append(python_scripts)

# Helper class 
import helper_class as helper 
import aerodynamics_class as aero
import box_class as box  
import box_plots 

# User Inputs 
data_path         = '/p/lustre1/liza1/dns_margot'
pickle_path       = os.path.join(data_path, 'pickle')
temp_path         = os.path.join(data_path, 'temp_data')  
fluct_pickle_path = os.path.join(data_path, 'fluct_pickle')
rms_pickle_path   = os.path.join(data_path, 'rms_pickle')
results_path      = '/p/lustre1/liza1/dns_results' 
nx                = 1439
ny                = 89
nz                = 638 
U_init            = 3000        #[m/s] 
T_init            = 216.66      #[K] 
RHO_init          = 0.18874     #[kg/m3] 
P_init            = 11737       #[Pa] 
X_init            = 11.5970E-2  #[m]
x_                = int(nx/2) 
y_                = int(ny/2) 
z_                = int(nz/2)
xn                = 'x2'
yn                = 'y2'
energy_name       = f'energy_spectrum_{xn}{yn}'
vanDriest_name    = f'van_driest_{xn}' 
procced_flag      = True 
new_data_flag     = False   
delta_coarse      = 30 

# x_ = [0.06, 0.105, 0.140] => [0, nx/2, 1290]
# y_ = [0.001, 0.00269] => [3, ny/2] 

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
mu_init      = aero.sutherland_law(T_init)
re_init      = U_init * RHO_init / mu_init

# Downstream properties, assumes a shock at 45 
T_2   = T_init * oblique_dict['T_ratio']      #[K]
sos_2 = aero.speed_of_sound(T_2)              #[m/s]
U_2   = oblique_dict['mach_2'] * sos_2        #[m/s]
Rho_2 = RHO_init * oblique_dict['Rho_ratio']  #[kg/m3] 
M_2   = U_2 / sos_2                           #[ ]

# Print statments  
print(f'{x_} = {xn}, {y_} = {yn}') 
print(f'The unit Re is {re_init:.6E} [1/m]')
print(f'The post shock temperature is {T_2:.3} [K]')  
print(f'The post shock mach is {M_2:.2} [ ]')  

# Loading time steps and grid_3D pickle files   
time_steps = helper.pickle_manager(pickle_name_file='time_steps', 
                                   pickle_path=pickle_path)
grid_3D    = helper.pickle_manager(pickle_name_file='grid_3D', 
                                   pickle_path=pickle_path)
mean_grid  = box.mean_positions(grid_3D) 
time_len   = len(time_steps) 

# String Locations for plotting 
yz_str = box.str_locations(mean_grid, x=None, y=y_, z=z_)
xy_str = box.str_locations(mean_grid, x=x_, y=y_, z=None)
x_str  = box.str_locations(mean_grid, x=x_, y=None, z=None)
y_str  = box.str_locations(mean_grid, x=None, y=y_, z=None)

dict_out = { }
for count, val in enumerate(time_steps):
    # Loading dict_3D, fluctuation_3D and rms_2D
    if new_data_flag:
        field_3D = helper.pickle_manager(pickle_name_file=f'{val}_dict3D', 
                                     pickle_path=pickle_path)
        '''
        rms_2D   = helper.pickle_manager(pickle_name_file=f'{val}_rms2D', 
                                     pickle_path=rms_pickle_path)
        '''
        fluct_3D = helper.pickle_manager(pickle_name_file=f'{val}_fluct3D', 
                                     pickle_path=fluct_pickle_path)
    if procced_flag:
        field_3D = helper.pickle_manager(pickle_name_file=f'{val}_dict3D', 
                                     pickle_path=pickle_path)
        proc_3D  = helper.pickle_manager(pickle_name_file=f'{val}_processed',
                                         pickle_path=results_path)
        '''
        fluct_3D = helper.pickle_manager(pickle_name_file=f'{val}_fluct3D', 
                                     pickle_path=fluct_pickle_path)
        '''

        # Calculate energy cascade and Van Driest 
        van_driest     = proc_3D['vanDriest'] 
        energy_cascade = box.energy_spectrum(field_3D['Ux'][x_,y_,:], 
                                             field_3D['Uy'][x_,y_,:],
                                             field_3D['Uz'][x_,y_,:],
                                             n_elements=nz, 
                                             n_bins=2)
        # Initialize at the first iteration
        if count == 0:
            van_driest_matrix            = { }
            velocity_matrix              = np.empty([len(time_steps), nx]) 
            velocity_thickness_matrix    = np.empty([len(time_steps), nx]) 
            temperature_matrix           = np.empty([len(time_steps), nx]) 
            temperature_thickness_matrix = np.empty([len(time_steps), nx]) 
            energy_spectrum_matrix       = np.empty([len(time_steps), 
                                                    np.shape(energy_cascade)[0]]) 
            # Van Driest values 
            for k in van_driest.keys():
                if k.split('_')[1] == 'w':
                    van_driest_matrix[k] = np.empty([len(time_steps), nx]) 
                elif k.split('_')[1] == 'tau':
                    van_driest_matrix[k] = np.empty([len(time_steps), nx]) 
                else:
                    van_driest_matrix[k] = np.empty([len(time_steps),  
                                       np.shape(van_driest[k][x_,:])[0]]) 

        # Create Matrices 
        energy_spectrum_matrix[count] = energy_cascade 
        velocity_matrix[count]        = proc_3D['velocityEdge']['mean_edge_field'] 
        temperature_matrix[count]     = proc_3D['temperatureEdge']['mean_edge_field'] 
        velocity_thickness_matrix[count]    = proc_3D['velocityEdge']['mean_edge_thickness'] 
        temperature_thickness_matrix[count] = proc_3D['temperatureEdge']['mean_edge_thickness'] 
        # Iterates through Van Driest 
        for k in van_driest.keys():
            if k.split('_')[1] == 'w':
                van_driest_matrix[k][count] = van_driest[k]  
            elif k.split('_')[1] == 'tau':
                van_driest_matrix[k][count] = van_driest[k]  
            else:
                van_driest_matrix[k][count] = van_driest[k][x_,:] 

    if new_data_flag:
        # Empty dictionary 
        dict_out = { }
        # Calculate edge values  
        temperature_edge = box.edge_properties(field_3D['T'], grid_3D['Y'],
                                               freestream_value=T_2)
        velocity_edge    = box.edge_properties(field_3D['Ux'], grid_3D['Y'],
                                               freestream_value=U_2)
        # Calculate Van Driest transformation 
        s12_mean   = box.mean_fields(field_3D['GRADV_12'])
        Ux_mean    = box.mean_fields(field_3D['Ux'])
        Y_mean     = box.mean_fields(grid_3D['Y'])
        rho_mean   = box.mean_fields(field_3D['RHO'])
        mu_mean    = box.mean_fields(field_3D['MU'])
        T_mean     = box.mean_fields(field_3D['T'])
        van_driest = box.van_driest(s12_mean, Ux_mean, Y_mean, 
                                    rho_mean, mu_mean, T_mean)  

        # Add data to the dictionary 
        dict_out['velocityEdge']    = velocity_edge 
        dict_out['temperatureEdge'] = temperature_edge 
        dict_out['vanDriest']       = van_driest

        helper.pickle_manager(pickle_name_file=f'{val}_processed',
                              pickle_path=results_path,
                              data_to_save=dict_out)

if procced_flag:
    # Perform time average 
    velocity_mean         = box.time_average(velocity_matrix)
    velocity_thickness    = box.time_average(velocity_thickness_matrix)
    temperature_mean      = box.time_average(temperature_matrix)
    temperature_thickness = box.time_average(temperature_thickness_matrix)
    energy_spectrum_mean  = box.time_average(energy_spectrum_matrix)
    van_driest_mean       = { }
    # Iterates through van driest dictionary  
    for k in van_driest.keys():
        van_driest_mean[k] = box.time_average(van_driest_matrix[k])

    # Return Dictionary 
    dict_out =  { 'velocity_mean'          : velocity_mean, 
                  'velocity_thickness'     : velocity_thickness, 
                  'temperature_mean'       : temperature_mean,
                  'temperature_thickness'  : temperature_thickness,
                  'energy_spectrum'        : energy_spectrum_mean,
                  'energy_spectrum_matrix' : energy_spectrum_matrix, 
                  'van_driest'             : van_driest_mean}

    helper.pickle_manager(pickle_name_file=f'time_average_{xn}{yn}',
                          pickle_path=results_path,
                          data_to_save=dict_out)
    # Plots 
    y_plus_str = f'$y^+$ = {van_driest_mean["y_plus"][y_]/10:.3}'
    y_plus  = van_driest_mean['yi_plus']
    u_plus  = van_driest_mean['ui_plus']

    box.plot_van_driest(y_plus, u_plus, 
                x_str, saving_path=results_path, fig_name=vanDriest_name) 
    box_plots.boundary_layers(velocity_thickness, temperature_thickness, 
                             mean_grid['mean_x'], saving_path=results_path) 
    box_plots.energy_cascade(energy_spectrum_mean[2:], xy_str, y_plus_str, 
                            shifting_factor=2E7,saving_path=results_path, 
                            fig_name=energy_name)

