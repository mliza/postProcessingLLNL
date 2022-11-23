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
data_path         = '/p/lustre1/liza1/dns_results'
pickle_path       = os.path.join(data_path, 'box_pickle')
fluct_pickle_path = os.path.join(data_path, 'fluct_pickle')
rms_pickle_path   = os.path.join(data_path, 'rms_pickle')
results_path      = os.path.join(data_path, 'results')  
pickle_results    = os.path.join(data_path, 'sub_pickle')  
coars_pickle_path = os.path.join(data_path, 'coarser_pickle') 
nx                = 1439
ny                = 89
nz                = 638 
U_init            = 3000        #[m/s] 
T_init            = 216.66      #[K] 
RHO_init          = 0.18874     #[kg/m3] 
P_init            = 11737       #[Pa] 
X_init            = 11.5970E-2  #[m]
'''
x_                = int(0)  
y_                = int(ny/2) 
z_                = int(nz/2)
xn                = 'x1'
yn                = 'y2'
'''
x_                = int(sys.argv[1]) 
y_                = int(sys.argv[2]) 
z_                = int(nz/2)
xn                = sys.argv[3]
yn                = sys.argv[4]
time_avg_flag     = True
new_data_flag     = False
coarser_flag      = False #NEED TO BE MOVE 
f_width           = 34   #NEED TO BE MOVE 

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
oblique_dict = aero.oblique_shock_relations(mach_init, shock_angle_deg=26)  
mu_init      = aero.sutherland_law(T_init)
re_init      = U_init * RHO_init / mu_init

# Downstream properties, assumes a shock at 30 
T_2   = T_init * oblique_dict['T_ratio']      #[K]
P_2   = P_init * oblique_dict['P_ratio']      #[Pa]
sos_2 = aero.speed_of_sound(T_2)              #[m/s]
U_2   = oblique_dict['mach_2'] * sos_2        #[m/s]
Rho_2 = RHO_init * oblique_dict['Rho_ratio']  #[kg/m3] 
M_2   = U_2 / sos_2                           #[ ]
T_2 = 1510
U_2 = 2531 


# Print statments  
print(f'{x_} = {xn}, {y_} = {yn}') 
'''
print(f'The unit Re is {re_init:.6E} [1/m]')
print(f'The post shock mach is {M_2:.2} [ ]')  
print(f'The post shock pressure is {P_2:.2} [Pa]')  
print(f'The post shock temperature is {T_2:.3} [K]')  
print(f'The post shock density is {Rho_2:.3} [kg/m^3]')  
'''


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
    if time_avg_flag:
        field_3D = helper.pickle_manager(pickle_name_file=f'{val}_dict3D', 
                                     pickle_path=pickle_path)
        rms_2D   = helper.pickle_manager(pickle_name_file=f'{val}_rms2D', 
                                     pickle_path=rms_pickle_path)
        proc_3D  = helper.pickle_manager(pickle_name_file=f'{val}_processed',
                                     pickle_path=pickle_results)
        fluct_3D = helper.pickle_manager(pickle_name_file=f'{val}_fluct3D', 
                                     pickle_path=fluct_pickle_path)

        # Calculate energy cascade and Van Driest 
        van_driest     = proc_3D['vanDriest'] 
        energy_cascade = box.energy_spectrum(field_3D['Ux'][x_,y_,:], 
                                             field_3D['Uy'][x_,y_,:],
                                             field_3D['Uz'][x_,y_,:],
                                             n_elements=nz, 
                                             n_bins=2)
        # Dilatation, Shear and rotation
        dilatation = box.mean_fields(field_3D['dilatation_norm'])['mean_xy']
        rotation   = box.mean_fields(field_3D['rotation_norm'])['mean_xy']
        shear      = box.mean_fields(field_3D['shear_norm'])['mean_xy']
        DIL        = box.mean_fields(field_3D['DIL'])['mean_xy']
        rho        = box.mean_fields(field_3D['RHO'])['mean_xy']
        

        # Initialize at the first iteration
        if count == 0:
            time_len                     = len(time_steps)
            van_driest_matrix            = { }
            Ux_rms_matrix                = np.empty([time_len, ny])
            Uy_rms_matrix                = np.empty([time_len, ny])
            Uz_rms_matrix                = np.empty([time_len, ny])
            Mt_matrix                    = np.empty([time_len, ny])  
            M_matrix                     = np.empty([time_len, ny])
            dilatation_matrix            = np.empty([time_len, ny])
            rotation_matrix              = np.empty([time_len, ny])
            DIL_matrix                   = np.empty([time_len, ny])
            shear_matrix                 = np.empty([time_len, ny])
            velocity_matrix              = np.empty([time_len, nx]) 
            velocity_thickness_matrix    = np.empty([time_len, nx]) 
            temperature_matrix           = np.empty([time_len, nx]) 
            rho_matrix                   = np.empty([time_len, nx]) 
            temperature_thickness_matrix = np.empty([time_len, nx]) 
            energy_spectrum_matrix       = np.empty([time_len, 
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
        Mt_matrix[count]              = rms_2D['Mt'][x_]
        M_matrix[count]               = rms_2D['M'][x_]
        Ux_rms_matrix[count]          = rms_2D['Ux'][x_]
        Uy_rms_matrix[count]          = rms_2D['Uy'][x_]
        Uz_rms_matrix[count]          = rms_2D['Uz'][x_]
        dilatation_matrix[count]      = dilatation[x_]
        rotation_matrix[count]        = rotation[x_]
        rho_matrix[count]             = rho[x_,:]
        DIL_matrix[count]             = DIL[x_,:]
        shear_matrix[count]           = shear[x_,:]
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
        field_3D = helper.pickle_manager(pickle_name_file=f'{val}_dict3D', 
                                         pickle_path=pickle_path)
        # Empty dictionary 
        dict_out = { }
        # Calculate edge values  
        # TESTING ML #
        '''
        temperature_edge = box.edge_properties(field_3D['T'], grid_3D['Y'],
                                               freestream_value=T_2)
        velocity_edge    = box.edge_properties(field_3D['Ux'], grid_3D['Y'],
                                               freestream_value=U_2)
        '''
        temperature_edge = box.edge_properties_mean(box.mean_fields(field_3D['T'])['mean_xy'], 
                                               mean_grid['mean_y'],
                                               freestream_value=T_2)
        velocity_edge    = box.edge_properties_mean(box.mean_fields(field_3D['Ux'])['mean_xy'], 
                                               mean_grid['mean_y'],
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
                              pickle_path=pickle_results,
                              data_to_save=dict_out)
    # Creating coarser field 
    if coarser_flag:
        field_3D = helper.pickle_manager(pickle_name_file=f'{val}_dict3D', 
                                         pickle_path=pickle_path)
        dict_out = { }
        for key in field_3D:
            dict_out[key] = box.coarser_field(field_3D[key], f_width=f_width)

    # Save in a pickle file 
    helper.pickle_manager(pickle_name_file=f'{val}_coarse_box_{f_width}',
                          pickle_path=coars_pickle_path,
                          data_to_save=dict_out)

if time_avg_flag:
    # Perform time average 
    velocity_mean         = box.time_average(velocity_matrix)
    velocity_thickness    = box.time_average(velocity_thickness_matrix)
    temperature_mean      = box.time_average(temperature_matrix)
    rho_mean              = box.time_average(rho_matrix)
    temperature_thickness = box.time_average(temperature_thickness_matrix)
    energy_spectrum_mean  = box.time_average(energy_spectrum_matrix)
    Mt_mean               = box.time_average(Mt_matrix)
    M_mean                = box.time_average(M_matrix)
    Ux_rms_mean           = box.time_average(Ux_rms_matrix)
    Uy_rms_mean           = box.time_average(Uy_rms_matrix)
    Uz_rms_mean           = box.time_average(Uz_rms_matrix)
    dilatation_mean       = box.time_average(dilatation_matrix)
    rotation_mean         = box.time_average(rotation_matrix)
    shear_mean            = box.time_average(shear_matrix)
    DIL_mean              = box.time_average(DIL_matrix)
                
    van_driest_mean       = { }
    # Iterates through van driest dictionary  
    for k in van_driest.keys():
        van_driest_mean[k] = box.time_average(van_driest_matrix[k])

    # Location dictionary
    location = { 'index' : [x_, y_, z_],
                 'x_loc' : mean_grid['mean_x'][x_], 
                 'y_loc' : mean_grid['mean_y'][y_], 
                 'y_p'   : np.round(van_driest_mean['y_plus'][y_]/10,3), 
                 'z_loc' : mean_grid['mean_z'][z_] } 

    # Return Dictionary 
    dict_out =  { 'location'               : location, 
                  'velocity_mean'          : velocity_mean, 
                  'velocity_thickness'     : velocity_thickness, 
                  'temperature_mean'       : temperature_mean,
                  'rho_mean'               : rho_mean, 
                  'temperature_thickness'  : temperature_thickness,
                  'energy_spectrum'        : energy_spectrum_mean,
                  'van_driest'             : van_driest_mean,
                  'Mt'                     : Mt_mean, 
                  'M_rms'                  : M_mean,
                  'Ux_rms'                 : Ux_rms_mean,
                  'Uy_rms'                 : Uy_rms_mean,
                  'Uz_rms'                 : Uz_rms_mean,
                  'rotation_mean'          : rotation_mean,
                  'dilatation_mean'        : dilatation_mean,
                  'DIL'                    : DIL_Mean, 
                  'shear_mean'             : shear_mean }
                  
    helper.pickle_manager(pickle_name_file=f'time_average_{xn}{yn}',
                          pickle_path=pickle_results,
                          data_to_save=dict_out)
    # Plots 
    y_plus_str = f'$y^+$ = {van_driest_mean["y_plus"][y_]/10:.3}'
    y_plus   = van_driest_mean['y_plus']
    u_plus   = van_driest_mean['u_plus']
    rho_plus = van_driest_mean['rho_plus'] 
    T_plus   = van_driest_mean['T_plus'] 

    box.plot_van_driest(y_plus, u_plus, 
                x_str, saving_path=results_path, fig_name=f'van_driest_{xn}') 
    box_plots.y_plus(y_plus, Mt_mean, 'M_t', x_str, saving_path=results_path,
                     fig_name=f'Mt_{xn}') 
    box_plots.y_plus(y_plus, T_plus, 'T^+', xy_str, saving_path=results_path,
                     fig_name=f'T_plus_{xn}') 
    box_plots.y_plus(y_plus, rho_plus, '\\rho^+', xy_str, saving_path=results_path,
                     fig_name=f'rho_plus_{xn}') 
    box_plots.boundary_layers(velocity_thickness, temperature_thickness, 
                             mean_grid['mean_x'], saving_path=results_path) 
    box_plots.energy_cascade(energy_spectrum_mean[2:], xy_str, y_plus_str, 
                        shifting_factor=2E7,saving_path=results_path, 
                        fig_name=f'energy_spectrum_{xn}{yn}')
