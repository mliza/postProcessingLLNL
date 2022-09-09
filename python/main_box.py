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
import matplotlib
import matplotlib.pyplot as plt 
import pandas as pd 
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
papers_path  = '../../plate_data/papers_data'
pino_path    = os.path.join(papers_path, 'pino_martin') 
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

# Calculate freestream conditions 
sos_init     = aero.speed_of_sound(T_init) 
mu_init      = aero.sutherland_law(T_init) 
mach_init    = U_init / sos_init 
normal_dict  = aero.normal_shock_relations(mach_init) 
oblique_dict = aero.oblique_shock_relations(mach_init, shock_angle_deg=45)  

# Downstream properties, assumes a normal shock wave  
T_2   = T_init * oblique_dict['T_ratio']      #[K]
sos_2 = aero.speed_of_sound(T_2)              #[m/s]
U_2   = oblique_dict['mach_2'] * sos_2        #[m/s]
Rho_2 = RHO_init * oblique_dict['Rho_ratio']  #[kg/m3] 

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
        fluctuations_3D['K']  = (0.5 * ( 
                        box.mean_fields(fluctuations_3D['Ux']**2)['mean_xy'] + 
                        box.mean_fields(fluctuations_3D['Uy']**2)['mean_xy'] + 
                        box.mean_fields(fluctuations_3D['Uz']**2)['mean_xy'] )) 
        fluctuations_3D['Mt'] = (np.sqrt(fluctuations_3D['K']) / 
                                box.mean_fields(data_in3D['SoS'])['mean_xy'])

        # Saving 3D arrays 
        helper.pickle_manager(pickle_name_file='fluctuations_dict_3D', 
                              pickle_path=pickle_path,
                              data_to_save=fluctuations_3D)

    # Adding fluctuation dictionary 
    if rms_flag:
        # Key values 
        keys   = list(fluct_3D.keys()) 
        rms_3D = { } 
        # Calculate Reynolds decomposition 
        for i in keys:
            rms_3D[i] = np.sqrt((fluct_3D[i])**2) 
        
        # Saving 3D arrays 
        helper.pickle_manager(pickle_name_file='rms_dict_3D', 
                              pickle_path=pickle_path,
                              data_to_save=rms_3D)

    # List keys 
    fluct_keys = list(fluct_3D.keys()) 
    data_keys  = list(data_in1D.keys()) 

    # Calculate mean fields  
    loc_mean  = box.mean_positions(data_in3D) 
    X_mean    = box.mean_fields(data_in3D['X'])
    Ux_mean   = box.mean_fields(data_in3D['Ux'])
    rho_mean  = box.mean_fields(data_in3D['RHO'])
    T_mean    = box.mean_fields(data_in3D['T'])
    y_mean    = box.mean_fields(data_in3D['Y'])
    mu_mean   = box.mean_fields(data_in3D['MU'])
    M_mean    = box.mean_fields(data_in3D['M'])
    s12_mean  = box.mean_fields(data_in3D['GRADV_12'])
    Uy_mean   = box.mean_fields(data_in3D['Uy'])
    Uz_mean   = box.mean_fields(data_in3D['Uz'])
    Ux_rms    = box.mean_fields(rms_3D['Ux'])
    Uy_rms    = box.mean_fields(rms_3D['Uy'])
    Uz_rms    = box.mean_fields(rms_3D['Uz'])
    T_rms     = box.mean_fields(rms_3D['T'])
    M_rms     = box.mean_fields(rms_3D['M'])
    grid_dict = { 'X' : data_in3D['X'], 
                  'Y' : data_in3D['Y'],
                  'Z' : data_in3D['Z'] } 
    
    x_        = int(nz/2) 
    y_        = int(ny/2) 
    z_        = int(nz/2)
    str_title = box.str_locations(loc_mean, y=y_, z=z_)
    y_val     = loc_mean['mean_y'][y_] 
    x_val     = loc_mean['mean_x'][x_]
    y_str     = f'y = {y_val:.3} $[m]$'
    x_str     = f'x = {x_val:.3} $[m]$'


    # Longitudinal correlation 
    f_correlation = box.correlation_function(data_in3D['X'][:,y_,z_], 
                                             fluct_3D['Ux'][:,y_,z_],
                                             fluct_3D['Ux'][:,y_,z_],
                                             autocorrelation_len=80) 

    # Transversal correlation 
    g_correlation = box.correlation_function(data_in3D['X'][:,y_,z_], 
                                             fluct_3D['Uy'][:,y_,z_],
                                             fluct_3D['Uy'][:,y_,z_],
                                             autocorrelation_len=70) 

    # Looping energy cascade  
    K_3D = 0.5 * ( fluct_3D['Ux']**2 + fluct_3D['Uy']**2 + fluct_3D['Uz']**2 ) 
                                         
    energy_box = np.empty([nz, 720])
    for k in range(nz):
        energy_box[k] = box.energy_spectrum(data_in3D['X'][:,y_,k], 
                                            K_3D[:,y_,k], 
                                            n_bins=2)
    energy_cascade = np.empty(720)
    for i in range(720):
        energy_cascade[i] = np.mean(energy_box[:,i])
    


    # Van Driest transformation 
    van_driest = box.van_driest(s12_mean, Ux_mean, y_mean, rho_mean, mu_mean)  

    # Calculate wall and edge values  
    temperature_edge = box.edge_properties(data_in3D['T'], data_in3D['Y'],
                                           freestream_value=T_2)
    density_edge     = box.edge_properties(data_in3D['RHO'],data_in3D['Y'],
                                           freestream_value=Rho_2)
    velocity_edge    = box.edge_properties(data_in3D['Ux'],data_in3D['Y'],
                                           freestream_value=U_2)
    temperature_wall = box.wall_properties(data_in3D['T'], data_in3D['Y'])
    density_wall     = box.wall_properties(data_in3D['RHO'], data_in3D['Y'])
    velocity_wall    = box.wall_properties(data_in3D['Ux'], data_in3D['Y'])
    y_plus           = van_driest['mean_y_plus'] 

    # Calculate normalized fields  
    u_rms = Ux_rms['mean_y'] / np.mean(van_driest['u_tau']) 
    v_rms = Uy_rms['mean_y'] / np.mean(van_driest['u_tau'])
    w_rms = Uz_rms['mean_y'] / np.mean(van_driest['u_tau'])
    t_rms = T_rms['mean_y']  / T_mean['mean_y']

    # Plots correlation
    # Energy Cascade plot  
    p_factor = -5/3 
    shifting_factor = 0.85E11
    wave_vector = np.array(range(10, int(np.floor(2 * len(energy_cascade)) / 15))) 
    plt.title(y_str)
    plt.loglog(energy_cascade, 'o-', markerfacecolor='lightgray', 
               linewidth=3, color='k') 
    slope = shifting_factor* wave_vector ** p_factor 
    plt.plot(wave_vector, slope, '--', linewidth='2.5') 
    plt.grid('-.')
    plt.xlim(left=0.9)
    plt.ylabel('$E(\kappa)$')
    plt.xlabel('$\kappa\;\;[1/m]$') 
    plt.tight_layout()
    plt.savefig(f'{saving_path}/energy_spectrum.png', dpi=300) 
    plt.close() 

    # Longitudinal Correlation  
    plt.plot(f_correlation['radius'], f_correlation['norm_correlation'], 'o-', 
             markerfacecolor='lightgray', linewidth='3', color='k') 
    plt.title(str_title)
    correlation_str = f'$L_{{11}} = {f_correlation["integral"]:.3}\n\lambda_{{f}} = {f_correlation["taylor"]:.3}$'
    plt.text(f_correlation['radius'][-1], 1, correlation_str, 
             horizontalalignment='right', verticalalignment='top', 
             bbox=dict(facecolor='white', edgecolor='lightgray',
                                      boxstyle='round,pad=0.2'))
    plt.xlabel('$radius\;\;[m]$')
    plt.ylabel('$correlation\;\;[\;]$')
    plt.grid('-.') 
    plt.tight_layout()
    plt.savefig(f'{saving_path}/longitudinal_correlation.png', dpi=300) 
    plt.close() 

    # Transversal Correlation  
    plt.plot(g_correlation['radius'], g_correlation['norm_correlation'], 'o-', 
             markerfacecolor='lightgray', linewidth='3', color='k') 
    plt.title(str_title)
    plt.xlabel('$radius\;\;[m]$')
    plt.ylabel('$correlation\;\;[\;]$')
    plt.grid('-.') 
    plt.tight_layout()
    plt.savefig(f'{saving_path}/transversal_correlation.png', dpi=300) 
    plt.close() 

    # RMS plots 
    # Plot t_rms
    testing_file = os.path.join(pino_path, 'tRMSnormalized.csv')
    df = pd.read_csv(testing_file)
    y_pino = np.array(df['y_plus']) 
    u_pino  = np.array(df['t_rms']) 
    plt.plot(y_pino, u_pino, '--', linewidth=2, label='pino')

    plt.plot(y_plus, t_rms, linewidth=2, label='$T_{rms}$')
    plt.legend()
    plt.grid('-.')
    plt.xlabel('$y^+$')
    plt.ylabel('$T_{rms}/<T>$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/temperature_rms.png', dpi=300) 
    plt.close() 

    # Plot u, v, w rms 
    # LOADING RMS 
    '''
    testing_file = os.path.join(pino_path, 'uRMSnormalized.csv')
    df = pd.read_csv(testing_file)
    y_pino = np.array(df['y_plus']) 
    u_pino  = np.array(df['u_rms']) 
    plt.plot(y_pino, u_pino, '--', linewidth=2, label='pino')
    '''

    plt.plot(y_plus, u_rms, linewidth=2, label='$u_{rms}/u_{\\tau}$')
    plt.plot(y_plus, v_rms, linewidth=2, label='$v_{rms}/u_{\\tau}$')
    plt.plot(y_plus, w_rms, linewidth=2, label='$w_{rms}/u_{\\tau}$')
    plt.legend()
    plt.grid('-.')
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$U_{rms}/u_{\\tau}$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/velocity_rms.png', dpi=300) 
    plt.close() 

    # Plot Mt and M fluctuations 
    plt.plot(y_plus, M_rms['mean_xy'][x_,:], linewidth=2, label='$M^{\prime}$')
    plt.plot(y_plus, rms_3D['Mt'][x_,:], linewidth=2, label='$M_t^{\prime}$')
    plt.title(x_str)
    plt.legend()
    plt.grid('-.')
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$M^{\prime} \;\;&\;\; M^{\prime}_t\;\;\;[\;]$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/mach_fluctuations.png', dpi=300) 
    plt.close() 

    # Plot Mt 
    plt.plot(y_plus, fluct_3D['Mt'][x_,:], 'o-',
             markerfacecolor='lightgray', linewidth='3', color='k') 
    plt.title(x_str)
    plt.grid('-.')
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$M_t\;\;[\;]$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/machTurbulent.png', dpi=300) 
    plt.close() 

    # Plot M 
    plt.plot(y_plus, M_mean['mean_y'], 'o-', 
             markerfacecolor='lightgray', linewidth='3', color='k') 
    plt.grid('-.')
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$M\;\;[\;]$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/mach.png', dpi=300) 
    plt.close() 

    # Plot K 
    plt.plot(y_plus, fluct_3D['K'][x_,:], 'o-', 
             markerfacecolor='lightgray', linewidth='3', color='k') 
    plt.title(x_str)
    plt.grid('-.')
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$TKE\;\;[m/s]$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/turbulentKineticEnergy.png', dpi=300) 
    plt.close() 

    # Plot T 
    plt.plot(y_plus, T_mean['mean_y'], 'o-',
             markerfacecolor='lightgray', linewidth='3', color='k') 
    plt.grid('-.')
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$T\;\;[K]$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/temperature.png', dpi=300) 
    plt.close() 

    # Generate plots  
    # Plot boundary layer planes 
    box.plot_boundary_layers(velocity_edge, temperature_edge,
                             Ux_mean['mean_y'], T_mean['mean_y'], 
                             loc_mean, 
                             velocity_freestream=U_2,
                             temperature_freestream=T_2,
                             saving_path=saving_path) 

    box.plot_van_driest(van_driest, testing_path=pino_path, 
                        saving_path=saving_path)

    box.plot_contour(data_in3D, grid_dict, grid_x='X', grid_y='Z', field='M', 
                     slice_cut=40, slice_direction='Y', 
                     levels=500, saving_path=saving_path) 
    box.plot_contour(fluct_3D, grid_dict, grid_x='X', grid_y='Z', field='M', 
                     slice_cut=40, slice_direction='Y', 
                     levels=500, saving_path=saving_path) 
    box.plot_contour(data_in3D, grid_dict, grid_x='X', grid_y='Z', field='RHO', 
                     slice_cut=40, slice_direction='Y', 
                     levels=500, saving_path=saving_path) 
    box.plot_contour(data_in3D, grid_dict, grid_x='X', grid_y='Z', field='VORTMAG', 
                     slice_cut=40, slice_direction='Y', 
                     levels=500, saving_path=saving_path) 
    box.plot_contour(data_in3D, grid_dict, grid_x='X', grid_y='Z', field='T', 
                     slice_cut=40, slice_direction='Y', 
                     levels=500, saving_path=saving_path) 
    box.plot_contour(data_in3D, grid_dict, grid_x='X', grid_y='Z', field='Ux', 
                     slice_cut=50, slice_direction='Y', 
                     levels=700, saving_path=saving_path) 
