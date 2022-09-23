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
import box_plots 
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
U_init       = 3000        #[m/s] 
T_init       = 216.66      #[K] 
RHO_init     = 0.18874     #[kg/m3] 
P_init       = 11737       #[Pa] 
X_init       = 11.5970E-2  #[m]
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
mu_init      = aero.sutherland_law(T_init)
re_init      = U_init * RHO_init / mu_init
print(f'The unit Re is {re_init:.6E} [1/m]')


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

    # Calculate mean fields  
    loc_mean  = box.mean_positions(data_in3D) 
    X_mean    = box.mean_fields(data_in3D['X'])
    Ux_mean   = box.mean_fields(data_in3D['Ux'])
    rho_mean  = box.mean_fields(data_in3D['RHO'])
    T_mean    = box.mean_fields(data_in3D['T'])
    Y_mean    = box.mean_fields(data_in3D['Y'])
    mu_mean   = box.mean_fields(data_in3D['MU'])
    M_mean    = box.mean_fields(data_in3D['M'])
    s12_mean  = box.mean_fields(data_in3D['GRADV_12'])
    Uy_mean   = box.mean_fields(data_in3D['Uy'])
    Uz_mean   = box.mean_fields(data_in3D['Uz'])
    Ux_rms    = rms_3D['Ux'] 
    Uy_rms    = rms_3D['Uy'] 
    Uz_rms    = rms_3D['Uz'] 
    T_rms     = rms_3D['T'] 
    M_rms     = rms_3D['M'] 
    TKE_mean  = box.mean_fields(fluct_3D['TKE'])
    grid_dict = { 'X' : data_in3D['X'], 
                  'Y' : data_in3D['Y'],
                  'Z' : data_in3D['Z'] } 
    
    # Titles 
    x_     = int(nx/2) 
    y_     = int(ny/2) 
    z_     = int(nz/2)
    yz_str = box.str_locations(loc_mean, x=None, y=y_, z=z_)
    xy_str = box.str_locations(loc_mean, x=x_, y=y_, z=None)
    y_str  = box.str_locations(loc_mean, x=None, y=y_, z=None)  
    x_str  = box.str_locations(loc_mean, x=x_, y=None, z=None) 

    # Longitudinal correlation 
    f_correlation = box.correlation_function(data_in3D['X'][:,y_,z_], 
                                             fluct_3D['Ux'][:,y_,z_],
                                             fluct_3D['Ux'][:,y_,z_],
                                             autocorrelation_len=100) 

    # Transversal correlation 
    g_correlation = box.correlation_function(data_in3D['Z'][x_,y_,:], 
                                             fluct_3D['Ux'][x_,y_,:],
                                             fluct_3D['Ux'][x_,y_,:],
                                             autocorrelation_len=100) 

    # Looping energy cascade  
    energy_cascade = box.energy_spectrum(data_in3D['X'][x_,y_,:],
                                         fluct_3D['TKE'][x_,y_,:],
                                         n_bins = 2) 


    # Van Driest transformation 
    van_driest = box.van_driest(s12_mean, Ux_mean, Y_mean, 
                                rho_mean, mu_mean, T_mean)  
    y_plus_str = f'$y^+$ = {van_driest["y_plus"][x_,y_]/10:.3}'


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

    ## PLOTS ##  
    box_plots.energy_cascade(energy_cascade, xy_str, y_plus_str, p_factor=-5/3, 
                             shifting_factor=8.62E11,saving_path=saving_path)
    box_plots.correlation(f_correlation, 'f', xy_str, saving_path=saving_path)
    box_plots.correlation(g_correlation, 'g', xy_str, saving_path=saving_path)
    box_plots.boundary_layers(velocity_edge, temperature_edge,
                             loc_mean, saving_path=saving_path) 
    IPython.embed(colors='Linux') 

    # Calculate normalized fields  
    u_rms = Ux_rms[x_,:] / van_driest['u_tau'][x_] 
    v_rms = Uy_rms[x_,:] / van_driest['u_tau'][x_]
    w_rms = Uz_rms[x_,:] / van_driest['u_tau'][x_]
    t_rms = T_rms[x_,:]  / T_mean['mean_xy'][x_,:]
    y_wall   = Y_mean['mean_xy'][x_,:] / np.max(Y_mean['mean_xy'][x_,:]) 
    y_plus   = van_driest['y_plus'][x_,:] 
    T_plus   = van_driest['T_plus'][x_,:]
    rho_plus = van_driest['rho_plus'][x_,:]

    # RMS plots 
    # Plot t_rms
    plt.plot(y_plus, t_rms, 
             'o-' ,markerfacecolor='lightgray', linewidth='3', color='k',        
             label='$T_{rms}$') 
    plt.title(x_str)
    plt.xscale('log') 
    plt.legend()
    plt.grid('-.')
    plt.xlabel('$y^+$')
    plt.ylabel('$T_{rms}/<T>$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/temperature_rms.png', dpi=300) 
    plt.close() 

    # Plot u, v, w rms 
    # LOADING RMS 
    testing_M5  = os.path.join(pino_path, 'uRMS_3fig6_M5.csv')
    testing_M8  = os.path.join(pino_path, 'uRMS_3fig6_M8.csv')
    testing_M12 = os.path.join(pino_path, 'uRMS_3fig6_M12.csv')
    df_M5  = pd.read_csv(testing_M5)
    df_M8  = pd.read_csv(testing_M8)
    df_M12 = pd.read_csv(testing_M12)
    y_pino_M5  = np.array(df_M5['y_wall']) 
    u_pino_M5  = np.array(df_M5['u_rms']) 
    y_pino_M8  = np.array(df_M8['y_wall']) 
    u_pino_M8  = np.array(df_M8['u_rms']) 
    y_pino_M12 = np.array(df_M12['y_wall']) 
    u_pino_M12 = np.array(df_M12['u_rms']) 
    plt.plot(y_wall, u_rms,
             'o-' ,markerfacecolor='lightgray', linewidth='3', color='k',        
             label='MARGOT, M10')
    plt.plot(y_pino_M5, u_pino_M5, '--', label='Duan, M5')
    plt.plot(y_pino_M8, u_pino_M8, '--', label='Duan, M8')
    plt.plot(y_pino_M12, u_pino_M12, '--', label='Duan, M12')
    plt.title(x_str)
    plt.legend()
    plt.grid('-.')
    plt.xlabel('$y/\\delta$')
    plt.ylabel('$u_{rms}/u_{\\tau}$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/u_rms.png', dpi=300) 
    plt.close() 

    testing_M5 = os.path.join(pino_path, 'vRMS_3fig6_M5.csv')
    testing_M8 = os.path.join(pino_path, 'vRMS_3fig6_M8.csv')
    testing_M12 = os.path.join(pino_path, 'vRMS_3fig6_M12.csv')
    df_M5 = pd.read_csv(testing_M5)
    df_M8 = pd.read_csv(testing_M8)
    df_M12 = pd.read_csv(testing_M12)
    y_pino_M5  = np.array(df_M5['y_wall']) 
    v_pino_M5  = np.array(df_M5['v_rms']) 
    y_pino_M8  = np.array(df_M8['y_wall']) 
    v_pino_M8  = np.array(df_M8['v_rms']) 
    y_pino_M12 = np.array(df_M12['y_wall']) 
    v_pino_M12 = np.array(df_M12['v_rms']) 
    plt.plot(y_wall, v_rms, 
             'o-' ,markerfacecolor='lightgray', linewidth='3', color='k',        
             label='MARGOT, M10')
    plt.plot(y_pino_M5, v_pino_M5, '--', label='Duan, M5')
    plt.plot(y_pino_M8, v_pino_M8, '--', label='Duan, M8')
    plt.plot(y_pino_M12, v_pino_M12, '--', label='Duan, M12')
    plt.title(x_str)
    plt.legend()
    plt.grid('-.')
    plt.xlabel('$y/\\delta$')
    plt.ylabel('$v_{rms}/u_{\\tau}$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/v_rms.png', dpi=300) 
    plt.close() 

    testing_M5 = os.path.join(pino_path, 'wRMS_3fig6_M5.csv')
    testing_M8 = os.path.join(pino_path, 'wRMS_3fig6_M8.csv')
    testing_M12 = os.path.join(pino_path, 'wRMS_3fig6_M12.csv')
    df_M5 = pd.read_csv(testing_M5)
    df_M8 = pd.read_csv(testing_M8)
    df_M12 = pd.read_csv(testing_M12)
    y_pino_M5  = np.array(df_M5['y_wall']) 
    w_pino_M5  = np.array(df_M5['w_rms']) 
    y_pino_M8  = np.array(df_M8['y_wall']) 
    w_pino_M8  = np.array(df_M8['w_rms']) 
    y_pino_M12 = np.array(df_M12['y_wall']) 
    w_pino_M12 = np.array(df_M12['w_rms']) 
    plt.plot(y_wall, w_rms,
             'o-' ,markerfacecolor='lightgray', linewidth='3', color='k',        
             label='MARGOT, M10')
    plt.plot(y_pino_M5, w_pino_M5, '--', label='Duan, M5')
    plt.plot(y_pino_M8, w_pino_M8, '--', label='Duan, M8')
    plt.plot(y_pino_M12, w_pino_M12, '--', label='Duan, M12')
    plt.title(x_str)
    plt.legend()
    plt.grid('-.')
    plt.xlabel('$y/\\delta$')
    plt.ylabel('$w_{rms}/u_{\\tau}$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/w_rms.png', dpi=300) 
    plt.close() 

    # Plot M fluctuations 
    testing_M5 = os.path.join(pino_path, 'machNumberFluctuation_1fig11_M5.csv')
    testing_M8 = os.path.join(pino_path, 'machNumberFluctuation_3fig14_M8.csv')
    testing_M12 = os.path.join(pino_path, 'machNumberFluctuation_3fig14_M12.csv')
    df_M5 = pd.read_csv(testing_M5)
    df_M8 = pd.read_csv(testing_M8)
    df_M12 = pd.read_csv(testing_M12)
    y_M5 = np.array(df_M5['y_plus']) 
    mach_M5  = np.array(df_M5['mach']) 
    y_M8 = np.array(df_M8['y_plus']) 
    mach_M8  = np.array(df_M8['mach']) 
    y_M12 = np.array(df_M12['y_plus']) 
    mach_M12  = np.array(df_M12['mach']) 
    plt.plot(y_plus, M_rms[x_,:], 
             'o-' ,markerfacecolor='lightgray', linewidth='3', color='k',        
             label='MARGOT, M10')
    plt.plot(y_M5, mach_M5, '--', label='Pino, M5')
    plt.plot(y_M8, mach_M8, '--', label='Duan, M8')
    plt.plot(y_M12, mach_M12, '--', label='Duan, M12')
    plt.title(x_str)
    plt.legend()
    plt.grid('-.')
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$M^{\prime}_{rms}\;\;\;[\;]$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/mach_fluctuations.png', dpi=300) 
    plt.close() 

    # Plot Mt 
    plt.plot(y_plus, fluct_3D['Mt'][x_,:], 'o-',
             markerfacecolor='lightgray', linewidth='3', color='k', label='MARGOT, M10') 
    testing_M5  = os.path.join(pino_path, 'turbulentMachNumber_3fig14_M5.csv')
    testing_M8  = os.path.join(pino_path, 'turbulentMachNumber_3fig14_M8.csv')
    testing_M12 = os.path.join(pino_path, 'turbulentMachNumber_3fig14_M12.csv')
    df_M5 = pd.read_csv(testing_M5)
    df_M8 = pd.read_csv(testing_M8)
    df_M12 = pd.read_csv(testing_M12)
    y_M5 = np.array(df_M5['y_plus']) 
    mach_M5  = np.array(df_M5['turbulent_mach']) 
    y_M8 = np.array(df_M8['y_plus']) 
    mach_M8  = np.array(df_M8['turbulent_mach']) 
    y_M12 = np.array(df_M12['y_plus']) 
    mach_M12  = np.array(df_M12['turbulent_mach']) 
    plt.plot(y_M5, mach_M5, '--', label='Duan, M5')
    plt.plot(y_M8, mach_M8, '--', label='Duan, M8')
    plt.plot(y_M12, mach_M12,'--',label='Duan, M12')
    plt.title(x_str)
    plt.legend()
    plt.grid('-.')
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$M_t\;\;[\;]$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/machTurbulent.png', dpi=300) 
    plt.close() 

    # Plot M 
    plt.plot(y_plus, M_mean['mean_xy'][x_,:], 'o-', 
             markerfacecolor='lightgray', linewidth='3', color='k') 
    plt.title(x_str)
    plt.grid('-.')
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$M\;\;[\;]$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/mach.png', dpi=300) 
    plt.close() 

    # Plot K 
    plt.plot(y_plus, TKE_mean['mean_xy'][x_,:], 'o-', 
             markerfacecolor='lightgray', linewidth='3', color='k') 
    plt.title(xy_str)
    plt.grid('-.')
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$TKE\;\;[m/s]$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/turbulentKineticEnergy.png', dpi=300) 
    plt.close() 

    # Plot T 
    plt.plot(y_plus, T_plus, 'o-',
             markerfacecolor='lightgray', linewidth='3', color='k') 
    plt.title(x_str)
    plt.grid('-.')
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$T^+$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/temperature_plus.png', dpi=300) 
    plt.close() 

    # Plot Rho 
    plt.plot(y_plus, rho_plus, 'o-',
             markerfacecolor='lightgray', linewidth='3', color='k') 
    plt.title(x_str)
    plt.grid('-.')
    plt.xscale('log')
    plt.xlabel('$y^+$')
    plt.ylabel('$\\rho^+$')
    plt.tight_layout()
    plt.savefig(f'{saving_path}/density_plus.png', dpi=300) 
    plt.close() 

    # Generate plots  
    box.plot_van_driest(van_driest, x_, x_str,  testing_path=pino_path, 
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
