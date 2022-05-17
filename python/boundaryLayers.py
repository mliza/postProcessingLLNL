#!/opt/homebrew/bin/python3
'''
    Date:   01/19/2022
    Author: Martin E. Liza
    File:   boundaryLayers.py
    Def:    fluids properties on SI

    Author		    Date		Revision
    --------------------------------------------
    Martin E. Liza	01/19/2022	Initial version.
'''
import pandas as pd 
import numpy as np
import pyvista as pv 
import matplotlib.pyplot as plt 
import os 
import IPython 

# Paths 
script_path  = os.getcwd() 
project_path = os.path.dirname(os.getcwd())
data_path    = os.path.join(project_path, 'plate_data', 'visit_data') 
temp_path    = os.path.join(data_path, 'temperature')
vel_path     = os.path.join(data_path, 'velocityX')
pvtu_path    = '/Users/martin/Desktop/workingFiles/plate_2/SOLUT.0757000.pvtu'

# Uses py-vista functionality 
def data_cutter(pvtu_path):
    pvtu_in        = pv.read(pvtu_path) 
    data_names     = pvtu_in.array_names 
    cell_centers   = pvtu_in.cell_centers()
    cell_positions = np.array(cell_centers.points) 
    orthogonal_slice = pvtu_in.slice_orthogonal(z=0.00085) 
    # p = pv.Plotter() 
    # p.add_mesh(orthogonal_slice)  
    # p.show_grid() 
    # p.show()  

# Visit data Loader 
def data_loader(data_path_in, x_value, y_value):
    curve_list = sorted(os.listdir(data_path_in)) 
    data_dict  = { }
    for i, value in enumerate(curve_list): 
        probe_name = value.split('.')[0] 
        data_dict[probe_name] = { } 
        probe_file = os.path.join(data_path_in, value) 
        data_in    = pd.read_csv(probe_file, sep=" ", header=None)
        data_dict[probe_name][f'{x_value}'] = np.array(data_in.iloc[:][0])
        data_dict[probe_name][f'{y_value}'] = np.array(data_in.iloc[:][1])
    return data_dict 

# Boundary Layer data Parser 
def boundary_layer_thickness(data_dict, variable, freestream_condition):
    probe_numbers = list(data_dict.keys())
    probe_x_axis  = np.arange(9, 16) * 1E-2 
    edge_variable = [ ]
    thickness_y   = [ ]
    wall_variable = [ ]
    # Plot all data and find plates where BL = 0.99U_freestream
    plt_1 = plt.figure(f'{variable}_1') 
    for i in probe_numbers:
        abs_min = np.abs(data_dict[i][f"{variable}"] / 
                         (0.99 * freestream_condition))  
        indx = np.where(abs_min == np.min(abs_min))[0][0]
        thickness_y.append(data_dict[i]['Y'][indx])
        edge_variable.append(data_dict[i][f'{variable}'][indx])
        # Forward difference for tau_wall 
        '''
        tau_wall = ( (data_dict[i]['Ux'][1] - data_dict[i]['Ux'][0]) / 
                     (data_dict[i]['Y'][1] - data_dict[i]['Y'][0]) )
        data_dict[i]['tau_wall'] = tau_wall 
        '''
        wall_variable.append(data_dict[i][f'{variable}'][0])
        plt.plot(data_dict[i][f'{variable}'], data_dict[i]['Y'] * 1E3, 
                '-', label=f'{i}', linewidth=2)
    # Plot figure 1
    plt.legend() 
    plt.grid('-.')
    if (variable == 'Ux'):
        plt.xlabel('Ux [m/s]')
    if (variable == 'T'):
        plt.xlabel('T [K]')
    plt.ylabel('Thickness [mm]')
    # Plot figure 2 
    plt_2 = plt.figure(f'{variable}_2') 
    plt.plot(probe_x_axis * 1E3, np.array(thickness_y) * 1E3, 
            '-*', linewidth=2) 
    plt.grid('-.') 
    plt.xlabel('Streamwise [mm]')
    plt.ylabel('Thickness [mm]')
    #plt.show() 
    # Return dictionary 
    IPython.embed(colors='Linux') 
    dict_out = { 'bl_thickness_y'   : thickness_y, 
                 f'edge_{variable}' : edge_variable, 
                 f'wall_{variable}' : wall_variable,
                 'wall_position_x'  : probe_x_axis, 
                 'wall_probe_data'  : data_dict }
    return dict_out  

# Sutherland's Law
def sutherland_law(temperature): 
# https://doc.comsol.com/5.5/doc/com.comsol.help.cfd/cfd_ug_fluidflow_high_mach.08.27.html 
    # From table for air 
    mu_ref  = 1.716e-5 #[kg/ms]
    T_ref   = 273      #[K]
    S_const = 111      #[K] 
    mu = mu_ref * ( (temperature / T_ref)**(3/2) * ( 
                    (T_ref + S_const) / (temperature + S_const) ) ) 
    return mu 

# Calculate Reynolds number 
def reynolds_number(velocity, density, dynamic_viscosity, 
                    characteristic_length=1): 
    reynolds_number = ( (density * velocity * characteristic_length) / 
                            dynamic_viscosity )
    return reynolds_number

# Calculate Speed of Sound
def speed_of_sound(temperature):
    gamma            = 1.4   #[ ]
    gas_constant_air = 286.9 #[J/kg] 
    speed_of_sound = np.sqrt(gamma * gas_constant_air * temperature)
    return speed_of_sound 

#def wall_law(temperature_dict, velocity_dict):

# Boundary Layer Thickness approximation for a hypersonics BL
def boundary_layer_hypersonics(mach_number, reynolds_number): 
    # Modern Compressible Flow With Historical Perspective, 4th by John Anderson
    # Section 15.2.3 Viscous Interaction
    boundary_layer_thickness = mach_number**2 / np.sqrt(reynolds_number) 
    return boundary_layer_thickness


if __name__ == "__main__":
    velocity          = 3000.0  #[m/s]
    velocity_2        = 2600.0  #[m/s] 
    density           = 0.18874 #[kg/m3]
    pressure          = 11737   #[Pa] 
    temperature       = 216.66  #[K] 
    temperature_2     = 2164    #[K] 
    char_length       = 158e-3  #[m], total plate length  
    # Calculations 
    dynamic_viscosity = sutherland_law(temperature)  
    reynolds_number   = reynolds_number(velocity, density, dynamic_viscosity,
                        characteristic_length=char_length) 
    mach_number       = velocity / speed_of_sound(temperature)
    bl_thickness      = boundary_layer_hypersonics(mach_number, reynolds_number) 
    temperature_dict  = data_loader(temp_path, 'Y', 'T')
    velocity_dict     = data_loader(vel_path, 'Y', 'Ux')
    vel_bl            = boundary_layer_thickness(velocity_dict, 'Ux', velocity_2)
    temp_bl           = boundary_layer_thickness(temperature_dict, 'T', temperature_2)
    plt.show() 


                        


