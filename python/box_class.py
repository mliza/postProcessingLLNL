#!/usr/local/bin/python3
'''
    Date:   07/25/2022
    Author: Martin E. Liza
    File:   box_class.py 
    Def:    Functions used to post process binary 
            box probes output by margot.

    Author           Date         Revision
    -----------------------------------------------------
    Martin E. Liza   07/19/2022   Initial Version.
    Martin E. Liza   07/25/2022   Added mean fields and 
                                  Reynolds decomposition.
'''
import numpy as np 
import matplotlib.pyplot as plt
import IPython
import sys 
import os 
from dataclasses import dataclass, field  
from scipy.io import FortranFile
from scipy.optimize import curve_fit 
# My own stuff 
scripts_path   = os.environ.get('SCRIPTS')
python_scripts = os.path.join(scripts_path, 'Python') 
sys.path.append(python_scripts) 
import helper_class as helper 

# Probe Class 
@dataclass 
class Box():
# Initialize variables 
    nx: int
    ny: int
    nz: int
    n_total : int = field(init=False)

 # Initialize variables 
    def __post_init__(self):
        self.n_total = self.nx * self.ny * self.nz 
    
# Loads mapping fortran data and saves a pickle file if pickle_path given 
    def mapping_reader(self, mapping_data_in, pickle_path=None):
        mapping_out = np.empty([self.n_total, 4], dtype=int)
        f_in        = FortranFile(mapping_data_in, 'r') 
        for n in range(self.n_total):
            mapping_out[n] = f_in.read_ints() 
        f_in.close() 
        # Save as a pickle file 
        if pickle_path is None:
            return mapping_out
        else:
            helper_scripts = helper.Helper() 
            helper_scripts.pickle_manager(pickle_name_file='mapping', 
                                          pickle_path=pickle_path, 
                                          data_to_save=mapping_out)

# Split the 1D array into a 3D array 
    def split_plot3D(self, array_1D, mapping):
        array_3D = np.empty([self.nx, self.ny, self.nz])
        for n in range(self.n_total):
            # Mapping = [n, i, j, k] 
            i = mapping[n][0]
            j = mapping[n][1]
            k = mapping[n][2]
            array_3D[i,j,k] = array_1D[n]  
        return array_3D

# NOTE: Rename in C++ SHEAR for GRADU
# Return a dictionary with gradient fields 
    def gradient_fields(self, array_dict_1D):  
        omega_x    = 1/2 * (array_dict_1D['GRADV_23'] - 
                            array_dict_1D['GRADV_32'])  
        omega_y    = 1/2 * (array_dict_1D['GRADV_31'] - 
                            array_dict_1D['GRADV_13']) 
        omega_z    = 1/2 * (array_dict_1D['GRADV_12'] - 
                            array_dict_1D['GRADV_21']) 
        u_x        = array_dict_1D['Ux'] 
        u_y        = array_dict_1D['Uy'] 
        u_z        = array_dict_1D['Uz'] 
        vort_mag   = np.sqrt(omega_x**2 + omega_y**2 + omega_z**2)
        u_mag      = np.sqrt(u_x**2 + u_y**2 + u_z**2) 
        dilatation = (array_dict_1D['GRADV_11'] + 
                      array_dict_1D['GRADV_22'] + 
                      array_dict_1D['GRADV_33']) 
        enstrophy = 2 * vort_mag * vort_mag 

        # Equation from donzis (missing mu multiply results by mu)
        #disipation_solenoidal  = vort_mag**2 
        #dissipation_dilatation = 4/3 * dilatation  
        gradient_dict = { 'VortX'     : omega_x, 
                          'VortY'     : omega_y, 
                          'VortZ'     : omega_z,  
                          'VORTMAG'   : vort_mag,  
                          'DIL'       : dilatation, 
                          'ENSTROPHY' : enstrophy, 
                          'UMAG'      : u_mag}
        return gradient_dict 

# Return fluctuation fields 
    def reynolds_decomposition(self, array_field1D):
        fluctuation = array_field1D - np.mean(array_field1D)
        return fluctuation

# Boundary layer thickness
    def boundary_layer_thickness(self, array_field3D, array_height3D,
                                freestream_value):
        boundary_plane_dict = { }
        temp_field  = np.empty([self.nx, self.nz]) 
        temp_height = np.empty([self.nx, self.nz]) 
        cut_value   = 0.99 * freestream_value 
        # Find positions at 0.99 freestream 
        for i in range(self.nx):
            for k in range(self.nz): 
                indx             = np.abs(array_field3D[i,:,k] - 
                                          cut_value).argmin() 
                temp_field[i,k]  = array_field3D[i,:,k][indx] 
                temp_height[i,k] = array_height3D[i,:,k][indx] 
        # Calculate the mean (compress on z, assume frozen flow)  
        temp_mean_height = np.empty(self.nx)
        for i in range(self.nx):
            temp_mean_height[i] = np.mean(temp_height[i,:])
        # Dictionary to return 
        boundary_plane_dict['field']          = temp_field 
        boundary_plane_dict['thickness']      = temp_height
        boundary_plane_dict['mean_thickness'] = temp_mean_height
        boundary_plane_dict['mean_field']     = self.mean_fields(array_field3D)['mean_y'] 

        return boundary_plane_dict 

# Calculates fluctuation fields in a given 3D data set, 
# assume frozen flow hypothesis on z, and returns a 
# 2D field as a function of x and y.  
    def mean_fields(self, array_3D):
        mean_yz = np.empty([self.ny, self.nz]) 
        mean_xz = np.empty([self.nx, self.nz]) 
        mean_xy = np.empty([self.nx, self.ny]) 
        mean_x  = np.empty(self.nx)
        mean_y  = np.empty(self.ny)
        mean_z  = np.empty(self.nz)
        # Array y-z
        for k in range(self.nz):
            for j in range(self.ny):
                mean_yz[j,k] = np.mean(array_3D[:,j,k]) 
        # Array x-z
        for k in range(self.nz):
            for i in range(self.nx):
                mean_xz[i,k] = np.mean(array_3D[i,:,k]) 
        # Array x-y
        for j in range(self.ny):
            for i in range(self.nx):
                mean_xy[i,j] = np.mean(array_3D[i,j,:]) 
        # Mean x 
        for i in range(self.nx): 
            mean_x[i] = np.mean(mean_xy[i,:])
        # Mean y 
        for j in range(self.ny):
            mean_y[j] = np.mean(mean_yz[j,:])
        # Mean z 
        for k in range(self.nz):
            mean_z[k] = np.mean(mean_xz[:,k])

        dict_out = { 'mean_xy' : mean_xy,
                     'mean_yz' : mean_yz,
                     'mean_xz' : mean_xz, 
                     'mean_x'  : mean_x,
                     'mean_y'  : mean_y,
                     'mean_z'  : mean_z }
        return dict_out 

# Reynolds Decomposition 
    def reynolds_decomposition(self, array_1D):
        decomposition_1D = array_1D - np.mean(array_1D)
        return decomposition_1D 

# Reynolds stress structure parameters
    def reynolds_stress_structure_parameters(self, fluctuation_1D):
        u_x     = fluctuation_1D['Ux']
        u_y     = fluctuation_1D['Uy']
        u_z     = fluctuation_1D['Uz']
        k       = fluctuation_1D['K'] 
        rssp_1D = (2 * (u_x * u_y + u_y * u_z + u_x * u_z) + k) / k
        return rssp_1D
            

# Wall shear-stress 
    def wall_shear_stress(self, velocity_boundary_dict, grid_mean_dict):
        y_mean = grid_mean_dict['mean_y'] 
        u_mean = velocity_boundary_dict['mean_field'] 
        # IPython.embed(colors='Linux') 

# Fitting function
    def smoothing_function(self, data_in, box_pts):
        box         = np.ones(box_pts)/box_pts 
        data_smooth = np.convolve(data_in, box, mode='same')
        return data_smooth 

# Plot boundary Layers 
    def plot_boundary_layers(self, velocity_boundary_dict, 
            temperature_boundary_dict, grid_mean_dict, 
            velocity_freestream, temperature_freestream,
            saving_path=None):
        # Loading variables 
        n_box           = 50  
        temp_mean_thick = temperature_boundary_dict['mean_thickness'] * 10**3
        vel_mean_thick  = velocity_boundary_dict['mean_thickness'] * 10**3
        temp_mean_field = temperature_boundary_dict['mean_field'] 
        vel_mean_field  = velocity_boundary_dict['mean_field'] 
        x_mean          = grid_mean_dict['mean_x'] * 10**2
        y_mean          = grid_mean_dict['mean_y'] * 10**3
        v_color         = 'mediumturquoise' 
        t_color         = 'darkorange'
        temp_smooth     = self.smoothing_function(temp_mean_thick, n_box) 
        vel_smooth      = self.smoothing_function(vel_mean_thick, n_box) 
        n_box           /= 2
        n_box           = int(n_box) 
        # Plotting figures  
        fig, (ax1, ax2) = plt.subplots(1,2, figsize=(8,5))
        # Plot thickness 
        ax1.plot(x_mean, vel_mean_thick, 'o', markersize=3,
                markerfacecolor='lightgrey', markeredgecolor='k') 
        ax1.plot(x_mean[n_box:-n_box], vel_smooth[n_box:-n_box], color=v_color, 
                 linestyle='-', linewidth=1.5, label='Ux')
        ax1.plot(x_mean, temp_mean_thick, 'o', markersize=3, 
                markerfacecolor='lightgrey', markeredgecolor='k')
        ax1.plot(x_mean[n_box:-n_box], temp_smooth[n_box:-n_box], color=t_color, 
                 linestyle='-', linewidth=1.5, label='T') 
        ax1.legend()
        ax1.set_ylabel('y-axis [mm]')
        ax1.set_xlabel('x-axis [cm]')
        ax1.grid('-.') 
        ax1.legend() 
        # Plot value 
        l1 = ax2.plot(vel_mean_field, y_mean, color=v_color, 
                linestyle='-', linewidth=3, label=f'Ux={velocity_freestream:.2f}[m/s]')
        ax2.set_xlabel('Ux [m/s]', color=v_color)
        ax21 = ax2.twiny() 
        l2 = ax21.plot(temp_mean_field, y_mean, color=t_color, 
                linestyle='-', linewidth=3, label=f'T={temperature_freestream:.2f}[K]')
        ax21.set_xlabel('T [K]', color=t_color)
        ax2.set_ylabel('y-axis [mm]')
        ax2.grid('-.') 
        ax2.legend(handles=l1+l2, loc='upper center') 
        
        if saving_path == None:
            plt.show() 
        if saving_path != None:
            fig.tight_layout()
            fig.savefig(f'{saving_path}/boundary_layers.png', dpi=300)
            plt.close() 

# Plot boundary surface 
    def plot_boundary_surface(self, boundary_plane_dict, 
                              grid_mean_dict):  
        # Loading data 
        X_mean        = grid_mean_dict['mean_x'] 
        Z_mean        = grid_mean_dict['mean_z'] 
        height        = boundary_plane_dict['thickness']
        boundary_mean = boundary_plane_dict['mean_thickness']  
        X,Y           = np.meshgrid(Z_mean, X_mean)

        #IPython.embed(colors='Linux') 
        #fig = plt.figure(figsize=(8,8))
        fig = plt.figure()
        ax  = fig.add_subplot(111, projection='3d') 
        ax.plot_surface(X, Y, height) 
            
# Plot line for 2 variables  
    def plot_lineXY(self, array_dict_3D, var_x, var_y, 
                    x_dim=None, y_dim=None, z_dim=None, saving_path=None):
        if x_dim is None:
            x_axis    = array_dict_3D[var_x][:, y_dim, z_dim] 
            y_axis    = array_dict_3D[var_y][:, y_dim, z_dim] 
            label_str = f'y={y_dim}, z={y_dim}'
        if y_dim is None:
            x_axis    = array_dict_3D[var_x][x_dim, :, z_dim] 
            y_axis    = array_dict_3D[var_y][x_dim, :, z_dim] 
            label_str = f'x={x_dim}, z={y_dim}'
        if z_dim is None:
            x_axis    = array_dict_3D[var_x][x_dim, y_dim, :] 
            y_axis    = array_dict_3D[var_y][x_dim, y_dim, :]
            label_str = f'x={x_dim}, z={y_dim}'
        # Legend, title 
        plt.plot(x_axis, y_axis, color='k',  linestyle='-', linewidth=3,
                 label=label_str) 
        plt.grid('-.') 
        plt.legend() 
        plt.xlabel(f'{var_x}')
        plt.ylabel(f'{var_y}')
        plt.title(f'{var_y} vs. {var_x}') 
        # Saving if needed 
        if saving_path == None:
            plt.show() 
        if saving_path != None:
            plt.savefig(f'{saving_path}/{var_x}_{var_y}.png', dpi=300)
            plt.close() 

# Making contour plots 
    def plot_contour(self, array_dict_3D, grid_x, grid_y, 
                field, slice_cut, slice_direction,
                levels=6, cmap='inferno', saving_path=None): 
        # Create slides 
        slice_direction = slice_direction.upper() 
        if slice_direction == 'X':
            x_plane     = array_dict_3D[grid_x][slice_cut,:,:]
            y_plane     = array_dict_3D[grid_y][slice_cut,:,:]
            z_plane     = array_dict_3D[field][slice_cut,:,:]
            slice_value = array_dict_3D[slice_direction][:,-1,-1][slice_cut] 
        if slice_direction == 'Y':
            x_plane     = array_dict_3D[grid_x][:,slice_cut,:]
            y_plane     = array_dict_3D[grid_y][:,slice_cut,:]
            z_plane     = array_dict_3D[field][:,slice_cut,:]
            slice_value = array_dict_3D[slice_direction][-1,:,-1][slice_cut] 
        if slice_direction == 'Z':
            x_plane     = array_dict_3D[grid_x][:,:,slice_cut]
            y_plane     = array_dict_3D[grid_y][:,:,slice_cut]
            z_plane     = array_dict_3D[field][:,:,slice_cut]
            slice_value = array_dict_3D[slice_direction][-1,-1,:][slice_cut] 
        # Plotting 
        plt.contourf(x_plane, y_plane, z_plane, 
                    levels=levels, cmap=cmap)

        plt.xlabel(f'{grid_x} [m]')
        plt.ylabel(f'{grid_y} [m]')
        plt.title(f'{field}, at {slice_direction}={slice_value:.3E} [m]') 
        plt.colorbar() 
        # Saving if needed 
        if saving_path == None:
            plt.show() 
        if saving_path != None:
            plt.savefig(f'{saving_path}/contour{grid_x}{grid_y}_{field}.png', 
                            bbox_inches='tight', dpi=300)
            plt.close() 
