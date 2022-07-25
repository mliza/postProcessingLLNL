#!/usr/local/bin/python3
'''
    Date:   06/24/2021
    Author: Martin E. Liza
    File:   box_class.py Def:               

    Author           Date         Revision
    ----------------------------------------------------------------------- 
    Martin E. Liza   07/19/2022   Initial Version.
'''
import numpy as np 
import matplotlib.pyplot as plt
import IPython
import sys 
import os 
from dataclasses import dataclass, field  
from scipy.io import FortranFile
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
    def gradient_fields(self, data_dict):  
        omega_x   = 1/2 * (data_dict['GRADV_23'] - data_dict['GRADV_32'])  
        omega_y   = 1/2 * (data_dict['GRADV_31'] - data_dict['GRADV_13']) 
        omega_z   = 1/2 * (data_dict['GRADV_12'] - data_dict['GRADV_21']) 
        vort_mag  = np.sqrt(omega_x**2 + omega_y**2 + omega_z**2)
        dilatation =  (data_dict['GRADV_11'] + 
                      data_dict['GRADV_22'] + 
                      data_dict['GRADV_33']) 
        enstrophy = 2 * vort_mag * vort_mag 

        # Equation from donzis (missing mu multiply results by mu)
        #disipation_solenoidal  = vort_mag**2 
        #dissipation_dilatation = 4/3 * dilatation  
        gradient_dict = { 'VortX'     : omega_x, 
                          'Vorty'     : omega_y, 
                          'Vortz'     : omega_z,  
                          'VORTMAG'   : vort_mag,  
                          'DIL'       : dilatation, 
                          'ENSTROPHY' : enstrophy }
        return gradient_dict 

# Calculate mean fields in a given direction
    def mean_fields(self, data_in, given_direction):
        if given_direction == 'x':
            mean_field = np.empty([self.nx, self.ny]) 
            for i in range(self.nx):
                for j in range(self.ny):
                    mean_field[i,j] = np.mean(data_in[i,j,:])
            
# Plot line for 2 variables  

    def plot_lineXY(self, data_in, var_x, var_y, x_dim=None, y_dim=None, z_dim=None,
                    saving_path=None):
        if x_dim is None:
            plt.plot(data_in[var_x][:, y_dim, z_dim], 
                    data_in[var_y][:, y_dim, z_dim], '-o', 
                    label=f'y={y_dim}, z={z_dim}')
        if y_dim is None:
            plt.plot(data_in[var_x][x_dim, :, z_dim], 
                    data_in[var_y][x_dim, :, z_dim], '-o',
                    label=f'x={x_dim}, z={z_dim}')
        if z_dim is None:
            plt.plot(data_in[var_x][x_dim, y_dim, :], 
                    data_in[var_y][x_dim, y_dim, :], '-o',
                    label=f'x={x_dim}, y={y_dim}')
        # Legend, title 
        plt.grid('-.') 
        plt.legend() 
        plt.xlabel(f'{var_x}')
        plt.ylabel(f'{var_y}')
        plt.title(f'{var_y} vs. {var_x}') 
        # Saving if needed 
        if saving_path == None:
            plt.show() 
        if saving_path != None:
            plt.savefig(f'{saving_path}/{var_x}_{var_y}.png')
            plt.close() 

# Making contour plots 
    def plot_contour(self, data_in, grid_x, grid_y, 
                field, slice_cut, slice_direction,
                levels=6, cmap='inferno', saving_path=None): 
        # Create slides 
        slice_direction = slice_direction.upper() 
        if slice_direction == 'X':
            x_plane     = data_in[grid_x][slice_cut,:,:]
            y_plane     = data_in[grid_y][slice_cut,:,:]
            z_plane     = data_in[field][slice_cut,:,:]
            slice_value = data_in[slice_direction][:,-1,-1][slice_cut] 
        if slice_direction == 'Y':
            x_plane     = data_in[grid_x][:,slice_cut,:]
            y_plane     = data_in[grid_y][:,slice_cut,:]
            z_plane     = data_in[field][:,slice_cut,:]
            slice_value = data_in[slice_direction][-1,:,-1][slice_cut] 
        if slice_direction == 'Z':
            x_plane     = data_in[grid_x][:,:,slice_cut]
            y_plane     = data_in[grid_y][:,:,slice_cut]
            z_plane     = data_in[field][:,:,slice_cut]
            slice_value = data_in[slice_direction][-1,-1,:][slice_cut] 
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
            plt.savefig(f'{saving_path}/contour{grid_x}{grid_y}_{field}.png', bbox_inches='tight', dpi=300)
            plt.close() 





# Reynolds Decomposition 
    def reynolds_decomposition(self, data_dict): 
        u_x = data_dict['Ux']
        u_y = data_dict['Uy']
        u_z = data_dict['Uz']





