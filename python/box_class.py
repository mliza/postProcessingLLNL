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
        omega_x   = 1/2 * (array_dict_1D['GRADV_23'] - array_dict_1D['GRADV_32'])  
        omega_y   = 1/2 * (array_dict_1D['GRADV_31'] - array_dict_1D['GRADV_13']) 
        omega_z   = 1/2 * (array_dict_1D['GRADV_12'] - array_dict_1D['GRADV_21']) 
        vort_mag  = np.sqrt(omega_x**2 + omega_y**2 + omega_z**2)
        dilatation =  (array_dict_1D['GRADV_11'] + 
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
                          'ENSTROPHY' : enstrophy }
        return gradient_dict 

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
            
# Plot line for 2 variables  
    def plot_lineXY(self, array_dict_3D, var_x, var_y, 
                    x_dim=None, y_dim=None, z_dim=None, saving_path=None):
        if x_dim is None:
            plt.plot(array_dict_3D[var_x][:, y_dim, z_dim], 
                    array_dict_3D[var_y][:, y_dim, z_dim], '-o', 
                    label=f'y={y_dim}, z={z_dim}')
        if y_dim is None:
            plt.plot(array_dict_3D[var_x][x_dim, :, z_dim], 
                    array_dict_3D[var_y][x_dim, :, z_dim], '-o',
                    label=f'x={x_dim}, z={z_dim}')
        if z_dim is None:
            plt.plot(array_dict_3D[var_x][x_dim, y_dim, :], 
                    array_dict_3D[var_y][x_dim, y_dim, :], '-o',
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
            plt.savefig(f'{saving_path}/contour{grid_x}{grid_y}_{field}.png', bbox_inches='tight', dpi=300)
            plt.close() 
