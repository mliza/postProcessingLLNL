#!/usr/bin/python3
'''
    Date:   04/09/2022
    Author: Martin E. Liza
    File:   box_parser.py
    Def:

    Author		    Date		Revision 
    ----------------------------------------------------
    Martin E. Liza	04/09/2022	Initial version.
    Martin E. Liza  05/18/2022  Added the data_split function
                                easier to implement in fortran and,
                                removed the previous which used zip.
    Martin E. Liza  05/20/2022  Added contour plots and some testings plots. 
'''
import numpy as np
import os 
import IPython 
import pickle 
import matplotlib  
import matplotlib.pyplot as plt 
from scipy.io import FortranFile
# Mine 
import helper_class as helper 

# Loading data mapping
def data_mapping(abs_path_in):
    data_in = os.path.join(abs_path_in, 'mappingVector.dat')
    f_in = FortranFile(data_in, 'r')
    x_mapping = f_in.read_ints(dtype=np.int32) 
    y_mapping = f_in.read_ints(dtype=np.int32)
    z_mapping = f_in.read_ints(dtype=np.int32)
    f_in.close() 
    return [x_mapping, y_mapping, z_mapping] 

# Data split (easier to implement in fortran)  
def data_split(dict_in, nx, ny, nz, mapping_path):
    N = nx* ny* nz 
    mapping = data_mapping(mapping_path) 
    mapping = np.array(mapping).transpose() 
    dict_out = { } 
    for key in dict_in: 
        for n in range(N): 
            if n == 0:
                dict_out[key] = np.empty([nx, ny, nz]) 
            i = mapping[n][0] 
            j = mapping[n][1] 
            k = mapping[n][2] 
            dict_out[key][i][j][k] = dict_in[key][n] 
    return dict_out


# Plot line for 2 variables  
def plot_lineXY(data_in, var_x, var_y, x_dim=None, y_dim=None, z_dim=None):
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
    plt.grid('-.') 
    plt.legend() 
    plt.xlabel(f'{var_x}')
    plt.ylabel(f'{var_y}')

# Plot line for 1 variable 
def plot_line(data_in, var_y, x_dim=None, y_dim=None, z_dim=None):
    if x_dim is None:
        plt.plot(data_in[var_y][:, y_dim, z_dim], '-o', 
                label=f'y={y_dim}, z={z_dim}')
    if y_dim is None:
        plt.plot(data_in[var_y][x_dim, :, z_dim], '-o', 
                label=f'x={x_dim}, z={z_dim}')
    if z_dim is None:
        plt.plot(data_in[var_y][x_dim, y_dim, :], '-o', 
                label=f'x={x_dim}, y={y_dim}')
    plt.grid('-.') 
    plt.legend() 
    plt.xlabel('Iterations')
    plt.ylabel(f'{var_y}')

# For Testing
def x_axis_test(x_grid, nx, ny, nz):   
    for i in range(ny * nz - 1):
        plt.plot(x_grid[i * nx:(i + 1) * nx], 'o-')   
    plt.grid('-.')
    plt.xlabel('nx [ ]')
    plt.ylabel('X [m]')

# Making contour plots 
def contour(data_in, grid_x, grid_y, field, slice_cut, slice_direction,
            levels=6, cmap='inferno'): 
    # Create slides 
    slice_direction = slice_direction.upper() 
    if slice_direction == 'X':
        x_plane     = data_out[grid_x][slice_cut,:,:]
        y_plane     = data_out[grid_y][slice_cut,:,:]
        z_plane     = data_out[field][slice_cut,:,:]
        slice_value = data_out[slice_direction][:,-1,-1][slice_cut] 
    if slice_direction == 'Y':
        x_plane     = data_out[grid_x][:,slice_cut,:]
        y_plane     = data_out[grid_y][:,slice_cut,:]
        z_plane     = data_out[field][:,slice_cut,:]
        slice_value = data_out[slice_direction][-1,:,-1][slice_cut] 
    if slice_direction == 'Z':
        x_plane     = data_out[grid_x][:,:,slice_cut]
        y_plane     = data_out[grid_y][:,:,slice_cut]
        z_plane     = data_out[field][:,:,slice_cut]
        slice_value = data_out[slice_direction][-1,-1,:][slice_cut] 

    # Plotting 
    plt.contourf(x_plane, y_plane, z_plane, 
                levels=levels, cmap=cmap)
    plt.xlabel(f'{grid_x} [m]')
    plt.ylabel(f'{grid_y} [m]')
    plt.title(f'{field}, at {slice_direction}={slice_value} [m]') 
    plt.colorbar() 

    
    

if __name__ =="__main__":
    path_in      = '../../plate_data/data_12'
    path_temp    = os.path.join(path_in, 'temp_data')
    path_pickle  = os.path.join(path_in, 'pickle')
    writing_flag = False   
    nx     = 5#1359  #32
    ny     = 3#89    #18
    nz     = 2#638   #16
    var_in = ['X', 'Y', 'Z', 'Ux', 'Uy', 'Uz', 'P', 'T', 'DIL', 'RHO', 'GRADRHOMAG'] 
    var_in = ['X', 'Y', 'Z', 'Ux', 'Uy', 'Uz', 'T'] 

    # TESTING ASCII Vs. BIN
    '''
    ascii_path   = os.path.join(path_in, 'smallBOX_ASCII') 
    f            = open(os.path.join(ascii_path, 'T.xyz'), 'r')
    lines        = f.read().splitlines() 
    f.close() 
    x_ascii      = np.fromstring(lines[2], dtype=np.longdouble, sep=' ') 
    y_ascii      = np.fromstring(lines[3], dtype=np.longdouble, sep=' ') 
    z_ascii      = np.fromstring(lines[4], dtype=np.longdouble, sep=' ') 
    '''
    
    # Writing flag 
    helper = helper.Helper() 
    if writing_flag: 
        dict_in = { }
        for i in var_in:
            dict_in[i] = helper.data_loader(i, path_temp) 
        helper.pickle_manager(pickle_name='data_in', pickle_path=path_pickle,
                        data_in=dict_in)  
        dict_out = data_split(dict_in, nx, ny, nz, mapping_path=path_temp)
        helper.pickle_manager(pickle_name='data_out', pickle_path=path_pickle,
                        data_in=dict_out)  

    # Loading flag 
    if not writing_flag:
        data_in = helper.pickle_manager(pickle_name='data_in', 
                                 pickle_path=path_pickle)  
        data_out = helper.pickle_manager(pickle_name='data_out', 
                                 pickle_path=path_pickle)  
        IPython.embed(colors='Linux') 
        contour(data_out, grid_x='X', grid_y='Z', field='DIL',
                slice_cut=1, slice_direction='y', levels=20, cmap='bwr') 
        #%matplotlib auto
        plt.ion() 
        #plane_xz(data_in, nx, nz, y_val=80) 

