#!/usr/bin/python3
'''
    Date:   04/09/2022
    Author: Martin E. Liza
    File:   box_parser.py
    Def:

    Author		    Date		Revision 
    ------------------------------------------------------------------------
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
def mapping_reader(abs_path_in, N):
    data_in  = os.path.join(abs_path_in, 'mappingVector.dat')
    data_out = np.empty([N, 4], dtype=int) 
    f_in = FortranFile(data_in, 'r')
    for n in range(N): 
        data_out[n] = f_in.read_ints() 
    f_in.close() 
    return data_out  

# Data split (easier to implement in fortran)  
def data_split(dict_in, nx, ny, nz, mapping_path):
    N = nx* ny* nz 
    mapping = mapping_reader(mapping_path, N) 
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
def plot_lineXY(data_in, var_x, var_y, x_dim=None, y_dim=None, z_dim=None,
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
def plot_testing(data_in, nx, ny, nz, saving_path, val_fix=0):
    # Plot Xy, fix z at 0
    for j in range(ny):
        plt.plot(data_in['X'][:,j,val_fix], 'o-') 
    plt.grid('-.')
    plt.xlabel('Iterations')
    plt.ylabel('X-y [m]')
    plt.title(f'X, changing on Y at fixed Z={val_fix}') 
    plt.savefig(f'{saving_path}/Xy_z{val_fix}.png')
    plt.close() 
        
    # Plot Xz, fix y at 0
    for k in range(nz):
        plt.plot(data_in['X'][:,val_fix,k], 'o-') 
    plt.grid('-.')
    plt.xlabel('Iterations')
    plt.ylabel('X-z [m]')
    plt.title(f'X, changing on Z at fixed Y={val_fix}') 
    plt.savefig(f'{saving_path}/Xz_y{val_fix}.png')
    plt.close() 

    # Plot Yx, fix z at 0
    for i in range(nx):
        plt.plot(data_in['Y'][i,:,val_fix], 'o-') 
    plt.grid('-.')
    plt.xlabel('Iterations')
    plt.ylabel('Y-x [m]')
    plt.title(f'Y, changing on X at fixed Z={val_fix}') 
    plt.savefig(f'{saving_path}/Yx_z{val_fix}.png')
    plt.close() 
        
    # Plot Yz, fix x at 0
    for k in range(nz):
        plt.plot(data_in['Y'][val_fix,:,k], 'o-') 
    plt.grid('-.')
    plt.xlabel('Iterations')
    plt.ylabel('Y-z [m]')
    plt.title(f'Y, changing on Z at fixed X={val_fix}') 
    plt.savefig(f'{saving_path}/Yz_x{val_fix}.png')
    plt.close() 

    # Plot Zx, fix y at 0
    for i in range(nx):
        plt.plot(data_in['Z'][i,val_fix,:], 'o-') 
    plt.grid('-.')
    plt.xlabel('Iterations')
    plt.ylabel('Z-x [m]')
    plt.title(f'Z, changing on X at fixed Y={val_fix}') 
    plt.savefig(f'{saving_path}/Zx_y{val_fix}.png')
    plt.close() 
        
    # Plot Zy, fix x at 0
    for j in range(ny):
        plt.plot(data_in['Z'][val_fix,j,:], 'o-') 
    plt.grid('-.')
    plt.xlabel('Iterations')
    plt.ylabel('Z-y [m]')
    plt.title(f'Z, changing on Y at fixed X={val_fix}') 
    plt.savefig(f'{saving_path}/Zy_x{val_fix}.png')
    plt.close() 
        


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
    path_in      = '../../plate_data/data_15'
    path_temp    = os.path.join(path_in, 'temp_data')
    path_pickle  = os.path.join(path_in, 'pickle')
    saving_path  = os.path.join(path_in, 'results') 
    writing_flag = False 
    nx     = 1439 
    ny     = 85  
    nz     = 638 
    var_in = ['X', 'Y', 'Z', 'Ux', 'Uy', 'Uz', 'RHO', 'P', 'T', 'DIL', 
              'GRADRHOMAG', 'RHOE', 'VORTMAG'] 
    
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
        plot_testing(data_out, nx, ny, nz, saving_path, val_fix=-1)
        plot_testing(data_out, nx, ny, nz, saving_path, val_fix=0)
        plot_testing(data_out, nx, ny, nz, saving_path, val_fix=10)
        plot_testing(data_out, nx, ny, nz, saving_path, val_fix=5)
        plot_testing(data_out, nx, ny, nz, saving_path, val_fix=30)
        plot_testing(data_out, nx, ny, nz, saving_path, val_fix=20)

        IPython.embed(colors='Linux') 
        plot_lineXY(data_out, 'X', 'T', y_dim=2, z_dim=8, saving_path=saving_path)
        contour(data_out, grid_x='X', grid_y='Y', field='T',
                slice_cut=30, slice_direction='Z', levels=500, cmap='bwr') 
        #%matplotlib auto
        plt.ion() 
        #plane_xz(data_in, nx, nz, y_val=80) 

