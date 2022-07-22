#!/opt/homebrew/bin/python3 
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
import IPython 
import pickle 
import matplotlib  
import matplotlib.pyplot as plt 
from scipy.io import FortranFile
import time 
# Mine 
import sys 
import os 
scripts_path   = os.environ.get('SCRIPTS')
python_scripts = os.path.join(scripts_path, 'Python')
sys.path.insert(1, python_scripts) 
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
def data_split(dict_in, nx, ny, nz, mapping):
    N = nx* ny* nz 
    dict_out = { } 
    for key in dict_in: 
        print(f'Working on {key}:')
        start = time.time() 
        for n in range(N): 
            if n == 0:
                dict_out[key] = np.empty([nx, ny, nz]) 
            i = mapping[n][0] 
            j = mapping[n][1] 
            k = mapping[n][2] 
            dict_out[key][i][j][k] = dict_in[key][n] 
        end = time.time() 
        elapse_time = end - start
        print(f'Ending {key}, it took {elapse_time:.4}s\n')
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
        plt.close() 


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
    path_in      = '../../plate_data/data_11'
    path_temp    = os.path.join(path_in, 'temp_data')
    path_pickle  = os.path.join(path_in, 'pickle')
    saving_path  = os.path.join(path_in, 'results') 
    writing_flag = False   
    mapping_flag = False  
    nx     = 1439 
    ny     = 85  
    nz     = 638 
    var_in = [ 'X', 'Y', 'Z', 
               'Ux', 'Uy', 'Uz',  
               'T', 'RHO', 'P',  
               'RHOE', 'GRADRHOMAG',
               'GRADV_11', 'GRADV_12', 'GRADV_13', 
               'GRADV_21', 'GRADV_22', 'GRADV_23', 
               'GRADV_31', 'GRADV_32', 'GRADV_33'] 
    
    # Writing flag 
    helper = helper.Helper() 
    if writing_flag: 
        if mapping_flag: 
            print('Starting Mapping:')
            start = time.time() 
            mapping = mapping_reader(path_temp, nx*ny*nz) 
            end = time.time() 
            elapse_time = end - start
            print(f'Ending mapping, it took {elapse_time:.4}s\n')
            helper.pickle_manager(pickle_name_file='mapping', pickle_path=path_pickle,
                                 data_to_save=mapping)
        mapping = helper.pickle_manager(pickle_name_file='mapping', 
                    pickle_path=path_pickle) 
        dict_in = { }
        for i in var_in:
            dict_in[i] = helper.fortran_data_loader(f'{i}.dat', path_temp) 
        helper.pickle_manager(pickle_name_file='data_in', pickle_path=path_pickle,
                        data_to_save=dict_in)  
        dict_out = data_split(dict_in, nx, ny, nz, mapping)
        helper.pickle_manager(pickle_name_file='data_out', pickle_path=path_pickle,
                        data_to_save=dict_out)  

    # Loading flag 
    if not writing_flag:
        data_in = helper.pickle_manager(pickle_name_file='data_in', 
                                 pickle_path=path_pickle)  
        data_out = helper.pickle_manager(pickle_name_file='data_out', 
                                 pickle_path=path_pickle)  
        P_DIL = data_in['P'] * (data_in['GRADV_11'] + data_in['GRADV_22'] + data_in['GRADV_33']) 
        #helper.pickle_dict_add(P_DIL, 'P_DIL', path_pickle, 'data_in', 'data_in_test') 
        #TEST = helper.pickle_manager(pickle_name_file='data_in_test',  
        #                         data_to_save=path_pickle)  


        IPython.embed(colors='Linux')
        contour(data_out, grid_x='X', grid_y='Y', field='T',
                slice_cut=30, slice_direction='Z', levels=500, cmap='bwr') 
        plot_testing(data_out, nx, ny, nz, saving_path, val_fix=20)
        plot_lineXY(data_out, 'X', 'T', y_dim=2, z_dim=8, saving_path=saving_path)
        plot_lineXY(data_out, 'T', 'Y', x_dim=2, z_dim=8, saving_path=saving_path)
        plot_lineXY(data_out, 'X', 'Ux', y_dim=2, z_dim=8, saving_path=saving_path)
        plot_lineXY(data_out, 'Ux', 'Y', x_dim=2, z_dim=8, saving_path=saving_path)
        plot_lineXY(data_out, 'Ux', 'Y', x_dim=-1, z_dim=-1, saving_path=saving_path)
        plot_lineXY(data_out, 'T', 'Y', x_dim=-1, z_dim=-1, saving_path=saving_path)
        #%matplotlib auto
        plt.ion() 
        #plane_xz(data_in, nx, nz, y_val=80) 

