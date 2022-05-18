#!/opt/homebrew/bin/python3
'''
    Date:   04/09/2022
    Author: Martin E. Liza
    File:   parser.py
    Def:

    Author		    Date		Revision 
    ----------------------------------------------------
    Martin E. Liza	04/09/2022	Initial version.
    Martin E. Liza  05/18/2022  Added the data_split function
                                easier to implement in fortran and,
                                removed the previous which used zip.
'''
import numpy as np
import os 
import IPython 
import pickle 
import matplotlib  
import matplotlib.pyplot as plt 
from scipy.io import FortranFile

# Loading data 
def data_loader(variable_in, abs_path_in): 
    data_in = os.path.join(abs_path_in, f'{variable_in}.dat')
    f_in = FortranFile(data_in, 'r')
    data = f_in.read_reals(dtype=np.float64)
    f_in.close() 
    return data

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
    for n in range(N): 
        if n == 0:
            for key in dict_in: 
                dict_out[key] = np.empty([nx, ny, nz]) 
        for key in dict_in: 
            i = mapping[n][0] 
            j = mapping[n][1] 
            k = mapping[n][2] 
            dict_out[key][i][j][k] = dict_in[key][n] 
    return dict_out

# Save data as pickle and loads data as a pickle 
def pickle_manager(pickle_name, pickle_path, data_in=None):
    # Loads pickle file 
    if data_in is None:  
        file_in   = os.path.join(f'{pickle_path}',f'{pickle_name}.pickle') 
        pickle_in = open(file_in, 'rb')
        return pickle.load(pickle_in)
    # Creates pickle file  
    else:
        file_out   = os.path.join(f'{pickle_path}',f'{pickle_name}.pickle') 
        pickle_out = open(file_out, 'wb') 
        pickle.dump(data_in, pickle_out)
        pickle_out.close() 

# Plot plane 
def plot_line(data_in, val_x, val_y, x_dim=None, y_dim=None, z_dim=None):
    if x_dim is None:
        plt.plot(data_in[val_x][:, y_dim, z_dim], 
                data_in[val_y][:, y_dim, z_dim], '-o')
    if y_dim is None:
        plt.plot(data_in[val_x][x_dim, :, z_dim], 
                data_in[val_y][x_dim, :, z_dim], '-o')
    if z_dim is None:
        plt.plot(data_in[val_x][x_dim, y_dim, :], 
                data_in[val_y][x_dim, y_dim, :], '-o')
    plt.grid('-.') 

def plane_xz(data_in, nx, nz, y_val):
    X, Y = np.meshgrid(data_in['X'][0,:,0], data_in['Y'][0,:,0]) 
    ax = plt.axes(projection='3d') 
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')
    #https://matplotlib.org/3.5.0/gallery/mplot3d/surface3d.html 
    for k in range(nz): 
        ax.plot3D(data_in['X'][:,y_val,k], data_in['Y'][:,y_val,k], 
                 data_in['Z'][:,y_val,k], 'o')

if __name__ =="__main__":
    path_in     = '../../plate_data/data_9'
    path_temp   = os.path.join(path_in, 'temp_data')
    path_pickle = os.path.join(path_in, 'pickle')
    writing_flag = False 
    nx     = 1359  
    ny     = 89    
    nz     = 638   
    var_in = ['X', 'Y', 'Z', 'Ux', 'Uy', 'Uz', 'P', 'T', 'DIL', 'RHO', 'GRADRHOMAG'] 
    
    # Writing flag 
    if writing_flag is True: 
        dict_in = { }
        for i in var_in:
            dict_in[i] = data_loader(i, path_temp) 
        pickle_manager(pickle_name='data_in', pickle_path=path_pickle,
                        data_in=dict_in)  
        dict_out = data_split(dict_in, nx, ny, nz, mapping_path=path_temp)
        pickle_manager(pickle_name='data_out', pickle_path=path_pickle,
                        data_in=dict_out)  

    # Loading flag 
    if writing_flag is False:
        data_in = pickle_manager(pickle_name='data_in', 
                                 pickle_path=path_pickle)  
        data_out = pickle_manager(pickle_name='data_out', 
                                 pickle_path=path_pickle)  
        IPython.embed(colors='Linux') 
        #%matplotlib auto
        plt.ion() 
        plane_xz(data_in, nx, nz, y_val=80) 

