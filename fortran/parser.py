#!/opt/homebrew/bin/python3
'''
    Date:   04/09/2022
    Author: Martin E. Liza
    File:   parser.py
    Def:

    Author		    Date		Revision 
    ----------------------------------------------------
    Martin E. Liza	04/09/2022	Initial version.
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

# Data split 
def data_split(dict_in, nx, ny, nz):
    N = range(nx* ny* nz) 
    mapping  = data_mapping(path_in) 
    dict_out = { } 
    # Using zip to iterates through 2 loops at the same time 
    for (i, j, k, n) in zip(mapping[0], mapping[1], mapping[2], N):
        # Initialize arrays only at the beginning 
        if n == 0:
            for key in dict_in: 
                dict_out[key] = np.empty([nx, ny, nz]) 
        for key in dict_in: 
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
def plane_xz(data_in, nx, nz, y_val):
    X, Y = np.meshgrid(data_in['X'][0,:,0], data_in['Y'][0,:,0]) 
    IPython.embed(colors='Linux') 
    ax = plt.axes(projection='3d') 
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')
    #https://matplotlib.org/3.5.0/gallery/mplot3d/surface3d.html 
    for k in range(nz): 
        ax.plot3D(data_in['X'][:,y_val,k], data_in['Y'][:,y_val,k], 
                 data_in['Z'][:,y_val,k], 'o')



if __name__ =="__main__":
    path_in = '../dataOut'
    writing_flag = False
    nx     = 1359  
    ny     = 89    
    nz     = 638   
    var_in = ['X', 'Y', 'Z', 'Ux', 'Uy', 'Uz', 'P', 'T', 'DIL', 'RHO', 'GRADRHOMAG'] 
    
    # Writing flag 
    if writing_flag is True: 
        dict_in = { }
        for i in var_in:
            dict_in[i] = data_loader(i, path_in) 
        dict_out = data_split(dict_in, nx, ny, nz)
        pickle_manager(pickle_name='data', pickle_path=path_in, data_in=dict_out)  

    # Loading flag 
    if writing_flag is False:
        data_in = pickle_manager(pickle_name='data', pickle_path=path_in)  
        plane_xz(data_in, nx, nz, y_val=80) 



