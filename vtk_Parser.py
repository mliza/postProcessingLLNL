#!/usr/local/bin/python3.9
'''
    Date:   10/26/2021
    Author: Martin E. Liza
    File:   vtk_Parser.py
    Def:    clips data at a defined location and returns a dictionary  

    Author		    Date		Revision
    ----------------------------------------------------
    Martin E. Liza	10/26/2021	Initial version.
'''
import numpy as np 
import pyvista as pv
import matplotlib.pyplot as plt 
import pickle
import os 
import IPython 
import time 

# Inputs, create a region where the data will be clipped (approx values) 
workgingFiles_path = '/Users/martin/Desktop/workingFiles/dump_172000'
x1, x2     = 9e-3, 1.2e-2
y1, y2     = 1e-2, 1.2e-2
z1, z2     = 1e-3, 1.2e-3
pvtu_in    = os.path.join(workingFiles_path, 'SOLUT3.172000.pvtu') 
pickle_out = 'clipped_data'
vtk_out    = 'clipped_data'
bounds     = [x1,x2, y1,y2, z1,z2] 
#path_to_save = 'data_out' 
path_to_save = 'text_out' 

def data_to_binary(bounds, pvtu_in, pickle_out, vtk_out): 
# pyvista functions 
    mesh          = pv.read(pvtu_in) # reads the file  
    data_names    = mesh.array_names # array names in mesh  
    cell_center   = mesh.cell_centers() 
    cell_position = np.array(cell_center.points) 
    # Remove CYCLE and TIME from data_names 
    data_names.remove('CYCLE')
    data_names.remove('TIME')
    X = cell_center.points[:,0] 
    Y = cell_center.points[:,1] 
    Z = cell_center.points[:,2] 
    
    IPython.embed(colors='Linux')
    ax = plt.axes(projection='3d') 
    ax.plot3D(X[0:5], Y[0:5], Z[0:5], '*-') 
    #Y = Y.sort() #test  
'''
    # Loops thought and saves data 
    for i in data_names: 
        print(f'Starting {i}\n')
        if (cell_center[i].ndim == 1):
            data_to_save = np.insert(cell_position, 3, 
                                     cell_center[i], axis=1)
            #np.savetxt(os.path.join(path_to_save, f'{i}.out'), data_to_save)
            data_to_save.tofile(os.path.join(path_to_save, f'{i}.bin')) 
        else: 
            for j in range(3):
                positions = ['X', 'Y', 'Z']
                data_to_save = np.insert(cell_position, 3, 
                        cell_center[i][:,j], axis=1)
                #np.savetxt(os.path.join(path_to_save, f'{i}.out'), 
                #        data_to_save)
                data_to_save.tofile(os.path.join(path_to_save, 
                            f'{i}-{positions[j]}.bin')) 
        print(f'Ending {i}\n\n')
'''

'''
    clipped_data = mesh.clip_box(bounds, invert=False)
    clipped_data.save(f'{vtk_out}.vtk') # Save clipped data as a vtk  
    data_dictionary = { }
# Return Dictionary with data 
    for i in data_names:
        data_dictionary[i] = clipped_data[i] 

# Fix dictionary 
# Create x,y,z positions 
    cell_center   = clipped_data.cell_centers() 
    cell_position = cell_center.points.T 
    position_x    = cell_position[0] 
    position_y    = cell_position[1]
    position_z    = cell_position[2] 

# Re populate dictionaries 
    data_dictionary['BOUNDS']  = clipped_data.bounds 
    data_dictionary['X']       = position_x
    data_dictionary['Y']       = position_y
    data_dictionary['Z']       = position_z
    data_dictionary['RHOU']    = data_dictionary['RHOU'].T 
    data_dictionary['U']       = data_dictionary['U'].T
    data_dictionary['VORT']    = data_dictionary['VORT'].T
    data_dictionary['GRADRHO'] = data_dictionary['GRADRHO'].T

# Save data as a pickle file 
    file_out = open(f'{pickle_out}.pickle', "wb") 
    pickle.dump(data_dictionary, file_out)
    file_out.close()

    return data_dictionary
'''


if __name__ == "__main__": 
    start_time = time.time() 
    data_dict  = data_to_binary(bounds, pvtu_in, pickle_out, vtk_out) 
    total_time = round(time.time()  - start_time, 2)
    print(f"----- {total_time} [s] -----")  


