#!/opt/homebrew/bin/python3
'''
    Date:   02/24/2022
    Author: Martin E. Liza
    File:   subDataProcess.py
    Def:

    Author		    Date		Revision
    ----------------------------------------------------
    Martin E. Liza	02/24/2022	Initial version.
'''
import numpy as np 
import dataParser as dp 
import IPython 
import pickle 
import os 

# Paths 
script_path = os.getcwd() 
directory_path = os.path.dirname(script_path)
data_path = os.path.join(directory_path, 'plate_data', 'data_7') 
saving_path  = os.path.join(data_path, 'pickle') 

# Loading parsing function 
line_data = dp.line_parser(data_path) 
line_keys = list(line_data.keys()) 
var_keys  = list(line_data[line_keys[0]].keys()) 
# Remove ITER and TIME for the for loops 
line_keys.sort() 
var_keys.remove('ITER')
var_keys.remove('TIME')
[time, space] = np.shape(line_data[line_keys[0]][var_keys[0]])
new_data = { }
for i in range(7): #iterates through x-position 
    new_data[f'l{i}'] = { }
    for k in var_keys: #iterates through variables 
        new_data[f'l{i}'][k] = [ ] 
        for t in range(time): #iterates through time steps 
            temp_time = [ ]
            for s in range(space): #iterates through space steps 
                temp_space = [ ]
                for j in range(6): #iterates through z-position  
                    temp_space.append(line_data[f'l{i}{j}'][k][t][s])  
                temp_time.append(np.mean(temp_space)) 
            new_data[f'l{i}'][k].append(np.array(temp_time))
        # Add iteration and time 
        new_data[f'l{i}'][k]  = np.array(new_data[f'l{i}'][k])
    new_data[f'l{i}']['ITER'] = line_data[f'l{i}{j}']['ITER'] 
    new_data[f'l{i}']['TIME'] = line_data[f'l{i}{j}']['TIME'] 
pickleOut  = open(f'{saving_path}/line_data.pickle', 'wb') 
pickle.dump(new_data, pickleOut)
pickleOut.close() 
