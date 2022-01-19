#!/opt/homebrew/bin/python3.9 
'''
    Date:   07/26/2021
    Author: Martin E. Liza
    File:   dataParser.py
    Def:    parses probe or line data 

    Author		    Date		Revision
    ----------------------------------------------------
    Martin E. Liza	07/26/2021	Initial version.
'''
import numpy as np
import pandas as pd 
import os
import re
import sys 
import pickle
import IPython

# Paths setup 
probe_flag     = True  
line_flag      = False 
script_path    = os.getcwd() 
directory_path = os.path.dirname(script_path) 
data_path      = os.path.join(directory_path, 'plate_data', 'data_1')
saving_path    = os.path.join(data_path, 'pickle')  

def probe_parser(data_path):
# Paths 
    probe_path     = os.path.join(data_path, 'PROBES')
    
# List all input files 
    total_data_in     = os.listdir(probe_path) 
    total_num_data_in = len(total_data_in) 
    file_in = [ ]

# Create a list with all the files starting with p and not ending with README 
    for i in range(total_num_data_in):
        file_name = total_data_in[i] 
        if re.search('^p', file_name):
            if not re.search('README$', file_name):
                file_in.append(file_name) 

# Loads and stores data as multi dictionary 
    data = { }
    for i in range(len(file_in)):
        file_in_abs_path = os.path.join(probe_path, file_in[i])
        pd_in = pd.read_csv(file_in_abs_path, sep='\t', header=None) 
        array_in = np.array(pd_in.iloc[:][0]) 
        [iteration, variable] = file_in[i].split('.') 
            
        # Creates a temporary vectors to store the data 
        temp = np.empty(np.shape(array_in)[0], dtype=np.float64)
        iteration_number = np.empty(np.shape(array_in)[0], dtype=np.float64)
        time = np.empty(np.shape(array_in)[0], dtype=np.float64)
        # Loops through the input array and stores into the temp 
        # and iteration_number vector 
        for j in range(np.shape(array_in)[0]):
            temp[j] = array_in[j].split()[3] 
            iteration_number[j] = array_in[j].split()[0] 
            time[j] = array_in[j].split()[1] 

        # Creates an empty dictionary for each iteration and add the  
        # iteration number to it 
        if not iteration in data: 
            data[iteration] = { } 
            data[iteration]['ITER'] = iteration_number 
            data[iteration]['TIME'] = time  

        # Stores the data 
        data[iteration][variable] = temp  
    return data 

def line_parser(data_path):
    probe_path     = os.path.join(data_path, 'PROBES')

# List all input files 
    total_data_in     = os.listdir(probe_path) 
    total_num_data_in = len(total_data_in) 
    file_in = [ ]

# Create a list with all the files starting with p and not ending with README 
    for i in range(total_num_data_in):
        file_name = total_data_in[i] 
        if re.search('^l', file_name):
            if not re.search('README$', file_name):
                file_in.append(file_name) 

# Loads and stores data as multi dictionary 
    data = { }
#array_in = [ ]
    for i in range(len(file_in)):
        file_in_abs_path = os.path.join(probe_path, file_in[i])
        pd_in = pd.read_csv(file_in_abs_path, sep='\t', header=None) 
        [iteration, variable] = file_in[i].split('.') 

        # Create a matrix from the variable, and vectors for time and iteration 
        for l in range(len(pd_in)):
            array_in = pd_in.iloc[l][0].split() 

            # Load data and converts it to array, float 64 
            iter_elem = array_in[0]
            time_elem = array_in[1]
            var_vect  = array_in[3:-1]

            # Only initialize the first element 
            if l == 0:
                iter_vect  = [ ]
                time_vect  = [ ]
                var_matrix = [ ] 

            # Populating vectors and matrix 
            iter_vect.append(iter_elem) 
            time_vect.append(time_elem) 
            var_matrix.append(var_vect) 

            # Convert strings to arrays/flot64 before finishing the iteration set 
            if l == (len(pd_in)-1):
                iter_vect  = np.array(iter_vect).astype(np.float64) 
                time_vect  = np.array(time_vect).astype(np.float64) 
                var_matrix = np.array(var_matrix).astype(np.float64) 

        # Creates an empty dictionary for each iteration and stores them in a  
        # multi-dimensional dictionary 
        if not iteration in data: 
            data[iteration] = { } 
            data[iteration]['ITER'] = iter_vect  
            data[iteration]['TIME'] = time_vect  
            
        data[iteration][variable] = var_matrix 
    return data 

if __name__=="__main__": 
# Probe Data  
    if probe_flag:
        data_probe = probe_parser(data_path) 
        pickleOut  = open(f'{saving_path}/probe_data.pickle', 'wb')
        pickle.dump(data_probe, pickleOut)
        pickleOut.close() 

# Line Data 
    if line_flag:
        data_line= line_parser(data_path) 
        pickleOut = open(f'{saving_path}/line_data.pickle', 'wb')
        pickle.dump(data_line, pickleOut)
        pickleOut.close() 
