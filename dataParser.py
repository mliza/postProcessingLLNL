#!/opt/homebrew/bin/python3 
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
probe_flag     = False 
line_flag      = True 
script_path    = os.getcwd() 
directory_path = os.path.dirname(script_path) 
data_path      = os.path.join(directory_path, 'plate_data', 'data_7')
saving_path    = os.path.join(data_path, 'pickle')  

def probe_parser(data_path):
# Paths 
    probe_path    = os.path.join(data_path, 'PROBES')
# List all input files 
    total_data_in = os.listdir(probe_path) 
    file_in       = [ ]
# Create a list with all the files starting with p and not ending with README 
    for i in total_data_in: 
        if re.search('^p', i):
            if not re.search('README$', i):
                file_in.append(i) 
# Loads and stores data as multi dictionary 
    data = { }
    for i in file_in:
        file_in_abs_path = os.path.join(probe_path, i)
        pd_in = pd.read_csv(file_in_abs_path, sep='\t', header=None) 
        array_in = np.array(pd_in.iloc[:][0]) 
        [iteration, variable] = i.split('.') 
        # Creates a temporary vectors to store the data 
        temp             = [ ] 
        iteration_number = [ ]
        time_in          = [ ]
        # Loops through the input array and stores into the temp 
        # and iteration_number vector 
        for j in array_in:
            iteration_number.append(int(j.split()[0])) 
            time_in.append(float(j.split()[1])) 
            temp.append(float(j.split()[3])) 
        # Creates an empty dictionary for each iteration and add the  
        # iteration number to it 
        if not iteration in data: 
            data[iteration] = { } 
            data[iteration]['ITER'] = np.array(iteration_number)
            data[iteration]['TIME'] = np.array(time_in) 
        # Stores the data 
        data[iteration][variable] = np.array(temp)  
    return data 

def line_parser(data_path):
    probe_path    = os.path.join(data_path, 'PROBES')
# List all input files 
    total_data_in = os.listdir(probe_path) 
    file_in       = [ ]
# Create a list with all the files starting with p and not ending with README 
    for i in total_data_in: 
        if re.search('^l', i):
            if not re.search('README$', i):
                file_in.append(i) 
# Loads and stores data as multi dictionary 
    data = { }
#array_in = [ ]
    for i in file_in: 
        file_in_abs_path = os.path.join(probe_path, i)
        pd_in    = pd.read_csv(file_in_abs_path, sep='\t', header=None) 
        pd_len   = len(pd_in)
        [iteration, variable] = i.split('.') 
        # Create a matrix from the variable, and vectors for time and iteration 
        for l in range(pd_len):
            temp_var = [ ]
            array_in = pd_in.iloc[l][0].split() 
            # Load data and converts it to array, float 64 
            iter_elem = int(array_in[0])
            time_elem = float(array_in[1])
            for m in range(3, len(array_in)): 
                temp_var.append(float(array_in[m])) 
            # Only initialize the first element 
            if l == 0:
                iter_vect  = [ ]
                time_vect  = [ ]
                var_matrix = np.empty(shape=(pd_len, len(temp_var))) 
                var_matrix = [ ] 
            # Populating vectors and matrix 
            iter_vect.append(iter_elem) 
            time_vect.append(time_elem) 
            var_matrix.append(temp_var) 
        # Store data  
        if not iteration in data: 
            data[iteration] = { } 
            data[iteration]['ITER'] = np.array(iter_vect) 
            data[iteration]['TIME'] = np.array(time_vect)  
        data[iteration][variable]   = np.array(var_matrix) 
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
        data_line = line_parser(data_path) 
        pickleOut = open(f'{saving_path}/line_data.pickle', 'wb')
        pickle.dump(data_line, pickleOut)
        pickleOut.close() 
