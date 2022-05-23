#!/opt/homebrew/bin/python3
'''
    Date:   05/23/2022
    Author: Martin E. Liza
    File:   helper_class.py
    Def:    

    Author		    Date		Revision
    ----------------------------------------------------
    Martin E. Liza	05/23/2022	Initial version.
'''

from dataclasses import dataclass
import pickle 
import os 
import numpy as np 
from scipy.io import FortranFile  

@dataclass 
class Helper:  
    pass 

# Loading Unformat fortran data  
    def data_loader(self, variable_in, abs_path_in): 
        data_in = os.path.join(abs_path_in, f'{variable_in}.dat')
        f_in = FortranFile(data_in, 'r')
        data = f_in.read_reals(dtype=np.float64)
        f_in.close() 
        return data

# Save data as pickle and loads data as a pickle 
    def pickle_manager(self, pickle_name, pickle_path, data_in=None):
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


