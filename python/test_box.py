import numpy as np 
import IPython 
import base_analysis 
import box 

pickle_path = '../../plate_data/data_12/pickle'
box = box.Box(pickle_path)
TT = box.working_data['data_out']['T'] * - 1 
TT = box.working_data['data_in']['T'] * - 1 
#box.add_new_variable('data_in', 'TT', TT, pickle_path, 'TEST') 
box.derivative(box.working_data['data_out']) 
