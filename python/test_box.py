import numpy as np 
import IPython 
import base_analysis 
import box 

pickle_path = '../../plate_data/data_13/pickle'
box = box.Box(pickle_path)
#box.add_new_variable('data_in', 'TT', TT, pickle_path, 'TEST') 
ds  = box.spatial_derivative(box.working_data['data_out']) 
dT  = box.derivative_3D(box.working_data['data_out']['T'])  
dUx = box.derivative_3D(box.working_data['data_out']['Ux'])  
dUy = box.derivative_3D(box.working_data['data_out']['Uy'])  
dUz = box.derivative_3D(box.working_data['data_out']['Uz'])  
J   = box.jacobian_3D(ds, dUx, dUy, dUz) 
