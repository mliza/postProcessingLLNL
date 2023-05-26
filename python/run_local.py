#!/opt/homebrew/bin/python3.9
'''
    Date:   09/27/2022
    Author: Martin E. Liza
    File:   main_box.py 
    Def:    Main file to process box data. 

    Author		    Date		Revision
    ----------------------------------------------------
    Martin E. Liza	07/22/2022	Initial version.
'''
import IPython 
import sys 
import os 
import numpy as np 
import matplotlib.pyplot as plt 
scripts_path   = os.environ.get('SCRIPTS')
python_scripts = os.path.join(scripts_path, 'Python')
sys.path.append(python_scripts)

# Helper class 
import helper_class as helper 
import aerodynamics_class as aero
import box_class as box  
import box_plots 

# User Inputs 
data_path         = '/Users/martin/Desktop/workingFiles/grid'
pino_path         = '../../plate_data/papers_data/pino_martin'
saving_path       = '/Users/martin/Desktop'
pickle_path       = os.path.join(data_path, 'pickle')
nx                = 1439
ny                = 89
nz                = 638 
U_init            = 3000        #[m/s] 
T_init            = 216.66      #[K] 
RHO_init          = 0.18874     #[kg/m3] 
P_init            = 11737       #[Pa] 
X_init            = 11.5970E-2  #[m]


# Loading my classes 
helper = helper.Helper()
aero   = aero.Aero()
box    = box.Box(nx=nx, ny=ny, nz=nz)

# Loading full fields and grid  
fields = helper.pickle_manager(pickle_name_file='temporal_average', 
                                  pickle_path=pickle_path)
grid = helper.pickle_manager(pickle_name_file='grid_3D', 
                                   pickle_path=pickle_path)

# Grid positions 
mean_grid  = box.mean_positions(grid) 

# Calculate Van Driest 
van_driest = box.van_driest(box.mean_fields(fields['GRADV_12'])['mean_xy'],
                            box.mean_fields(fields['Ux'])['mean_xy'],
                            box.mean_fields(grid['Y'])['mean_xy'],
                            box.mean_fields(fields['RHO'])['mean_xy'],
                            box.mean_fields(fields['MU'])['mean_xy'],
                            box.mean_fields(fields['T'])['mean_xy'])

# Norm [1/T] 
norm_tau = van_driest['nu_w']/van_driest['u_tau']**2 
x1 = int(nx/2)
x2 = int(1000)
x3 = int(1300)


# Strings title  
x1_str = box.str_locations(mean_grid, x=x1, y=None, z=None)
x2_str = box.str_locations(mean_grid, x=x2, y=None, z=None)
x3_str = box.str_locations(mean_grid, x=x3, y=None, z=None)

# Gradient fields  
dilatation      = box.mean_fields(fields['DIL'])['mean_xy'] 
dilatation_norm = box.mean_fields(fields['dilatation_norm'])['mean_xy']
rotation_norm   = box.mean_fields(fields['rotation_rate_norm'])['mean_xy'] 
shear_rate_norm = box.mean_fields(fields['shear_rate_norm'])['mean_xy'] 
pure_shear_norm = box.mean_fields(fields['pure_shear_norm'])['mean_xy'] 
Q_norm          = box.mean_fields(fields['Q'])['mean_xy'] 
vortmag         = box.mean_fields(fields['VORTMAG'])['mean_xy'] 

# Dilatation 
plt.plot(van_driest['y_plus'][x1,:], dilatation[x1,:], 
         linewidth=2.5, label=x1_str) 
plt.plot(van_driest['y_plus'][x2,:], dilatation[x2,:], 
         linewidth=2.5, label=x2_str) 
plt.plot(van_driest['y_plus'][x3,:], dilatation[x3,:], 
         linewidth=2.5, label=x3_str) 
plt.ylabel('$\\Theta\;\;[Hz]$')
plt.xlabel('$y^+$')
plt.xscale('log') 
plt.legend() 
plt.grid('.-')
plt.tight_layout()
plt.savefig(os.path.join(saving_path,'dilatation.png'), dpi=300)
plt.close() 

# Dilatation Norm 
plt.plot(van_driest['y_plus'][x1,:], dilatation_norm[x1,:] * norm_tau[x1]**2, 
         linewidth=2.5, label=x1_str) 
plt.plot(van_driest['y_plus'][x2,:], dilatation_norm[x2,:] * norm_tau[x2]**2, 
         linewidth=2.5, label=x2_str) 
plt.plot(van_driest['y_plus'][x3,:], dilatation_norm[x3,:] * norm_tau[x3]**2, 
         linewidth=2.5, label=x3_str) 
plt.ylabel('$(\\Theta\\Theta)^+$')
plt.xlabel('$y^+$')
plt.xscale('log') 
plt.legend() 
plt.grid('.-')
plt.tight_layout()
plt.savefig(os.path.join(saving_path,'normalized_dilatation.png'), dpi=300)
plt.close() 

# Rotation rate
plt.plot(van_driest['y_plus'][x1,:], vortmag[x1,:], 
         linewidth=2.5, label=x1_str) 
plt.plot(van_driest['y_plus'][x2,:], vortmag[x2,:],
         linewidth=2.5, label=x2_str) 
plt.plot(van_driest['y_plus'][x3,:], vortmag[x3,:], 
         linewidth=2.5, label=x3_str) 
plt.ylabel('$||\\Omega||\;\;[Hz]$')
plt.xlabel('$y^+$')
plt.xscale('log') 
plt.legend() 
plt.grid('.-')
plt.tight_layout()
plt.savefig(os.path.join(saving_path,'rotation_rate.png'), dpi=300)
plt.close() 

# Normalized rotation rate
plt.plot(van_driest['y_plus'][x1,:], rotation_norm[x1,:] * norm_tau[x1]**2,
         linewidth=2.5, label=x1_str) 
plt.plot(van_driest['y_plus'][x2,:], rotation_norm[x2,:] * norm_tau[x2]**2,
         linewidth=2.5, label=x2_str) 
plt.plot(van_driest['y_plus'][x3,:], rotation_norm[x3,:] * norm_tau[x3]**2, 
         linewidth=2.5, label=x3_str) 
plt.ylabel('$(\\Omega\\Omega)^+$')
plt.xlabel('$y^+$')
plt.xscale('log') 
plt.legend() 
plt.grid('.-')
plt.tight_layout()
plt.savefig(os.path.join(saving_path,'normalized_rotation_rate.png'), dpi=300)
plt.close() 

# Shearing rate  
plt.plot(van_driest['y_plus'][x1,:], np.sqrt(shear_rate_norm[x1,:]), 
         linewidth=2.5, label=x1_str) 
plt.plot(van_driest['y_plus'][x2,:], np.sqrt(shear_rate_norm[x2,:]), 
         linewidth=2.5, label=x2_str) 
plt.plot(van_driest['y_plus'][x3,:], np.sqrt(shear_rate_norm[x3,:]),
         linewidth=2.5, label=x3_str) 
plt.ylabel('$||S||\;\;[Hz]$')
plt.xlabel('$y^+$')
plt.xscale('log') 
plt.legend() 
plt.grid('.-')
plt.tight_layout()
plt.savefig(os.path.join(saving_path,'shearing_rate.png'), dpi=300)
plt.close() 

# Pure shearing norm 
plt.plot(van_driest['y_plus'][x1,:], pure_shear_norm[x1,:] * norm_tau[x1]**2,
         linewidth=2.5, label=x1_str) 
plt.plot(van_driest['y_plus'][x2,:], pure_shear_norm[x2,:] * norm_tau[x2]**2,
         linewidth=2.5, label=x2_str) 
plt.plot(van_driest['y_plus'][x3,:], pure_shear_norm[x3,:] * norm_tau[x3]**2,
         linewidth=2.5, label=x3_str) 
plt.ylabel('$(\\mathcal{T}\\mathcal{T})^+$')
plt.xlabel('$y^+$')
plt.xscale('log') 
plt.legend() 
plt.grid('.-')
plt.tight_layout()
plt.savefig(os.path.join(saving_path,'pure_shearing_norm.png'), dpi=300)
plt.close() 

# Q Norm
plt.plot(van_driest['y_plus'][x1,:], Q_norm[x1,:], linewidth=2.5, label=x1_str) 
plt.plot(van_driest['y_plus'][x2,:], Q_norm[x2,:], linewidth=2.5, label=x2_str) 
plt.plot(van_driest['y_plus'][x3,:], Q_norm[x3,:], linewidth=2.5, label=x3_str) 
plt.ylabel('$Q\;\;[Hz^2]$')
plt.xlabel('$y^+$')
plt.xscale('log') 
plt.legend() 
plt.grid('.-')
plt.tight_layout()
plt.savefig(os.path.join(saving_path,'Q_crit.png'), dpi=300)
plt.close() 
