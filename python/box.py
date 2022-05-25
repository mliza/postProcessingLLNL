#!/usr/local/bin/python3
'''
    Date:   06/24/2021
    Author: Martin E. Liza
    File:   probe.py Def:               

    Author           Date         Revision
    ----------------------------------------------------------------------- 
    Martin E. Liza   06/24/2021   Initial Version.
    Martin E. Liza   07/14/2021   Fixed auto_correlation. 
    Martin E. Liza   07/15/2021   Added sub sub_sampling function 
                                  and cleaned up code. 
    Martin E. Liza   07/16/2021   Deleted sub_sampling, cleaned up x_axis
                                  to only use <U-x> to calculate the x_axis.  
    Martin E. Liza   07/17/2021   Deleted length_scales, moved this capability
                                  to the base class.
    Martin E. Liza   07/27/2021   Added plot_results function.  
    Martin E. Liza   07/28/2021   Added low_pass_filter function.  
    Martin E. Liza   08/03/2021   Added boxcar_filter and removed 
                                  low_pass_filter.  
    Martin E. Liza   08/24/2021   Added legendre_interpolation function.
    Martin E. Liza   09/16/2021   Cleaned boxcar_filter and plots. 
    Martin E. Liza   10/08/2021   Cleaned death code, renamed plots and 
                                  added plot_scatter function.
'''
import numpy as np 
import IPython
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec 
import seaborn as sb
import os 
# From my baseAnalysis class 
from base_analysis import Base_Analysis

# Probe Class 
class Box(Base_Analysis): 
# Loads 'probe_data'
    flag_type   = None 
    flag_points = 'probe_points'

    # Calculate Derivatives
    def derivative(self, variable_3D): 
        [nx, ny, nz]  = np.shape(variable_3D['X'])
        dx = np.empty([nx-1, ny-1, nz-1]) 
        dy = np.empty([nx-1, ny-1, nz-1]) 
        dz = np.empty([nx-1, ny-1, nz-1]) 
        # X-constant first  
        for j in range(ny-1):
            for k in range(nz-1):
                dx[:,j,k] = np.diff(variable_3D['X'][:,j,k]) 
        # Y-constant first  
        for i in range(nx-1):
            for k in range(nz-1):
                dy[i,:,k] = np.diff(variable_3D['Y'][i,:,k]) 
        # Z-constant first  
        for i in range(nx-1):
            for j in range(ny-1):
                dz[i,j,:] = np.diff(variable_3D['Z'][i,j,:]) 
        
        IPython.embed(colors='Linux') 



    
