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

    # Calculates the derivative of a multidimensional array 
    def spatial_derivative(self, variable_3D): 
        [nx, ny, nz]  = np.shape(variable_3D['X'])
        d_total       = np.empty([nx-1, ny-1, nz-1]) 
        # Calculate dx  
        for j in range(ny-1):
            for k in range(nz-1):
                d_total[:,j,k] = np.diff(variable_3D['X'][:,j,k]) 
        # Calculate dy 
        for i in range(nx-1):
            for k in range(nz-1):
                d_total[i,:,k] = np.diff(variable_3D['Y'][i,:,k]) 
        # Calculate dz
        for i in range(nx-1):
            for j in range(ny-1):
                d_total[i,j,:] = np.diff(variable_3D['Z'][i,j,:]) 
        return d_total 

    # Calculates de derivative of a 3D field, for a Cartesian mesh  
    def derivative_3D(self, variable_3D):
        [nx, ny, nz] = np.shape(variable_3D) 
        d_total      = np.empty([nx-1, ny-1, nz-1]) 
        # Calculate dx 
        for j in range(ny-1):
            for k in range(nz-1):
                d_total[:,j,k] = np.diff(variable_3D[:,j,k]) 
        # Calculate dy 
        for i in range(nx-1):
            for k in range(nz-1):
                d_total[i,:,k] = np.diff(variable_3D[i,:,k]) 
        # Calculate dz 
        for i in range(nx-1):
            for j in range(ny-1):
                d_total[i,j,:] = np.diff(variable_3D[i,j,:]) 
        return d_total 

    # Calculates the Jacobian 
    # https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
    def jacobian_3D(self, ds, dUx, dUy, dUz): 
        [nx, ny, nz] = np.shape(ds)  
        # Initialize 
        dxdUx  = np.empty([nx, ny, nz]) 
        dxdUy  = np.empty([nx, ny, nz]) 
        dxdUz  = np.empty([nx, ny, nz]) 
        dydUx  = np.empty([nx, ny, nz]) 
        dydUy  = np.empty([nx, ny, nz]) 
        dydUz  = np.empty([nx, ny, nz]) 
        dzdUx  = np.empty([nx, ny, nz]) 
        dzdUy  = np.empty([nx, ny, nz]) 
        dzdUz  = np.empty([nx, ny, nz]) 
        # Calculate dUx 
        for j in range(ny-1):
            for k in range(nz-1):
                dxdUx[:,j,k] = dUx[:,j,k] / ds[:,j,k] 
                dxdUy[:,j,k] = dUy[:,j,k] / ds[:,j,k] 
                dxdUz[:,j,k] = dUz[:,j,k] / ds[:,j,k] 
        # Calculate dUy 
        for i in range(nx-1):
            for k in range(nz-1):
                dydUx[i,:,k] = dUx[i,:,k] / ds[i,:,k] 
                dydUy[i,:,k] = dUy[i,:,k] / ds[i,:,k]  
                dydUz[i,:,k] = dUz[i,:,k] / ds[i,:,k] 
        # Calculate dUz 
        for i in range(nx-1):
            for j in range(ny-1):
                dzdUx[i,j,:] = dUx[i,j,:] / ds[i,j,:]  
                dzdUy[i,j,:] = dUy[i,j,:] / ds[i,j,:]  
                dzdUz[i,j,:] = dUz[i,j,:] / ds[i,j,:]  
        IPython.embed(colors = 'Linux') 
'''  
        J = [ dxUx, dydUx, dzdUx,
              dxUy, dydUy, dzdUy,
              dxUz, dydUz, dzdUz ]
'''





        
        


    
   

        



    
