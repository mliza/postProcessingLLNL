#!/opt/homebrew/bin/python3
'''
    Date:   09/03/2022
    Author: Martin E. Liza
    File:   box_class.py 
    Def:    Functions used to post process binary 
            box probes output by margot.

    Author           Date         Revision
    ---------------------------------------------------------------------
    Martin E. Liza   07/19/2022   Initial Version.
    Martin E. Liza   07/25/2022   Added mean fields and 
                                  Reynolds decomposition.
    Martin E. Liza   08/17/2022   Added edge properties, wall
                                  properties and Van-Driest transform.
    Martin E. Liza   09/03/2022   Added energy spectrum, correlation
                                  and str_location functions. 
    Martin E. Liza   10/22/2022   Added coarser_field function.
'''
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import IPython
import sys 
import os 
from dataclasses import dataclass, field  
from scipy import integrate 
from scipy.fft import fft 
from scipy.io import FortranFile
from scipy.optimize import curve_fit 
# My own stuff 
scripts_path   = os.environ.get('SCRIPTS')
python_scripts = os.path.join(scripts_path, 'Python') 
sys.path.append(python_scripts) 
import helper_class as helper 

# Probe Class 
@dataclass 
class Box():
# Initialize variables 
    nx: int
    ny: int
    nz: int
    n_total : int = field(init=False)

 # Initialize variables 
    def __post_init__(self):
        self.n_total = self.nx * self.ny * self.nz 
    
# Loads mapping fortran data and saves a pickle file if pickle_path given 
    def mapping_reader(self, mapping_data_in, pickle_path=None):
        mapping_out = np.empty([self.n_total, 4], dtype=int)
        f_in        = FortranFile(mapping_data_in, 'r') 
        for n in range(self.n_total):
            mapping_out[n] = f_in.read_ints() 
        f_in.close() 
        # Save as a pickle file 
        if pickle_path is None:
            return mapping_out
        else:
            helper_scripts = helper.Helper() 
            helper_scripts.pickle_manager(pickle_name_file='mapping', 
                                          pickle_path=pickle_path, 
                                          data_to_save=mapping_out)

# Split the 1D array into a 3D array 
    def split_plot3D(self, array_1D, mapping):
        array_3D = np.empty([self.nx, self.ny, self.nz])
        for n in range(self.n_total):
            # Mapping = [n, i, j, k] 
            i = mapping[n][0]
            j = mapping[n][1]
            k = mapping[n][2]
            array_3D[i,j,k] = array_1D[n]  
        return array_3D

# Coarse mesh
    def coarser_field(self, field_3D, f_width):
        # Spacing vectors, with a filter width increment   
        fx  = range(0, self.nx, f_width)
        fy  = range(0, self.ny, f_width)
        fz  = range(0, self.nz, f_width)
        coarser_field = np.empty([len(fx)-1, len(fy)-1, len(fz)-1])
        # fi, fj, fk increment on the spacing vector 
        # i, j, k increments on the unfiltered field.
        for i, fi in enumerate(fx[:-1]):
            for j, fj in enumerate(fy[:-1]):
                for k, fk in enumerate(fz[:-1]):
                    coarser_field[i,j,k] = np.mean(field_3d[fi:fi+f_width,
                                                            fj:fj+f_width,
                                                            fk:fk+f_width])
        return coarser_field 

# Return a dictionary with gradient fields 
    def gradient_fields(self, array_dict3D):  
        # GRADV_ij = du_j/dx_i = d_i u_j  (not jacobian) 
        u_mag      = np.sqrt(array_dict3D['Ux']**2 + 
                             array_dict3D['Ux']**2 + 
                             array_dict3D['Uz']**2) 

        # Rotation_ij = 1/2 (grad_ji - grad_ij)
        rotation_zy = 1/2 * (array_dict3D['GRADV_23'] - array_dict3D['GRADV_32']) 
        rotation_xz = 1/2 * (array_dict3D['GRADV_31'] - array_dict3D['GRADV_13'])
        rotation_yx = 1/2 * (array_dict3D['GRADV_12'] - array_dict3D['GRADV_21'])

        # Strain_ij = 1/2 (grad_ji + grad_ij)
        strain_yx = 1/2 * (array_dict3D['GRADV_12'] + array_dict3D['GRADV_21']) 
        strain_xz = 1/2 * (array_dict3D['GRADV_31'] + array_dict3D['GRADV_13'])  
        strain_zy = 1/2 * (array_dict3D['GRADV_23'] + array_dict3D['GRADV_32']) 

        # Fibonacci Norm square ||A|| = sqrt(sum |a_ij|^2) 
        # Norms are fibonacy norm square  
        rotation_norm   = ( (2 * rotation_zy)**2 + 
                            (2 * rotation_xz)**2 + 
                            (2 * rotation_yx)**2 )
        shear_norm      = ( (2 * strain_yx)**2 + 
                            (2 * strain_xz)**2 + 
                            (2 * strain_zy)**2 )
        dilatation_norm = ( array_dict3D['GRADV_11']**2 + 
                            array_dict3D['GRADV_22']**2 + 
                            array_dict3D['GRADV_33']**2 ) 

        gradient_dict = { 'rotation_norm'   : rotation_norm, 
                          'shear_norm'      : shear_norm,
                          'dilatation_norm' : dilatation_norm,
                          'ratation_zy'     : rotation_zy,
                          'ratation_xz'     : rotation_xz,
                          'ratation_yx'     : rotation_yx,
                          'shear_yx'        : strain_yx,
                          'shear_xz'        : strain_xz,
                          'shear_zy'        : strain_zy,
                          'UMAG'            : u_mag,
                          'VORTMAG'         : np.sqrt(rotation_norm) }

        return gradient_dict 

# Return fluctuation fields 
    def mean_positions(self, dict_3D):
        # Loading data
        x = dict_3D['X']
        y = dict_3D['Y']
        z = dict_3D['Z']
        # Declaring empty arrays 
        x_mean = np.empty(self.nx) 
        y_mean = np.empty(self.ny) 
        z_mean = np.empty(self.nz) 
        # Calculating the means 
        for k in range(self.nz): 
            z_mean[k] = np.mean(z[:,:,k]) 
        for i in range(self.nx):
            x_mean[i] = np.mean(x[i,:,:])
        for j in range(self.ny):
            y_mean[j] = np.mean(y[:,j,:])
        for k in range(self.nz): 
            z_mean[k] = np.mean(z[:,:,k]) 
        # Returning dictionary 
        dict_out = { 'mean_x' : x_mean,
                     'mean_y' : y_mean,
                     'mean_z' : z_mean }
        return dict_out 

# Edge properties  
    def edge_properties(self, array_field3D, array_height3D, freestream_value):
        edge_dict      = { }
        edge_field     = np.empty([self.nx, self.nz]) 
        edge_thickness = np.empty([self.nx, self.nz]) 
        cut_value      = 0.99 * freestream_value 
        # Find positions at 0.99 freestream 
        for i in range(self.nx):
            for k in range(self.nz): 
                indx                = np.abs(array_field3D[i,:,k] - 
                                          cut_value).argmin() 
                edge_field[i,k]     = array_field3D[i,:,k][indx] 
                edge_thickness[i,k] = array_height3D[i,:,k][indx] 
        # Calculate the mean (compress on z, assume frozen flow)  
        mean_edge_thickness = np.empty(self.nx)
        mean_edge_field     = np.empty(self.nx)
        for i in range(self.nx):
            mean_edge_thickness[i] = np.mean(edge_thickness[i,:])
            mean_edge_field[i]     = np.mean(edge_field[i,:])
        # Dictionary to return 
        edge_dict['edge_field']          = edge_field 
        edge_dict['edge_thickness']      = edge_thickness 
        edge_dict['mean_edge_thickness'] = mean_edge_thickness
        edge_dict['mean_edge_field']     = mean_edge_field  

        return edge_dict 

# Calculates fluctuation fields in a given 3D data set, 
# assume frozen flow hypothesis on z, and returns a 
# 2D field as a function of x and y.  
    def mean_fields(self, array_3D):
        mean_yz = np.empty([self.ny, self.nz]) 
        mean_xz = np.empty([self.nx, self.nz]) 
        mean_xy = np.empty([self.nx, self.ny]) 
        mean_x  = np.empty(self.nx)
        mean_y  = np.empty(self.ny)
        mean_z  = np.empty(self.nz)
        # Array y-z
        for k in range(self.nz):
            for j in range(self.ny):
                mean_yz[j,k] = np.mean(array_3D[:,j,k]) 
        # Array x-z
        for k in range(self.nz):
            for i in range(self.nx):
                mean_xz[i,k] = np.mean(array_3D[i,:,k]) 
        # Array x-y
        for j in range(self.ny):
            for i in range(self.nx):
                mean_xy[i,j] = np.mean(array_3D[i,j,:]) 
        # Mean x 
        for i in range(self.nx): 
            mean_x[i] = np.mean(mean_xy[i,:])
        # Mean y 
        for j in range(self.ny):
            mean_y[j] = np.mean(mean_yz[j,:])
        # Mean z 
        for k in range(self.nz):
            mean_z[k] = np.mean(mean_xz[:,k])

        dict_out = { 'mean_xy' : mean_xy,
                     'mean_yz' : mean_yz,
                     'mean_xz' : mean_xz, 
                     'mean_x'  : mean_x,
                     'mean_y'  : mean_y,
                     'mean_z'  : mean_z }
        return dict_out 

# Reynolds Decomposition 
    def reynolds_decomposition(self, array_3D):
        decomp_3D = np.empty([self.nx, self.ny, self.nz])
        for i in range(self.nx):
            for j in range(self.ny):
                mean_var = np.mean(array_3D[i,j,:])
                for k in range(self.nz):
                    decomp_3D[i,j,k] = array_3D[i,j,k] - mean_var 
        return decomp_3D

# Auto correlation 
    def correlation_function(self, radius, field_1, field_2, correlation_len=50):
        field_len          = len(field_1) - 1
        numerator          = np.zeros(correlation_len) 
        denominator        = np.zeros(correlation_len) 
        radius             -= np.min(radius) 
        correlation_radius = np.linspace(0, np.max(radius), 
                                             correlation_len)
        for i in range(field_len):
            k = i
            for j in range(correlation_len):
                numerator[j]   += (field_1[i] * field_2[k])
                denominator[j] += field_1[i] * field_2[i] 
                k += 1
                if (k > field_len):
                    break 
        correlation_norm = numerator / denominator 
        # Dictionary to return 
        correlation_dict = { 'radius'          : correlation_radius,
                             'norm_correlation': correlation_norm, 
                             'correlation'     : numerator }
        return correlation_dict

# Length scales 
    def microscales(self, correlation_radius, correlation):
        delta_r    = np.mean(np.diff(correlation_radius)) 
        # 2nd order finite forward difference  
        derivative = (2 * autocorrelation[0] - 5 * autocorrelation[1] + 
                    4 * autocorrelation[2] - autocorrelation[3]) / delta_r**2 
        if derivative > 0:
            derivative = (autocorrelation[0] - 2 * autocorrelation[1] +  
                          autocorrelation[2]) / delta_r**2

        taylor_scale   = 1 / np.sqrt(-0.5 * derivative)
        integral_scale = np.abs(integrate.simpson(autocorrelation, dx=delta_r))
        # Dictionary to return
        microscale_dict = { 'integral' : integral_scale,
                            'taylor'   : taylor_scale}
        return microscale_dict 

# Energy Spectrum 
    def energy_spectrum(self, Ux, Uy, Uz, n_elements, n_bins=2):
        # Convert inputs to fourrier space 
        u_hat    = fft(Ux, n_elements)
        v_hat    = fft(Uy, n_elements)
        w_hat    = fft(Uz, n_elements)
        # Power Spectral Density 
        psd      = ( (u_hat * np.conj(u_hat) / n_elements).real + 
                     (v_hat * np.conj(v_hat) / n_elements).real + 
                     (w_hat * np.conj(w_hat) / n_elements).real ) 
        num_bins = int(np.floor(n_elements / n_bins) + 1)
        bin_vec  = range(1, num_bins + 1)
        pwf      = np.zeros(num_bins)
        log_freq = np.log10(bin_vec) 
        # Butterworth filter spectrum with 1/10th 
        for i in range(num_bins - 1):
            count   = 0 
            sum_pwd = 0
            log_freq_0 = log_freq[i] - 0.05 #initial 
            log_freq_t = log_freq[i] + 0.05 #final 
            for j in range(num_bins - 1):
                if (log_freq[j] > log_freq_0 and log_freq[j] < log_freq_t):
                    count   += 1
                    sum_pwd += psd[j]
            if count != 0:
                pwf[i] = sum_pwd / count 
        spectrum = np.append(0, pwf[:-1])
        return spectrum 

# Wall shear-stress 
    def van_driest(self, s12_mean, u_mean, y_mean, rho_mean, mu_mean, t_mean):
        # Three Point interpolation to calculate results at the wall
        rho_w = self.three_point_extrapolation(y_mean['mean_xy'], 
                                               rho_mean['mean_xy'], 
                                               y_loc=0.0) 
        mu_w  = self.three_point_extrapolation(y_mean['mean_xy'], 
                                               mu_mean['mean_xy'], 
                                               y_loc=0.0)
        T_w   = self.three_point_extrapolation(y_mean['mean_xy'], 
                                               t_mean['mean_xy'], 
                                               y_loc=0.0)
        S12_w = self.three_point_extrapolation(y_mean['mean_xy'], 
                                               s12_mean['mean_xy'], 
                                               y_loc=0.0)
        u_1   = self.three_point_extrapolation(y_mean['mean_xy'], 
                                               u_mean['mean_xy'], 
                                               y_loc=-y_mean['mean_xy'][:,0])
        # Second order Central Difference 
        S12_wi = (u_mean['mean_xy'][:,0] - u_1) / (2 * y_mean['mean_xy'][:,0]) 

        # Calculate wall properties 
        nu_w  = mu_w / rho_w  
        tau_w = -mu_w * S12_w 
        tau_wi = -mu_w * S12_wi 
        u_tau = np.sqrt(np.abs(tau_w / rho_w))  
        ui_tau = np.sqrt(np.abs(tau_wi / rho_w))  
        # Calculate van driest transformations and returns at each x-position
        y_plus   = np.empty([self.nx, self.ny]) 
        yi_plus   = np.empty([self.nx, self.ny]) 
        u_plus   = np.empty([self.nx, self.ny]) 
        ui_plus   = np.empty([self.nx, self.ny]) 
        T_plus   = np.empty([self.nx, self.ny]) 
        rho_plus = np.empty([self.nx, self.ny]) 
        for i in range(self.nx):
            y_plus[i,:]   = u_tau[i] * y_mean['mean_xy'][i,:] / nu_w[i] 
            yi_plus[i,:]  = ui_tau[i] * y_mean['mean_xy'][i,:] / nu_w[i] 
            u_plus[i,:]   = u_mean['mean_xy'][i,:] / u_tau[i] 
            ui_plus[i,:]  = u_mean['mean_xy'][i,:] / ui_tau[i] 
            T_plus[i,:]   = t_mean['mean_xy'][i,:] / T_w[i] 
            rho_plus[i,:] = rho_mean['mean_xy'][i,:] / rho_w[i]
        # Dictionary 
        van_driest_dict = { 'y_plus'      : y_plus, 
                            'yi_plus'     : yi_plus, 
                            'u_plus'      : u_plus, 
                            'ui_plus'     : ui_plus, 
                            'T_plus'      : T_plus,
                            'rho_plus'    : rho_plus,
                            'rho_w'       : rho_w, 
                            'mu_w'        : mu_w, 
                            'nu_w'        : nu_w,
                            'tau_w'       : tau_w, 
                            'u_tau'       : u_tau }

        return van_driest_dict 

# Fitting function
    def smoothing_function(self, data_in, box_pts):
        box         = np.ones(box_pts)/box_pts 
        data_smooth = np.convolve(data_in, box, mode='same')
        return data_smooth 

# Plotting locations 
    def str_locations(self, mean_loc_dict, x=None, y=None, 
                      z=None):  
        # Count how many None are as input 
        list  = [x, y, z]
        count_ = len([count for count in list if count is None])
        # Return 2D 
        if count_ == 1:
            if x is None:
                y_val = mean_loc_dict['mean_y'][y] 
                z_val = mean_loc_dict['mean_z'][z] 
                loc_str = f'y = {y_val:.3} $[m]$, z = {z_val:.3} $[m]$'
            if y is None:
                x_val = mean_loc_dict['mean_x'][x] 
                z_val = mean_loc_dict['mean_z'][z] 
                loc_str = f'x = {x_val:.3} $[m]$, z = {z_val:.3} $[m]$'
            if z is None:
                x_val = mean_loc_dict['mean_x'][x] 
                y_val = mean_loc_dict['mean_y'][y] 
                loc_str = f'x = {x_val:.3} $[m]$, y = {y_val:.3} $[m]$'
        # Return 1D 
        else:
            if y is None and z is None:
                x_val = mean_loc_dict['mean_x'][x] 
                loc_str = f'x = {x_val:.3} $[m]$'
            if x is None and z is None:
                y_val = mean_loc_dict['mean_y'][y] 
                loc_str = f'y = {y_val:.3} $[m]$'
            if x is None and y is None:
                z_val = mean_loc_dict['mean_z'][z] 
                loc_str = f'z = {z_val:.3} $[m]$'
        return loc_str 

# Legendre extrapolation 2nd order 
    def three_point_extrapolation(self, y_mean, f_mean, y_loc=0):
        # Legendre Constants 
        const_a  = ( (y_loc - y_mean[:,1]) * (y_loc - y_mean[:,2]) / ( 
                  (y_mean[:,0] - y_mean[:,1]) * (y_mean[:,0] - y_mean[:,2]) )) 

        const_b  = ( (y_loc - y_mean[:,2]) * (y_loc - y_mean[:,0]) / ( 
                  (y_mean[:,1] - y_mean[:,2]) * (y_mean[:,1] - y_mean[:,0]) )) 

        const_c  = ( (y_loc - y_mean[:,0]) * (y_loc - y_mean[:,1]) / ( 
                  (y_mean[:,2] - y_mean[:,0]) * (y_mean[:,2] - y_mean[:,1]) )) 

        # Extrapolation
        extrapolation = ( const_a * f_mean[:,0] + 
                          const_b * f_mean[:,1] + 
                          const_c * f_mean[:,2] )
        return extrapolation 

# Time average 
    def time_average(self, field_matrix, sub_sample_spacing=1):
        # field_matrix = [time, elements]
        [time_len, element_len] = np.shape(field_matrix)
        average_field           = np.empty(element_len)
        sub_sample_time         = range(0, time_len, sub_sample_spacing) 
        for i in range(element_len):
            average_field[i] = np.mean(field_matrix[sub_sample_time, i])

        return average_field 

        
        

# Van Driest plot 
    def plot_van_driest(self, y_plus, u_plus, title_in, testing_path=None, 
                        saving_path=None, fig_name=None): 
        plt.plot(y_plus, u_plus,
                 color='k', linestyle='-', linewidth=2, marker='o', markersize=5,
                 markerfacecolor='lightgrey', markeredgecolor='k', 
                 label='MARGOT, M10') 
        u_form = 2.44 * np.log(y_plus) + 5.2 
        u_mine = 2.44 * np.log(y_plus) + 8.72
        '''
        plt.plot(y_plus, u_form, 's', markersize=2.5, 
                 label='$2.44\,ln(y^+) + 5.2$') 
        plt.plot(y_plus[:6], y_plus[:6],'o',  markersize=2.5, label='$y^+$')  
        plt.plot(y_plus[5:-40], u_mine[5:-40], '-.', markersize=2.5, color='darkred', 
                 #label='$6.63\,ln(y^+) - 6.28$') 
                 label='$2.44\,ln(y^+) + 8.72$') 
        plt.legend() 
        '''
        plt.title(title_in)

        if testing_path != None: 
            testing_file = os.path.join(testing_path,
                        'vanDriestTransformation_1fig11.csv')
            df           = pd.read_csv(testing_file) 
            y_plus       = np.array(df['y_plus'])
            u_plus       = np.array(df['u_plus'])
            plt.plot(y_plus, u_plus, color='darkorange' ,
                     linewidth=1.5, linestyle='--', label='Pino, M5')
            plt.legend() 
        plt.xscale('log')
        plt.grid('-.')
        plt.xlabel('$y^+$')
        plt.ylabel('$u^+$')

        # Saving 
        if saving_path == None:
            plt.show() 

        if saving_path != None:
            plt.tight_layout()
            if fig_name == None:
                plt.savefig(f'{saving_path}/van_driest.png', dpi=300)
            if fig_name != None:
                plt.savefig(f'{saving_path}/{fig_name}.png', dpi=300)
            plt.close() 

## DELETE WHEN IT IS DONE ## 
# Plot boundary Layers 
    def plot_boundary_layers(self, velocity_boundary_dict, 
            temperature_boundary_dict, mean_velocity, mean_temperature, 
            grid_mean_dict, velocity_freestream, temperature_freestream,
            saving_path=None):
        # Loading variables 
        n_box           = 50  
        temp_mean_thick = temperature_boundary_dict['mean_edge_thickness'] * 10**3
        vel_mean_thick  = velocity_boundary_dict['mean_edge_thickness'] * 10**3
        temp_mean_field = mean_temperature 
        vel_mean_field  = mean_velocity  
        x_mean          = grid_mean_dict['mean_x'] * 10**2
        y_mean          = grid_mean_dict['mean_y'] * 10**3
        v_color         = 'mediumturquoise' 
        t_color         = 'darkorange'
        temp_smooth     = self.smoothing_function(temp_mean_thick, n_box) 
        vel_smooth      = self.smoothing_function(vel_mean_thick, n_box) 
        n_box           /= 2
        n_box           = int(n_box) 
        # Plotting figures  
        fig, (ax1, ax2) = plt.subplots(1,2, figsize=(8,5))
        # Plot thickness 
        ax1.plot(x_mean, vel_mean_thick, 'o', markersize=3,
                markerfacecolor='lightgrey', markeredgecolor='k') 
        ax1.plot(x_mean[n_box:-n_box], vel_smooth[n_box:-n_box], color=v_color, 
                 linestyle='-', linewidth=1.5, label='Ux')
        ax1.plot(x_mean, temp_mean_thick, 'o', markersize=3, 
                markerfacecolor='lightgrey', markeredgecolor='k')
        ax1.plot(x_mean[n_box:-n_box], temp_smooth[n_box:-n_box], color=t_color, 
                 linestyle='-', linewidth=1.5, label='T') 
        ax1.legend()
        ax1.set_ylabel('y-axis [mm]')
        ax1.set_xlabel('x-axis [cm]')
        ax1.grid('-.') 
        ax1.legend() 
        # Plot value 
        l1 = ax2.plot(vel_mean_field, y_mean, color=v_color, 
                linestyle='-', linewidth=3, label=f'Ux={velocity_freestream:.2f}[m/s]')
        ax2.set_xlabel('Ux [m/s]', color=v_color)
        ax21 = ax2.twiny() 
        l2 = ax21.plot(temp_mean_field, y_mean, color=t_color, 
                linestyle='-', linewidth=3, label=f'T={temperature_freestream:.2f}[K]')
        ax21.set_xlabel('T [K]', color=t_color)
        ax2.set_ylabel('y-axis [mm]')
        ax2.grid('-.') 
        ax2.legend(handles=l1+l2, loc='upper center') 
        
        if saving_path == None:
            plt.show() 
        if saving_path != None:
            fig.tight_layout()
            fig.savefig(f'{saving_path}/boundary_layers.png', dpi=300)
            plt.close() 

# Making contour plots 
    def plot_contour(self, data_dict3D, grid_dict3D, grid_x, grid_y, 
                field, slice_cut, slice_direction,
                levels=6, cmap='inferno', saving_path=None): 
        # Create slides 
        slice_direction = slice_direction.upper() 
        if slice_direction == 'X':
            x_plane     = grid_dict3D[grid_x][slice_cut,:,:]
            y_plane     = grid_dict3D[grid_y][slice_cut,:,:]
            z_plane     = data_dict3D[field][slice_cut,:,:]
            slice_value = grid_dict3D[slice_direction][:,-1,-1][slice_cut] 
        if slice_direction == 'Y':
            x_plane     = grid_dict3D[grid_x][:,slice_cut,:]
            y_plane     = grid_dict3D[grid_y][:,slice_cut,:]
            z_plane     = data_dict3D[field][:,slice_cut,:]
            slice_value = grid_dict3D[slice_direction][-1,:,-1][slice_cut] 
        if slice_direction == 'Z':
            x_plane     = grid_dict3D[grid_x][:,:,slice_cut]
            y_plane     = grid_dict3D[grid_y][:,:,slice_cut]
            z_plane     = data_dict3D[field][:,:,slice_cut]
            slice_value = grid_dict3D[slice_direction][-1,-1,:][slice_cut] 
        # Plotting 
        plt.contourf(x_plane, y_plane, z_plane, 
                    levels=levels, cmap=cmap)

        plt.xlabel(f'{grid_x} [m]')
        plt.ylabel(f'{grid_y} [m]')
        plt.title(f'{field}, at {slice_direction}={slice_value:.3E} [m]') 
        plt.colorbar() 
        # Saving if needed 
        if saving_path == None:
            plt.show() 
        if saving_path != None:
            plt.savefig(f'{saving_path}/contour{grid_x}{grid_y}_{field}.png', 
                            bbox_inches='tight', dpi=300)
            plt.close() 

# Plot line for 2 variables  
    def plot_lineXY(self, array_dict_3D, var_x, var_y, 
                    x_dim=None, y_dim=None, z_dim=None, saving_path=None):
        if x_dim is None:
            x_axis    = array_dict_3D[var_x][:, y_dim, z_dim] 
            y_axis    = array_dict_3D[var_y][:, y_dim, z_dim] 
            label_str = f'y={y_dim}, z={z_dim}'
        if y_dim is None:
            x_axis    = array_dict_3D[var_x][x_dim, :, z_dim] 
            y_axis    = array_dict_3D[var_y][x_dim, :, z_dim] 
            label_str = f'x={x_dim}, z={z_dim}'
        if z_dim is None:
            x_axis    = array_dict_3D[var_x][x_dim, y_dim, :] 
            y_axis    = array_dict_3D[var_y][x_dim, y_dim, :]
            label_str = f'x={x_dim}, y={y_dim}'
        # Legend, title 
        plt.plot(x_axis, y_axis, color='k',  linestyle='-', linewidth=3,
                 label=label_str) 
        plt.grid('-.') 
        plt.legend() 
        plt.xlabel(f'{var_x}')
        plt.ylabel(f'{var_y}')
        plt.title(f'{var_y} vs. {var_x}') 
        # Saving if needed 
        if saving_path == None:
            plt.show() 
        if saving_path != None:
            plt.savefig(f'{saving_path}/{var_x}_{var_y}.png', dpi=300)
            plt.close() 

# Plot mean_fields 
    def plot_mean_fields(self, x_axis, y_axis, x_str, y_str, saving_path=None):
        # Legend, title 
        plt.plot(x_axis, y_axis, color='k',  linestyle='-', linewidth=3)
        plt.grid('-.') 
        plt.xlabel(f'{x_str}')
        plt.ylabel(f'{y_str}')
        # Saving if needed 
        if saving_path == None:
            plt.show() 
        if saving_path != None:
            plt.savefig(f'{saving_path}/{x_str}_{y_str}.png', dpi=300)
            plt.close() 

# Plot boundary surface 
    def plot_boundary_surface(self, boundary_plane_dict, 
                              grid_mean_dict):  
        # Loading data 
        X_mean        = grid_mean_dict['mean_x'] 
        Z_mean        = grid_mean_dict['mean_z'] 
        height        = boundary_plane_dict['thickness']
        boundary_mean = boundary_plane_dict['mean_thickness']  
        X,Y           = np.meshgrid(Z_mean, X_mean)

        #fig = plt.figure(figsize=(8,8))
        fig = plt.figure()
        ax  = fig.add_subplot(111, projection='3d') 
        ax.plot_surface(X, Y, height) 
