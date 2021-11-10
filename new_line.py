#!/usr/local/bin/python3.9
'''
    Date:   06/24/2021
    Author: Martin E. Liza
    File:   line.py
    Def:               

    Author          Date        Revision
    ---------------------------------------------------- 
    Martin E. Liza  06/24/2021  Initial version.
    Martin E. Liza  07/15/2021  Fixed cross-correlation.  
    Martin E. Liza  07/20/2021  Fixed filter_decay and cleaned
                                up code. 
    Martin E. Liza  08/04/2021  Added plot_results function. 
'''
import numpy as np 
import matplotlib.pyplot as plt 
import IPython
import os 
# From my baseAnalysis class 
from baseAnalysis import Base_Analysis

# Line Class
class Line(Base_Analysis): 
    # Loads line data 
    flag_type   = 'line_data'
    flag_points = 'line_points' 

# Pre Process data 
    def pre_process(self, dataset_number, n_points): 
        time_axis  = self.working_data[dataset_number]['TIME'] 
        velocity_x = self.working_data[dataset_number]['U-X'] 
        [time_rows, spatial_columns] = np.shape(velocity_x) 
        # Make n_points times series 
        vector_series = np.arange(0, spatial_columns, n_points) 
        chosen_vel_x  = velocity_x[:,vector_series] 
        chosen_pos_x  = [ ]
        chosen_vel    = [ ]
        vel_reynolds_decomp = [ ]

        # Create fictitious position axis 
        for i in range(np.shape(chosen_vel_x)[1]): 
            mean_velocity = np.mean(chosen_vel_x[:,i])  
            pos_temp      = mean_velocity * time_axis 
            chosen_pos_x.append(pos_temp - np.min(pos_temp)) 
            vel_reynolds_decomp.append(chosen_vel_x[:,i] - mean_velocity) 
            chosen_vel.append(chosen_vel_x[:,i]) 

        # Dictionary  
        pre_process_dict = { 'reynolds_decomposition' : vel_reynolds_decomp,
                             'chosen_positions'       : vector_series, 
                             'raw_variable'           : chosen_vel, 
                             'time_radius_x'          : chosen_pos_x,
                             'velocity_x'             : chosen_vel_x } 

        return pre_process_dict 

# Calculates the mean 
    def mean_calculator(self, auto_correlation, length_scales, spe): 
        spe_keys         = list(spe.keys())
        correlation_keys = list(auto_correlation[spe_keys[0]].keys()) 
        scales_keys      = list(length_scales[spe_keys[0]].keys()) 
        spe_array        = len(spe[spe_keys[0]]) 
        corr_array       = len(auto_correlation[spe_keys[0]]
                                [correlation_keys[0]])

        # Return SPE 
        spe_return = [ ]
        for i in range(spe_array):  
            temp_val = [ ]
            for j in spe_keys: 
                temp_val.append(float(spe[j][i])) 
            spe_return.append(np.mean(temp_val)) 

        # Return Correlation 
        temp_corr_radius = [ ]
        temp_correlation = [ ]
        for i in range(corr_array): 
            temp_corr   = [ ] 
            temp_radius = [ ]
            for j in spe_keys:  
                temp_corr.append(auto_correlation[j]['correlation'][i])
                temp_radius.append(auto_correlation[j]
                        ['correlation_radius'][i])
            temp_corr_radius.append(np.mean(temp_radius))
            temp_correlation.append(np.mean(temp_corr))

        return_correlation = { 'correlation_radius' : temp_corr_radius,
                               'correlation'        : temp_correlation }
                
        # Length scales 
        temp_taylor_m   = [ ] 
        temp_integral_m = [ ] 
        temp_integral_k = [ ] 
        temp_cutoff_k   = [ ]
        temp_taylor_k   = [ ] 
        for i in spe_keys: 
            temp_taylor_m.append(float(length_scales[i]['taylor_m']))
            temp_integral_m.append(float(length_scales[i]['integral_m']))
            temp_integral_k.append(float(length_scales[i]['integral_k']))
            temp_cutoff_k.append(float(length_scales[i]['cutoff_k']))
            temp_taylor_k.append(float(length_scales[i]['taylor_k']))
        # Calculate the means 
        scales_return = { 'taylor_m'  : '{:.3e}'.format(
                                        np.mean(temp_taylor_m)),
                          'integral_m': '{:.3e}'.format(
                                        np.mean(temp_integral_m)),
                          'integral_k': '{:.3e}'.format(
                                        np.mean(temp_integral_k)),
                          'cutoff_k'  : '{:.3e}'.format(
                                        np.mean(temp_cutoff_k)),
                          'taylor_k'  : '{:.3e}'.format(
                                        np.mean(temp_taylor_k)) }

        return spe_return, scales_return, return_correlation 

        


 # Calculates x_axis
    def x_axis(self, dataset_number, sampling_rate=1):
        # Loading data 
        variable     = self.working_data[dataset_number]['U-X'] 
        sub_variable = self.sub_sampling(variable, sampling_rate) 
        # Load positions 
        # NOTE: x-axis is chosen to be x_3, if data needs to be collected 
        # in other directions, then this will need to change 
        data_location = self.location[dataset_number] 
        x1            = [data_location[0], data_location[3]] 
        x2            = [data_location[1], data_location[4]] 
        x3            = [data_location[2], data_location[5]] 
        x3_size       = np.shape(sub_variable)[1]  
        x_axis        = np.linspace(x3[0], x3[1], x3_size) 
        # Shift the axis, so it starts at 0.0 
        if (np.min(x_axis) < 0):  
            x_axis += np.abs(np.min(x_axis))
        if (np.min(x_axis) > 0): 
            x_axis -= np.min(x_axis)
        return x_axis 

# Calculates cross correlation 
    def cross_correlation(self, x_radius, fluctuation, cross_correlation_len): 
        # Auto Correlation First calculates it as a dictionary 
        auto_correlation_dict   = { }
        auto_correlation_radius = { }
        [rows, columns] = np.shape(fluctuation) 
        for i in range(rows): 
            auto_correlation_var = self.auto_correlation(x_radius, 
                                   fluctuation[i], cross_correlation_len) 
            auto_correlation_dict[f'{i}']   = auto_correlation_var['correlation']   
            auto_correlation_radius[f'{i}'] = auto_correlation_var['radius_x']
        # Calculates cross correlation 
        cor_columns = len(auto_correlation_dict['0']) 
        cor_temp = np.zeros(cor_columns) 
        rad_temp = np.zeros(cor_columns) 
        for i in range(cor_columns): 
            for j in range(rows):
                cor_temp[i] += auto_correlation_dict[f'{j}'][i] 
                rad_temp[i] += auto_correlation_radius[f'{j}'][i] 
        cross_correlation = cor_temp / rows 
        cross_radius      = rad_temp / rows  
        # Create a dictionary to return 
        crossCorrelation_dict = { 'radius_x'   : cross_radius, 
                                  'correlation': cross_correlation }
        return crossCorrelation_dict 

# Filter Decay 
    def line_filter_decay(self, dataset_number, dataset_variable, sampling_rate=1): 
        variable = self.working_data[dataset_number][dataset_variable] 
        sub_variable = self.sub_sampling(variable, sampling_rate)
        [rows, columns] = np.shape(sub_variable) 
        spe_dict = { } 
        # Sub sample and calculates the spe 
        for i in range(rows): 
            spe_dict[f'{i}'] = self.line_filter_decay(sub_variable[i])  
        # Shape of sampled data 
        sub_columns = len(spe_dict['0']) 
        spe_temp = np.zeros(sub_columns) 
        for j in range(sub_columns): 
            for k in range(rows): 
                spe_temp[j] += spe_dict[f'{k}'][j] 
        # Calculates the mean, and returns it  
        spe = spe_temp / rows  
        return spe  

# Plot results  
    def plot_results(self, dataset_number, dataset_variable, 
            cross_correlation=None, spe_var=[], sampling_rate=1, 
            correlation_len=50, inertial_shifting_factor=1e5, saving_path=None):    
    # Loading data 
        raw_variable = self.working_data[dataset_number][dataset_variable] 
        sub_variable = self.sub_sampling(raw_variable, sampling_rate)  
        sub_radius   = self.x_axis(dataset_number, dataset_variable, sampling_rate) 
        location     = self.location[dataset_number]
        x_1          = [location[0], location[3]]
        x_2          = [location[1], location[4]]
        x_3          = [location[2], location[5]]
        sampling_points = len(sub_variable)
        if (cross_correlation==None):
            fluctuation  = self.reynolds_decomposition(dataset_number, 
                           dataset_variable, sampling_rate)
            cross_correlation = self.cross_correlation(sub_radius, fluctuation, 
                                correlation_len) 
            corr_radius  = cross_correlation['radius_x']
            correlation  = cross_correlation['correlation']
        else:
            corr_radius  = cross_correlation['radius_x']
            correlation  = cross_correlation['correlation']
        # If there is not spe_var, calculate it 
        if (np.shape(spe_var)[0] == 0):
            spe = self.filter_decay(dataset_number, dataset_variable, sampling_rate)
        else: 
            spe = spe_var 
    # Length scales 
        scales_dict = self.length_scales(corr_radius, correlation) 
        integral = scales_dict['integral'] 
        taylor   = scales_dict['taylor'] 
    # Normalize length scales 
        dist_k  = np.array(range(len(spe))) 
        integral_k = dist_k[1] + (integral - sub_radius[1]) * (
                dist_k[-1] - dist_k[1]) / (sub_radius[-1] - sub_radius[1]) 
        taylor_k = dist_k[1] + (taylor - sub_radius[1]) * (
                dist_k[-1] - dist_k[1]) / (sub_radius[-1] - sub_radius[1]) 
    # Format order of magnitude 
        taylor       = "{:.3e}".format(taylor) 
        integral     = "{:.3e}".format(integral) 
        f_taylor_k   = "{:.3e}".format(taylor_k) 
        f_integral_k = "{:.3e}".format(integral_k)
    # Figure title  
        position     = f'x_1={x_1},\,x_2={x_2},\,x_3={x_3}' 
        var_string   = f'{dataset_number},\;{dataset_variable}'
        sampling_str = f'sampling\,rate = {sampling_rate}'
    # Creating plot  
        fig1 = plt.figure(1, figsize=(12, 6), dpi=95) 
        fig1.suptitle(f'${dataset_number},{dataset_variable}\;at\;{position}\;\;{sampling_str}$',
                fontsize=16) 
    # Correlation 
        plt.subplot(1,2,1)
        plt.plot(corr_radius, correlation, 'o-', 
                markerfacecolor='lightgray', linewidth='3', color='k') 
        if (dataset_variable == 'U-X'): 
            plt.text(corr_radius[-1], 1, 
             f'$L_{{11}}$= {integral}\n$\lambda_{{f}}$= {taylor}', 
             ha='right', va='top', 
             bbox=dict(facecolor='white', edgecolor='lightgray',
             boxstyle='round,pad=0.2'))
        elif (dataset_variable == 'U-Y'): 
            plt.text(corr_radius[-1], 1, 
             f'$L_{{22}}$= {integral}\n$\lambda_{{g}}$= {taylor}', 
             ha='right', va='top', 
             bbox=dict(facecolor='white', edgecolor='lightgray',
             boxstyle='round,pad=0.2'))
        else: 
            plt.text(corr_radius[-1], 1, 
             f'$L$= {integral}\n$\lambda$= {taylor}', 
             ha='right', va='top', 
             bbox=dict(facecolor='white', edgecolor='lightgray',
             boxstyle='round,pad=0.2'))
        plt.grid(linestyle='-.') 
        plt.ylabel('Correlation [ ]') 
        plt.xlabel('Radius [m]') 
    # SPE
        plt.subplot(1,2,2)
        #plt.axvline(x=dist_k[-1], color='g', linestyle='-.', label='$k_{max}$')
        #plt.axvline(x=dist_k[1], color='orange', linestyle='-.', label='$k_{min}$')
        if (dataset_variable == 'U-X'):
            plt.axvline(x=integral_k, color='r', linestyle='--',
                    linewidth='2', label=f'$L_{{11}}$= {f_integral_k}')
            plt.axvline(x=taylor_k, color='b', linestyle='--',
                    linewidth='2', label=f'$\lambda_{{f}}$= {f_taylor_k}')
        elif (dataset_variable == 'U-Y'):
            plt.axvline(x=integral_k, color='r', linestyle='--',
                    linewidth='2', label=f'$L_{{22}}$= {f_integral_k}')
            plt.axvline(x=taylor_k, color='b', linestyle='--',
                    linewidth='2', label=f'$\lambda_{{g}}$= {f_taylor_k}')
        else: 
            plt.axvline(x=integral_k, color='r', linestyle='--',
                    linewidth='2', label=f'$L$= {f_integral_k}')
            plt.axvline(x=taylor_k, color='b', linestyle='--',
                    linewidth='2', label=f'$\lambda$= {f_taylor_k}')
    # Plotting fiducial
        if (dataset_variable == 'P'):
            p_factor = -7/3 
        else: 
            p_factor = -5/3 
        wave_vector = np.array(range(10, int(np.floor(2 * sampling_points) / 15))) 
        fiducial = inertial_shifting_factor * wave_vector ** p_factor 
        plt.loglog(dist_k, spe, 'o-', markerfacecolor='lightgray', 
                   linewidth='3', color='k') 
        plt.plot(wave_vector, fiducial, linewidth='2')
        plt.grid(linestyle='-.') 
        plt.xlim(left=0.95) 
        plt.ylabel('E(k)')
        plt.xlabel('K-vector [1/m]') 
        plt.legend() 
    # Saving flag 
        if (saving_path == None):
            plt.show() 
        else: 
            plt.savefig(os.path.join(saving_path, 
                f'{dataset_number}_{dataset_variable}.png'))
            plt.close() 
