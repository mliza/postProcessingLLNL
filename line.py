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

# Temporal Data 
    def temporal_data(self, dataset_number, n_points, auto_correlation_len=50): 
        time_axis  = self.working_data[dataset_number]['TIME'] 
        velocity_x = self.working_data[dataset_number]['U-X'] 
        [time_rows, spatial_columns] = np.shape(velocity_x) 
        # Make n_points times series 
        spatial_sampling  = np.arange(0, spatial_columns, n_points) 
        # Create fictitious time series data  
        return_dict = { }
        for i in spatial_sampling: 
            return_dict[i] = { }
            mean_velocity = np.mean(velocity_x[i])  
            time_pos_temp = mean_velocity * time_axis 
            # Starts the time position at 0 
            position      = time_pos_temp - np.min(time_pos_temp)  
            fluctuation   = velocity_x[i] - mean_velocity
            correlation   = self.auto_correlation(position, fluctuation, 
                            auto_correlation_len)
            spe           = self.filter_decay(velocity_x[i]) 
            length_scale  = self.length_scales(correlation['correlation_radius'],
                            correlation['correlation'], fluctuation, spe)
            # Return Dictionary  
            return_dict[i]['time_radius']        = position 
            return_dict[i]['time_variable']      = velocity_x[i] 
            return_dict[i]['time_fluctuation']   = fluctuation 
            return_dict[i]['time_correlation']   = correlation
            return_dict[i]['time_spe']           = spe
            return_dict[i]['time_length_scales'] = length_scale

        return return_dict  

# Spatial Data 
    def spatial_data(self, dataset_number, n_points, auto_correlation_len=50): 
        time_axis  = self.working_data[dataset_number]['TIME'] 
        velocity_x = self.working_data[dataset_number]['U-X'] 
        [time_rows, spatial_columns] = np.shape(velocity_x) 
        spatial_location = self.location[dataset_number]
        radius = np.linspace(spatial_location[2], spatial_location[5],
                             spatial_columns) 
        # Radius_Z 
        radius_z = radius - np.min(radius)
        # Make n_points times series 
        temporal_sampling  = np.arange(0, time_rows, n_points) 
        return_dict = { }
        for i in temporal_sampling:  
            return_dict[i] = { }
            mean_velocity = np.mean(velocity_x[i])  
            fluctuation   = velocity_x[i] - mean_velocity
            correlation   = self.auto_correlation(radius_z, fluctuation, 
                            auto_correlation_len)
            spe           = self.filter_decay(velocity_x[i]) 
            length_scale  = self.length_scales(correlation['correlation_radius'],
                            correlation['correlation'], fluctuation, spe)
            # Return Dictionary  
            return_dict[i]['space_radius']        = radius_z  
            return_dict[i]['space_variable']      = velocity_x[i] 
            return_dict[i]['space_fluctuation']   = fluctuation 
            return_dict[i]['space_correlation']   = correlation
            return_dict[i]['space_spe']           = spe
            return_dict[i]['space_length_scales'] = length_scale

        return return_dict  

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

