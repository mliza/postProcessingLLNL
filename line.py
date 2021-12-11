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
    def temporal_data(self, dataset_number, 
            dataset_variable, n_points, auto_correlation_len=50): 
        # Loads data 
        data_var   = self.working_data[dataset_number][dataset_variable].T
        vel_x      = self.working_data[dataset_number]['U-X']
        time_axis  = self.working_data[dataset_number]['TIME'] 
        [time_rows, spatial_columns] = np.shape(vel_x) 
        # Make n_points times series 
        spatial_sampling  = np.arange(0, spatial_columns, n_points) 
        # Create fictitious time series data  
        return_dict = { }
        for i in spatial_sampling: 
            radius_x = time_axis * np.mean(vel_x[i]) 
            radius_x -= np.min(radius_x) 
            data_out = self.data_process(data_var[i], radius_x, 
                       auto_correlation_len)
            return_dict[i] = data_out 
            return_dict[i]['radius'] = radius_x 
        return return_dict  

# Spatial Data 
    def spatial_data(self, dataset_number,
            dataset_variable, n_points, auto_correlation_len=50): 
        # Loads data 
        data_var   = self.working_data[dataset_number][dataset_variable]
        [time_rows, spatial_columns] = np.shape(data_var) 
        radius_z = self.z_axis(dataset_number)  
        # Make n_points times series 
        temporal_sampling  = np.arange(0, time_rows, n_points) 
        return_dict = { }
        # Calculate absolute length scales 
        for i in temporal_sampling:  
            data_out = self.data_process(data_var[i], radius_z, 
                       auto_correlation_len)       
            return_dict[i] = data_out 
            return_dict[i]['radius'] = radius_z 
        return return_dict  

# Create constant cutoff_scale 
    def const_cutoff_k(self, dict_variable): 
        scales = self.length_scales( dict_variable['correlation']['correlation_radius'],
                                     dict_variable['correlation']['correlation'],
                                     dict_variable['fluctuation'],
                                     dict_variable['spe'] )
        return scales['cutoff_k'] 

# Calculates the mean 
    def data_cruncher(self, dict_in): 
        sampling_location = list(dict_in.keys()) 
        data_elements     = list(dict_in[sampling_location[0]].keys())  
        radius_len        = len(dict_in[sampling_location[0]]['radius'])
        correlation_keys  = list(dict_in[sampling_location[0]]
                                ['correlation'].keys())
        corr_len          = len(dict_in[sampling_location[0]]
                                ['correlation']['correlation'])
        # Radius and fluctuation
        radius      = [ ]
        fluctuation = [ ]
        for i in range(radius_len): 
            radius_temp      = [ ] 
            fluctuation_temp = [ ]
            for j in sampling_location: 
                radius_temp.append(dict_in[j]['radius'][i])
                fluctuation_temp.append(dict_in[j]['fluctuation'][i])
            # Append the mean of each element 
            radius.append(np.mean(radius_temp)) 
            fluctuation.append(np.mean(fluctuation_temp))
        # Correlation 
        correlation_radius = [ ]
        correlation        = [ ]
        for i in range(corr_len):
            corr_temp     = [ ] 
            corr_rad_temp = [ ]
            for j in sampling_location:
                corr_rad_temp.append(dict_in[j]['correlation']
                            ['correlation_radius'][i])
                corr_temp.append(dict_in[j]['correlation']
                            ['correlation'][i])
            # Append the mean of each element 
            correlation_radius.append(np.mean(corr_rad_temp))  
            correlation.append(np.mean(corr_temp)) 
        # Create dictionary 
        corr_temp_dict = { 'correlation_radius' : np.asarray(correlation_radius),
                           'correlation'        : np.asarray(correlation) }
        # Spe 
        spe_return = [ ]
        for i in range(len(dict_in[0]['spe'])):
            temp_spe = [ ] 
            for j in sampling_location: 
                temp_spe.append(dict_in[j]['spe'][i]) 
            spe_return.append(np.mean(temp_spe)) 

        # Create Return Dictionary 
        return_crunched_dat = { 'radius'            : np.asarray(radius), 
                                'fluctuation'       : np.asarray(fluctuation), 
                                'correlation'       : corr_temp_dict, 
                                'spe'               : np.asarray(spe_return), 
                                'sampling_location' : sampling_location }
        return return_crunched_dat  

# Calculate all filters 
    def filters(self, dict_in, window_size):
        sampling_elements = list(dict_in.keys()) 
        boxcar_dict       = { }
        legendre_dict     = { } 

        # Calculates boxcar filter and legendre interpolation 
        for i in sampling_elements: 
            boxcar           = self.boxcar_filter(dict_in[i]['radius'],
                                dict_in[i]['variable'], window_size) 
            legendre         = self.legendre_interpolation(boxcar) 
            boxcar_dict[i]   = boxcar  
            legendre_dict[i] = legendre      

        # Parameters for data crunching 
        legendre_keys = list(legendre.keys()) 
        boxcar_keys   = list(boxcar.keys())
        boxcar_keys.remove('window_size') 
        boxcar_return   = { } 
        legendre_return = { }

        # Iterates through boxcar_keys  
        for i in boxcar_keys:
            boxcar_vec   = [ ]
            # Iterates through boxcar_keys length  
            boxcar_len = range(len(boxcar[i])) 
            for j in boxcar_len:
                temp_boxcar = [ ]
                # Iterates through sampling elements
                for k in sampling_elements:
                    temp_boxcar.append(boxcar_dict[k][i][j]) 
                boxcar_vec.append(np.mean(temp_boxcar))
            boxcar_return[i] = boxcar_vec 

        # Iterates through legendre_keys  
        for i in legendre_keys:
            legendre_vec = [ ] 
            # Iterates through legendre_keys length  
            legendre_len = range(len(legendre[i])) 
            for j in legendre_len:
                temp_legendre = [ ]
                # Iterates through sampling elements
                for k in sampling_elements:
                    temp_legendre.append(legendre_dict[k][i][j]) 
                legendre_vec.append(np.mean(temp_legendre))
            legendre_return[i] = legendre_vec 

        # Return dictionary 
        filter_return = { 'boxcar'   : boxcar_return, 
                          'legendre' : legendre_return }

        return filter_return 


                

        
                
 # Calculates z_axis
    def z_axis(self, dataset_number, sampling_rate=1):
        # Loading data 
        variable     = self.working_data[dataset_number]['U-Z'] 
        sub_variable = self.sub_sampling(variable, sampling_rate) 
        # Load positions 
        # NOTE: x-axis is chosen to be x_3, if data needs to be collected 
        # in other directions, then this will need to change 
        data_location = self.location[dataset_number] 
        x1            = [data_location[0], data_location[3]] 
        x2            = [data_location[1], data_location[4]] 
        x3            = [data_location[2], data_location[5]] 
        x3_size       = np.shape(sub_variable)[1]  
        z_axis        = np.linspace(x3[0], x3[1], x3_size) 
        # Shift the axis, so it starts at 0.0 
        if (np.min(z_axis) < 0):  
            z_axis += np.abs(np.min(z_axis))
        if (np.min(z_axis) > 0): 
            z_axis -= np.min(z_axis)
        return z_axis 

# Plot correlation + spe (temporal)  
    def plot_correlation_spe(self, temporal_dict, spatial_dict, dataset, variable,
            time_sub_sampling, spatial_sub_sampling, saving_path=None): 
        # Load data for easy plotting 
        # Spatial Data 
        spatial_taylor_m   = spatial_dict['length_scales']['taylor_m'] 
        spatial_integral_m = spatial_dict['length_scales']['integral_m'] 
        spatial_integral_k = spatial_dict['length_scales']['integral_k'] 
        spatial_taylor_k   = spatial_dict['length_scales']['taylor_k'] 
        spatial_corr_rad   = spatial_dict['correlation']['correlation_radius'] 
        spatial_corr       = spatial_dict['correlation']['correlation'] 
        spatial_spe        = spatial_dict['spe'] 

        # Temporal Data 
        temporal_taylor_m   = temporal_dict['length_scales']['taylor_m'] 
        temporal_integral_m = temporal_dict['length_scales']['integral_m'] 
        temporal_integral_k = temporal_dict['length_scales']['integral_k'] 
        temporal_taylor_k   = temporal_dict['length_scales']['taylor_k'] 
        temporal_corr_rad   = temporal_dict['correlation']['correlation_radius'] 
        temporal_corr       = temporal_dict['correlation']['correlation'] 
        temporal_spe        = temporal_dict['spe'] 

        # Figure 
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        fig.set_size_inches(15, 8) 

        # Load all data 
        fig.suptitle(f'{dataset}, {variable}, time sampling={time_sub_sampling}, spatial sampling={spatial_sub_sampling}') 

         # Plotting Spatial Correlation 
        ax1.plot(spatial_corr_rad, spatial_corr, linewidth='2', color='k', 
                label=f'L={spatial_integral_m}, $\lambda$={spatial_taylor_m}')
        ax1.set_ylabel('Correlation, spatial series')
        ax1.grid('-.')
        ax1.legend(handlelength=0, handletextpad=0, fancybox=True) 
        ax1.set_xlim(left=0.0) 

         # Plotting Spatial Energy Spectrum
        ax2.loglog(spatial_spe, linewidth='3', color='k') 
        ax2.axvline(x=float(spatial_integral_k), color='r', linestyle='--', 
                linewidth='2', label=f'$L_k$={spatial_integral_k}') 
        ax2.axvline(x=float(spatial_taylor_k), color='b', linestyle='--', 
                linewidth='2', label=f'$\lambda_k$={spatial_taylor_k}') 
        ax2.set_ylabel('Energy Spectrum, spatial series')
        ax2.grid('-.')
        ax2.set_xlim(left=0.9) 
        ax2.legend() 

        # Plotting Time Correlation  
        ax3.plot(temporal_corr_rad, temporal_corr, linewidth='3', color='k',
            label=f'L={temporal_integral_m}, $\lambda$={temporal_taylor_m}')
        ax3.set_ylabel('Correlation, temporal series')
        ax3.set_xlabel('Radius [m]')
        ax3.grid('-.')
        ax3.legend(handlelength=0, handletextpad=0, fancybox=True) 
        ax3.set_xlim(left=0.0) 

         # Plotting time Energy Spectrum
        ax4.loglog(temporal_spe, linewidth='3', color='k') 
        ax4.axvline(x=float(temporal_integral_k), color='r', linestyle='--', 
                linewidth='2', label=f'$L_k$={temporal_integral_k}') 
        ax4.axvline(x=float(temporal_taylor_k), color='b', linestyle='--', 
                linewidth='2', label=f'$\lambda_k$={temporal_taylor_k}') 
        ax4.set_ylabel('Energy Spectrum, temporal series')
        ax4.set_xlabel('k-vector [1/m]')
        ax4.grid('-.')
        ax4.set_xlim(left=0.9) 
        ax4.legend() 
        if (saving_path == None):
            plt.show()
        else: 
            plt.savefig(os.path.join(saving_path, f'{dataset}_{variable}.png')) 
            plt.close() 

# String Title 
    def plot_title(self, dataset_number, dataset_variable):
    # Probe Locations  
        loc        = self.location[dataset_number] 
        location   = f'x_1=[{loc[0]}, {loc[3]}], x_2=[{loc[1]}, {loc[4]}], x_3=[{loc[2]}, {loc[5]}]'
        var_string = f'{dataset_number}, {dataset_variable}'  
        title_str  = f'{var_string} at {location}'
        return title_str 
