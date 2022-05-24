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
class Probe(Base_Analysis): 
# Loads 'probe_data'
    flag_type   = 'probe_data'
    flag_points = 'probe_points'

# Calculates the radius, might to rename this function to radius 
    def x_axis(self, dataset_number, sampling_flag=False):
    # Loads data 
        time_axis  = self.working_data[dataset_number]['TIME'] 
        velocity_x = self.working_data[dataset_number]['U-X'] 
    # Sub sample data 
        if sampling_flag:
            velocity_x  = self.sub_sampling(velocity_x)
            time_axis   = self.sub_sampling(time_axis)
    # Calculates radius 
        vel_x_len = len(velocity_x)
        mean      = np.mean(velocity_x) 
        time_mean = np.mean(np.diff(time_axis)) 
        delta_x   = np.abs(mean * time_mean)
        x_mean    = np.linspace(0, (delta_x * vel_x_len), vel_x_len)
        x_true_axis  = (time_axis * velocity_x) - (time_axis[0] * velocity_x[0])
        return_dict = { 'x_axis' : x_mean, 
                        'x_true' : x_true_axis }
        return return_dict

# String Title 
    def plot_title(self, dataset_number, dataset_variable):
    # Probe Locations  
        x_1               = self.location[dataset_number][0]  
        x_2               = self.location[dataset_number][1]  
        x_3               = self.location[dataset_number][2]  
        iterations        = self.working_data[dataset_number]['ITER'] 
        sub_iterations    = self.sub_sampling(iterations)
        sampling_rate     = str(int(sub_iterations[1] - sub_iterations[0])) 
        location     = f'x_1={x_1}, x_2={x_2}, x_3={x_3}'
        var_string   = f'{dataset_number}, {dataset_variable}'  
        sampling_str = f'sampling rate = {sampling_rate}' 
        title_str    = f'{var_string} at [{location}], sampling rate={sampling_rate}'
        return title_str 

# Define a constant cutoff
    def const_cutoff_k(self, dataset_number, probe_radius, correlation_lag=50): 
        # Cutoff constant 
        u_x    = self.working_data[dataset_number]['U-X'] 
        u_x    = self.sub_sampling(u_x)
        fluct  = self.reynolds_decomposition(u_x) 
        spe    = self.filter_decay(u_x) 
        auto_corr = self.auto_correlation(probe_radius, fluct, 
                auto_correlation_len=correlation_lag) 
        length_scales = self.length_scales(
                             auto_corr['correlation_radius'], 
                             auto_corr['correlation'],  
                             fluct, spe) 
        return length_scales['cutoff_k'] 
    
# Plot results  
    def plot_results(self, dataset_number, dataset_variable, radius, variable,  
                spe, auto_correlation_dict, length_scales_dict, moments_str_dict, 
                inertial_shifting_factor=1e6, histogram_bins=50, saving_path=None):    
    # Load dictionaries data 
        corr_radius = auto_correlation_dict['correlation_radius']
        correlation = auto_correlation_dict['correlation']
        integral_m  = length_scales_dict['integral_m']  
        taylor_m    = length_scales_dict['taylor_m']  
        integral_k  = length_scales_dict['integral_k']  
        taylor_k    = length_scales_dict['taylor_k']  
        iteration   = self.working_data[dataset_number]['ITER']
        dist_k      = range(len(spe))  
        sampling_points = iteration[1] - iteration[0] 
    # Split Plots  
        fig1 = plt.figure(2, figsize=(18,8), dpi=95) 
        gs  = gridspec.GridSpec(ncols=3, nrows=2) 
        gs.update(wspace=0.3, hspace=0.2) 
        ax1 = plt.subplot(gs[0,:1]) 
        ax2 = plt.subplot(gs[0,1:2]) 
        ax3 = plt.subplot(gs[0,2:3]) 
        ax4 = plt.subplot(gs[1,0:]) 
    # Figure title 
        title_str = self.plot_title(dataset_number, dataset_variable) 
        fig1.suptitle(f'{title_str}', fontsize=16) 
    # Intersection point
        taylor_rad  = np.linspace(float(taylor_m), 0) 
        shift_b     = np.linspace(1, float(taylor_m)) 
        parabola    = -(taylor_rad + shift_b)**2 + 1 
        try:
            max_indx = np.where(parabola < np.min(correlation))[0][-1] # index with the last negative value 
        except Exception:
            max_indx = 0

        taylor_rad = taylor_rad[max_indx:] 
        parabola   = parabola[max_indx:]
    # Plotting Correlation 
        ax1.plot(corr_radius, correlation, 'o-', 
                markerfacecolor='lightgray', 
                 linewidth='3', color='k') 
        ax1.plot(taylor_rad, parabola, '--', linewidth='3', color='b') 

        if (dataset_variable == 'U-X'): 
            ax1.text(corr_radius[-1], 1, 
                    f'$L_{{11}}$= {integral_m}\n$\lambda_{{f}}$= {taylor_m}', 
                    ha='right', va='top', 
                    bbox=dict(facecolor='white', 
                        edgecolor='lightgray', boxstyle='round,pad=0.2')) 
        elif (dataset_variable == 'U-Y'): 
            ax1.text(corr_radius[-1], 1, 
                    f'$L_{{22}}$= {integral_m}\n$\lambda_{{g}}$= {taylor_m}', 
                    ha='right', va='top', bbox=dict(facecolor='white', 
                        edgecolor='lightgray', boxstyle='round,pad=0.2')) 
        else: 
            ax1.text(corr_radius[-1], 1, 
                    f'$L$= {integral_m}\n$\lambda$= {taylor_m}',
                    ha='right', va='top', 
                    bbox=dict(facecolor='white', 
                        edgecolor='lightgray', boxstyle='round,pad=0.2')) 
        ax1.grid(linestyle='-.') 
        ax1.set_ylabel('Correlation [ ]') 
        ax1.set_xlabel('Radius [m]') 
    # Plotting Energy Spectrum 
        if (dataset_variable == 'U-X'): 
            ax2.axvline(x=float(integral_k), color='r', 
                    linestyle='--', linewidth='2', 
                    label=f'$L_{{11}}$= {integral_k}') 
            ax2.axvline(x=float(taylor_k), color='b', 
                    linestyle='--', linewidth='2', 
                    label=f'$\lambda_{{f}}$= {taylor_k}') 
        elif (dataset_variable == 'U-Y'): 
            ax2.axvline(x=float(integral_k), color='r', 
                    linestyle='--', linewidth='2', 
                    label=f'$L_{{22}}$= {integral_k}') 
            ax2.axvline(x=float(taylor_k), color='b', 
                    linestyle='--', linewidth='2', 
                    label=f'$\lambda_{{g}}$= {taylor_k}') 
        else: 
            ax2.axvline(x=float(integral_k), color='r', 
                    linestyle='--', linewidth='2', 
                    label=f'$L$= {integral_k}') 
            ax2.axvline(x=float(taylor_k), color='b', 
                    linestyle='--', linewidth='2', 
                    label=f'$\lambda$= {taylor_k}') 
                # Plotting fiducial 
        if (dataset_variable == 'P'): 
            p_factor = -7/3 
        else:
            p_factor = - 5/3 
        wave_vector = np.array(range(10, int(np.floor(2 * len(spe)) / 15))) 
        fiducial = inertial_shifting_factor * wave_vector ** p_factor 
        ax2.loglog(dist_k, spe, 'o-', markerfacecolor='lightgray', 
                linewidth='3', color='k') 
        ax2.plot(wave_vector, fiducial, linewidth='3')  
        ax2.grid(linestyle='-.') 
        ax2.set_xlim(left=0.9) 
        ax2.set_ylabel('E(k)') 
        ax2.set_xlabel('K-vector [1/m]') 
        ax2.legend() 

    # Plotting histogram 
        ax3.hist(variable, label=moments_str_dict['full_str'], 
                bins=histogram_bins, density=True, 
                color='lightcyan', edgecolor='k') 
        density, bins_loc = np.histogram(variable, density=True) 
        sb.kdeplot(variable, bw_adjust=1.3, fill=False, linewidth=2, 
                    color='crimson', ax=ax3) 
        ax3.grid(linestyle='-.') 
        ax3.set_xlabel(f'{dataset_variable}')
        ax3.set_ylabel('Density')
        ax3.legend(bbox_to_anchor=(0.8,0.7), loc='lower left', 
                facecolor='white', framealpha=1,
                handlelength=0, handletextpad=0, fancybox=True) 

        # Plotting original signal  
        ax4.plot(radius, variable, color='k', linewidth='3')
        ax4.axhline(y=np.mean(variable), color='darkorange', 
                    linestyle='--', linewidth='2',
                    label=f'$<{dataset_variable}>$') 
        ax4.set_xlim([radius[0]-1e-4, radius[-1]+1e-4])
        ax4.set_xlabel('Radius [m]')
        ax4.set_ylabel(f'{dataset_variable}') 
        ax4.grid(linestyle='-.') 
        ax4.legend() 
        # If saving flag, save 
        if (saving_path == None):
            plt.show() 
        else:
            plt.savefig(os.path.join(saving_path,
            f'results_{dataset_number}_{dataset_variable}.png'))
            plt.close() 

# Plot scatter filters  
    def plot_scatter(self, dataset_number, variable_x, variable_y, 
            probe_radius, const_cutoff_k, 
            correlation_lag=50, saving_path=None): 
        data_struct = { }
        variables   = [variable_x, variable_y] 

    # Title 
        # Probe Locations  
        x_1 = self.location[dataset_number][0]  
        x_2 = self.location[dataset_number][1]  
        x_3 = self.location[dataset_number][2]  
        # Sampling rate
        iterations     = self.working_data[dataset_number]['ITER'] 
        sub_iterations = self.sub_sampling(iterations)
        sampling_rate  = str(int(sub_iterations[1] - sub_iterations[0])) 
        location     = f'x_1={x_1},\,x_2={x_2},\,x_3={x_3}' 
        var_string   = f'{dataset_number},\;{variable_y}\;Vs.\;{variable_x}'  
        sampling_str = f'sampling\,rate = {sampling_rate}' 

    # Create structures for plotting  
        for i in variables: 
            var    = self.working_data[dataset_number][i] 
            var    = self.sub_sampling(var)
            data_struct[f'{i}_boxcar']   = self.boxcar_filter(probe_radius, 
                                           var, const_cutoff_k)
            data_struct[f'{i}_legendre'] = self.legendre_interpolation(
                                            data_struct[f'{i}_boxcar'])
        # Clean data, remove elements that no need to be plot  
            data_struct[f'{i}_boxcar'].pop('window_size') 
            data_struct[f'{i}_boxcar'].pop('full_radius') 
            data_struct[f'{i}_legendre'].pop('delta_radius') 
            data_struct[f'{i}_boxcar'].pop('les_radius') 
    # Generate plots 
        var_boxcar_keys = list(data_struct[f'{variables[0]}_boxcar'].keys())
        var_legend_keys = list(data_struct[f'{variables[0]}_legendre'].keys())
        var_keys        = var_boxcar_keys + var_legend_keys 
        # Creating figure for plot 
        fig = plt.figure(figsize=(17, 8)) 
        fig.subplots_adjust(hspace=0.4, wspace=0.4) 
        fig.suptitle(f'${var_string}\;at\;[{location}]\;\;{sampling_str}$',
                      fontsize=16) 
        for n in range(len(var_keys)):
            filter_name = 'boxcar'
            label_name  = var_keys[n] 
            label_type  = label_name.split('_')[0]
            ax = fig.add_subplot(2,3,n+1) 
            if (n >= 4):
                filter_name = 'legendre'
            ax.scatter(data_struct[f'{variables[0]}_{filter_name}'][label_name], 
                       data_struct[f'{variables[1]}_{filter_name}'][label_name],
                       facecolor='lightgray', edgecolors='black')   
            ax.set_xlabel(f'{variables[0]} {label_type}')
            ax.set_ylabel(f'{variables[1]} {label_type}')
            #ax.axis('equal') 
            #ax.set_aspect('equal', 'box') 
        # Saving flag 
        if (saving_path == None):
            plt.show() 
        else:
            plt.savefig(os.path.join(saving_path, 
                f'scatter_{dataset_number}_{variables[0]}_{variables[1]}.png'))
            plt.close() 
