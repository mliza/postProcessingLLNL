#!/usr/local/bin/python3.9
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
from scipy import signal 
from scipy.fft import fft, fftshift 
import scipy 
# From my baseAnalysis class 
from baseAnalysis import Base_Analysis

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
        if (sampling_flag == True):
            velocity_x  = self.sub_sampling(velocity_x)
            time_axis   = self.sub_sampling(time_axis)
    # Calculates radius 
        vel_x_len = len(velocity_x)
        mean      = np.mean(velocity_x) 
        delta_x   = np.abs(mean * (time_axis[2] - time_axis[1]))
        x_axis    = np.linspace(0, (delta_x * vel_x_len), vel_x_len)
        x_axis    = (time_axis * velocity_x) - (time_axis[0] * velocity_x[0])
        return x_axis

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
    # Plotting Correlation 
        ax1.plot(corr_radius, correlation, 'o-', 
                markerfacecolor='lightgray', 
                 linewidth='3', color='k') 
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

# Plot Legendre filter  
    def plot_legendre(self, dataset_number, dataset_variable, 
            boxcar_dictionary, legendre_dictionary, saving_path=None): 
    # Loading data  
        les_variable      = boxcar_dictionary['les_variable'] 
        les_radius        = boxcar_dictionary['les_radius']
        filtered_radius   = boxcar_dictionary['full_radius'] 
        filtered_variable = boxcar_dictionary['filtered_variable'] 
        window_size       = boxcar_dictionary['window_size']
        legendre_variable = legendre_dictionary['delta_variable']
        legendre_radius   = legendre_dictionary['delta_radius'] 
        iterations        = self.working_data[dataset_number]['ITER'] 
    # Title 
        title_str = self.plot_title(dataset_number, dataset_variable) 
    # Generate plots 
        fig, axs = plt.subplots(3, figsize=(15,7)) 
        fig.suptitle(f'{title_str}, window_size={window_size}', fontsize=16)
        # Plot 1st figure 
        axs[0].plot(filtered_radius, filtered_variable, label='filtered variable',
                linewidth=3, color='k') 
        axs[1].plot(les_radius, les_variable, 'o-', markerfacecolor='lightgray',
                    drawstyle='steps-mid', markersize=6, linewidth='3', 
                    color='k', label='boxcar sampling')
        axs[2].plot(legendre_radius, legendre_variable, 'o-', 
                    markerfacecolor='lightgray', 
                    drawstyle='steps-mid', markersize=6, linewidth='3', 
                    color='k', label='legendre')
        axs[0].grid('-.')
        axs[0].margins(x=0) 
        axs[0].set_ylabel(f'{dataset_variable}, filtered')
        axs[0].legend() 
        # Plot 1nd figure 
        axs[1].grid('-.')
        axs[1].margins(x=0) 
        axs[1].set_ylabel(f'{dataset_variable}, les')
        axs[1].legend() 
        # Plot 2nd figure 
        axs[2].grid('-.')
        axs[2].margins(x=0) 
        axs[2].set_xlabel('Radius [m]')
        axs[2].set_ylabel(f'{dataset_variable}, delta')
        axs[2].legend() 
    #  Save results  
        if (saving_path == None):
            plt.show() 
        else:
            plt.savefig(os.path.join(saving_path,
            f'legendre_{dataset_number}_{dataset_variable}.png'))
            plt.close() 

# Plot boxcar filter  
    def plot_boxcar(self, dataset_number, dataset_variable,
            boxcar_dictionary, moments_str_dict, saving_path=None): 
    # Loading data 
        resolved_eddies = boxcar_dictionary['filtered_variable'] 
        sgs_eddies      = boxcar_dictionary['sgs_variable']
        full_radius     = boxcar_dictionary['full_radius'] 
        full_variable   = boxcar_dictionary['full_variable']
        window_size     = boxcar_dictionary['window_size']
    # Plotting 
        fig = plt.figure() 
        plt.gcf().set_size_inches(17,8) 
        fig.tight_layout()  
        gs  = gridspec.GridSpec(3, 2, width_ratios=[5,2])
        ax_0 = plt.subplot(gs[0]) 
        ax_1 = plt.subplot(gs[1]) 
        ax_2 = plt.subplot(gs[2]) 
        ax_3 = plt.subplot(gs[3]) 
        ax_4 = plt.subplot(gs[4]) 
        ax_5 = plt.subplot(gs[5]) 
    # Figure title 
        title_str    = self.plot_title(dataset_number, dataset_variable) 
        fig.suptitle(f'{title_str}, window size = {window_size}', fontsize=16)
    # Create moment strings 
        full_str     = moments_str_dict['full_str'] 
        resolved_str = moments_str_dict['resolved_str'] 
        sgs_str      = moments_str_dict['sgs_str'] 
    # ax_0
        ax_0.plot(full_radius, full_variable, color='k', linewidth='3')
        ax_0.set_ylabel(f'{dataset_variable}, full')
        ax_0.grid('-.')
        ax_0.margins(x=0) 
    # Ax 1
        ax_1.hist(full_variable, label=full_str,  bins=50, density=True, 
                  color='lightcyan', edgecolor='k')
        sb.kdeplot(full_variable, bw_adjust=1.3, fill=False, 
                    linewidth=2, color='crimson', ax=ax_1) 
        ax_1.legend(bbox_to_anchor=(0.72,0.55), loc='lower left', 
                facecolor='white', framealpha=1,
                handlelength=0, handletextpad=0, fancybox=True) 
        ax_1.set_ylabel('Probability') 
    # Ax 2
        ax_2.plot(full_radius, resolved_eddies, color='k', linewidth='3')
        ax_2.set_ylabel(f'{dataset_variable}, filtered')
        ax_2.grid('-.')
        ax_2.margins(x=0) 
    # Ax 3
        ax_3.hist(resolved_eddies, label=resolved_str, bins=50, density=True, 
                  color='lightcyan', edgecolor='k')
        sb.kdeplot(resolved_eddies, bw_adjust=1.3, fill=False, 
                    linewidth=2, color='crimson', ax=ax_3) 
        ax_3.legend(bbox_to_anchor=(0.72,0.55), loc='lower left', 
                facecolor='white', framealpha=1,
                handlelength=0, handletextpad=0, fancybox=True) 
        ax_3.set_ylabel('Probability') 
    # Ax 4
        ax_4.plot(full_radius, sgs_eddies, color='k', linewidth='3')
        ax_4.set_ylabel(f'{dataset_variable}, sgs')
        ax_4.set_xlabel(f'Radius [m]')
        ax_4.grid('-.')
        ax_4.margins(x=0) 
    # Ax 5
        ax_5.hist(sgs_eddies, label=sgs_str,  bins=50, density=True, 
                  color='lightcyan', edgecolor='k')
        sb.kdeplot(sgs_eddies, bw_adjust=1.3, fill=False, 
                    linewidth=2, color='crimson', ax=ax_5) 
        ax_5.legend(bbox_to_anchor=(0.72,0.55), loc='lower left', 
                facecolor='white', framealpha=1,
                handlelength=0, handletextpad=0, fancybox=True) 
        ax_5.set_ylabel('Probability') 
        ax_5.set_xlabel(f'{dataset_variable}') 
        #fig.supxlabel('Radius [m]') 
        if (saving_path == None):
            plt.show() 
        else:
            plt.savefig(os.path.join(saving_path,
            f'boxcar_{dataset_number}_{dataset_variable}.png'))
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
            ax.axis('equal') 
            #ax.set_aspect('equal', 'box') 
        # Saving flag 
        if (saving_path == None):
            plt.show() 
        else:
            plt.savefig(os.path.join(saving_path, 
                f'scatter_{dataset_number}_{variables[0]}_{variables[1]}.png'))
            plt.close() 


