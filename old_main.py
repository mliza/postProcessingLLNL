#!/usr/local/bin/python3.9 
import numpy as np 
import pickle 
import os 
import IPython 
import math 
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec 
import seaborn as sb 
# My files 
import baseAnalysis
import probe 
#import line 
import new_line as line 

line_flag  = False    
probe_flag = False  
new_data_flag = False  
scatter_probe_flag = False  
probe_sampling_rate   = 1 
sub_sampling_flag = False
probe_correlation_lag = 100 
line_sampling_rate    = [1,1] 
line_correlation_lag  = 40 

# Paths 
pickle_path = '/Users/martin/Documents/Research/UoA/Projects/LLNL/data/data_5/pickle' 
save_path     = '/Users/martin/Desktop/workingFiles/results' 
probe_save    = os.path.join(save_path, 'probe', 'total')
boxcar_save   = os.path.join(save_path, 'probe', 'boxcar')
scatter_save  = os.path.join(save_path, 'probe', 'scatter')
legendre_save = os.path.join(save_path,'probe', 'legendre') 
line_save     = os.path.join(save_path, 'line')

# Loading class's instances  
base  = baseAnalysis.Base_Analysis(pickle_path)
probe = probe.Probe(pickle_path, sampling_rate=probe_sampling_rate) 
line  = line.Line(pickle_path, sampling_rate=probe_sampling_rate) 


# Keys 
probe_keys = list(probe.location.keys()) 
line_keys  = list(line.location.keys()) 
variables  = ['U-X', 'U-Y', 'U-Z', 'P', 'T', 'RHO', 'RHOE', 'GRADRHOMAG', 'DIL', 'VORTMAG', 'P-DIL', 'RHO-DIL'] 
variables  = ['U-X', 'U-Y', 'U-Z', 'P', 'T', 'RHO','DIL', 'P-DIL'] 

# Line Pre-Processing loop, returns the mean_cutoff for each  
time_sub_sampling    = 500 
spatial_sub_sampling = 120 
dict_pre_process = { } 
mean_cutoff_k    = { }
testing_len = np.shape(line.working_data['l0']['U-X'])[0] 
testing_vec = np.arange(0, testing_len, 1000)  

# Time Series 
for i in line_keys: 
    pre_process_dict  = line.pre_process(i, n_points=spatial_sub_sampling ) 
    radius_x          = pre_process_dict['time_radius_x'] 
    fluctuation       = pre_process_dict['reynolds_decomposition'] 
    variable          = pre_process_dict['velocity_x'] 
    dict_pre_process[i] = { }
    dict_pre_process[i]['time_correlation'] = { }
    dict_pre_process[i]['length_scales']    = { }
    dict_pre_process[i]['length_scales_k']  = { }
    dict_pre_process[i]['spe']              = { }
    mean_cutoff_k[i] = [ ] 
    cutoff_scale     = [ ] 
    for j in range(np.shape(radius_x)[0]): 
        auto_correlation = line.auto_correlation(radius_x[j], 
                        fluctuation[j], 60) 
        spe = line.filter_decay(variable[:,j])  
        length_scales = line.length_scales(auto_correlation['correlation_radius'], 
                                   auto_correlation['correlation'], 
                                   fluctuation[j], spe) 
        IPython.embed(colors='Linux') 
        # Appending to dictionaries 
        dict_pre_process[i]['time_correlation'][j] = auto_correlation 
        dict_pre_process[i]['length_scales'][j]    = length_scales
        dict_pre_process[i]['length_scales_k'][j]  = scales 
        dict_pre_process[i]['spe'][j] = spe  
        cutoff_scale.append(scales['cutoff_k']) 
    mean_cutoff_k[i] = np.mean(cutoff_scale) 

# Spatial series 
auto_correlation = { }
length_scales    = { }
spe              = { }
length_scales_k  = { }
variable = 'U-X'
for i in line_keys: 
    velocity_x       = line.working_data[i][variable] 
    spatial_location = line.location[i]  
    [time_rows, spatial_columns] = np.shape(velocity_x) 
    # Sub sample on time rows 
    time_vector   = np.arange(0, time_rows, time_sub_sampling) 
    temp_radius_z = np.linspace(spatial_location[2], 
                                   spatial_location[5],
                                   spatial_columns) 
    radius_z      = temp_radius_z - np.min(temp_radius_z)  
    fluctuation   = line.reynolds_decomposition(i, variable) 
    # Dictionaries for populating 
    auto_correlation[i] = { }
    length_scales[i]    = { } 
    spe[i]              = { } 
    length_scales_k[i]  = { } 
    for j in time_vector: 
        auto_correlation_temp = line.auto_correlation(radius_z, 
                        fluctuation[j], 15) 
        length_scales_temp = line.length_scales(
                auto_correlation_temp['radius_x'], 
                auto_correlation_temp['correlation'])  
        spe_temp      = line.base_filter_decay(velocity_x[j])  
        scales_k_temp = line.base_length_scales_k(fluctuation[j],
            auto_correlation_temp, length_scales_temp, spe_temp)  

        auto_correlation[i][j] = auto_correlation_temp 
        length_scales[i][j]    = length_scales_temp 
        spe[i][j]              = spe_temp 
        length_scales_k[i][j]  = scales_k_temp 
        
# PLOTING SPECTRUMS AND COrrelations 
spatial_vec = np.arange(0, spatial_columns, spatial_sub_sampling)
for i in line_keys:  
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2) 
    fig.set_size_inches(15,8) 
    fig.suptitle(f'{i}, {variable}') 
# Plotting spatial spatial  
    for j in time_vector:
     # Plotting Spatial Correlation 
        integral = "{:.3e}".format(length_scales[i][j]['integral']) 
        taylor = "{:.3e}".format(length_scales[i][j]['taylor'])
        ax1.plot(auto_correlation[i][j]['radius_x'], 
                auto_correlation[i][j]['correlation'], 
                linewidth='2', label=f'T={j}, L={integral}, $\lambda$={taylor}')
        ax1.set_ylabel('Correlation, spatial')
        ax1.grid('-.')
        ax1.legend() 
        ax1.set_xlim(left=0.0) 

     # Plotting Spatial Energy Spectrum
        integral_k = "{:.3e}".format(length_scales_k[i][j]['integral_k'])
        taylor_k = "{:.3e}".format(length_scales_k[i][j]['taylor_k'])
        ax2.loglog(spe[i][j], label=f'T={j}, $L_k$={integral_k}, $\lambda_k$={taylor_k}') 
        ax2.set_ylabel('Energy Spectrum, spatial')
        ax2.grid('-.')
        ax2.set_xlim(left=0.9) 
        ax2.legend() 
        #ax2.set_xlim(right=np.floor(spe['l0'][i][-1])) 

    # Plotting Time Series Energy Spectrum
    for k in range(len(spatial_vec)):  
            # Plotting Time Serires Correlation 
        integral = "{:.3e}".format(
                dict_pre_process[i]['length_scales'][k]['integral'])
        taylor = "{:.3e}".format(
                dict_pre_process[i]['length_scales'][k]['taylor'])
        ax3.plot(
            dict_pre_process[i]['time_correlation'][k]['radius_x'], 
            dict_pre_process[i]['time_correlation'][k]['correlation'],
            linewidth='2', label=f'X={spatial_vec[k]}, L={integral}, $\lambda$={taylor}')
        ax3.set_ylabel('Correlation, time series')
        ax3.set_xlabel('Radius [m]') 
        ax3.grid('-.')
        ax3.legend() 
        ax3.set_xlim(left=0.0) 

        # Plotting Time Series Energy Spectrum 
        integral_k = "{:.3e}".format(
                dict_pre_process[i]['length_scales_k'][k]['integral_k'])
        taylor_k = "{:.3e}".format(
                dict_pre_process[i]['length_scales_k'][k]['taylor_k'])
        ax4.loglog(dict_pre_process[i]['spe'][k], 
                label=f'X={spatial_vec[k]}, $L_k$={integral_k}, $\lambda_k$={taylor_k}') 
        ax4.grid('-.')
        ax4.set_ylabel('Energy Spectrum, time series')
        ax4.set_xlabel('K-vector [1/m]') 
        ax4.legend() 
        ax4.set_xlim(left=0.9) 
    plt.savefig(os.path.join(line_save, 'correlation', f'{i}.png')) 
    plt.close() 


boxcar_dict   = { } 
full_dict     = { }
filtered_dict = { }
sgs_dict      = { } 
box_radius    = { }
for i in line_keys: 
    boxcar_dict[i]   = { } 
    full_dict[i]     = { }
    filtered_dict[i] = { }
    sgs_dict[i]      = { } 
    box_radius[i]    = { }
    for j in variables: 
        variable = line.working_data[i][j] 
        radius   = line.x_axis(i, j) 
        cutoff_k = mean_cutoff_k[i]  
        boxcar_dict[i][j]   = { } 
        full_dict[i][j]     = { }
        filtered_dict[i][j] = { }
        sgs_dict[i][j]      = { } 
        box_radius[i][j]    = { }
        for k in testing_vec:
            boxcar = line.base_boxcar_filter(radius, variable[k],
                    cutoff_k)  
            boxcar_dict[i][j][k]   = boxcar  
            full_dict[i][j][k]     = boxcar['full_variable'] 
            filtered_dict[i][j][k] = boxcar['filtered_variable'] 
            sgs_dict[i][j][k]      = boxcar['sgs_variable'] 
            box_radius[i][j][k]    = boxcar['box_radius']

for i in line_keys:
    for j in variables: 
        fig = plt.figure() 
        plt.gcf().set_size_inches(16,7) 
        gs  = gridspec.GridSpec(3, 2, width_ratios=[3,1]) 
        ax_0 = plt.subplot(gs[0])
        ax_1 = plt.subplot(gs[1])
        ax_2 = plt.subplot(gs[2])
        ax_3 = plt.subplot(gs[3])
        ax_4 = plt.subplot(gs[4])
        ax_5 = plt.subplot(gs[5])
        location    = line.location[i] 
        x_1         = [location[0], location[3]]
        x_2         = [location[1], location[4]]
        x_3         = [location[2], location[5]]
        window_size = boxcar_dict[i][j][0]['window_size']
        pos_str     = f'x_1={x_1}, x_2={x_2}, x_3={x_3}' 
        fig.suptitle(f'{i}, {j} at {pos_str}, window_size = {window_size}') 
        for k in testing_vec:
            full_hist     = full_dict[i][j][k]
            filtered_hist = filtered_dict[i][j][k]
            sgs_hist      = sgs_dict[i][j][k] 
            radius        = box_radius[i][j][k]
        # Ax 0 
            ax_0.plot(radius, full_hist, linewidth='0.7', 
                      label=f'Time = {k}')
            ax_0.set_ylabel(f'{j}, full')
            ax_0.grid('-.') 
            ax_0.margins(x=0) 
            ax_0.legend() 
        # Ax 1
            ax_1.hist(full_hist, bins=50, density=True, 
                    color='lightcyan', edgecolor='k')
            sb.kdeplot(full_hist, bw_adjust=1.3, fill=False, 
                    linewidth=0.7,color='crimson', ax=ax_1) 

        # Ax 2
            ax_2.plot(radius, filtered_hist, linewidth='0.7', 
                      label=f'Time = {k}')
            ax_2.set_ylabel(f'{j}, filtered')
            ax_2.grid('-.') 
            ax_2.margins(x=0) 
            ax_2.legend() 
        # Ax 3
            ax_3.hist(filtered_hist, bins=50, density=True, 
                    color='lightcyan', edgecolor='k')
            sb.kdeplot(filtered_hist, bw_adjust=1.3, fill=False, 
                    linewidth=0.7, color='crimson', ax=ax_3) 
        # Ax 4
            ax_4.plot(radius, sgs_hist, linewidth='0.7',
                      label=f'Time = {k}')
            ax_4.set_ylabel(f'{j}, sgs')
            ax_4.grid('-.') 
            ax_4.margins(x=0) 
            ax_4.legend() 
            ax_4.set_xlabel('Radius [m]') 
        # Ax 5
            ax_5.hist(sgs_hist, bins=50, density=True, 
                    color='lightcyan', edgecolor='k')
            sb.kdeplot(sgs_hist, bw_adjust=1.3, fill=False, 
                    linewidth=0.7, color='crimson', ax=ax_5) 
        plt.savefig(os.path.join(line_save,'histogram', f'{i}_{j}.png')) 
        plt.close() 

# Add new data to structures if this is on 
if (new_data_flag == True):
# Add new data probe  
    for i in probe_keys: 
        pressure     = probe.working_data[i]['P'] 
        dilatation   = probe.working_data[i]['DIL'] 
        rho          = probe.working_data[i]['RHO'] 
        pressure_dil = pressure * dilatation 
        rho_dil      = rho * dilatation 
        probe.add_new_variable(i, 'RHO-DIL', rho_dil, 
                pickle_path, 'new_probe_data') 
# Add new data line 
    for i in line_keys: 
        pressure     = line.working_data[i]['P'] 
        dilatation   = line.working_data[i]['DIL'] 
        rho          = line.working_data[i]['RHO'] 
        pressure_dil = pressure * dilatation 
        rho_dil      = rho * dilatation 
        line.add_new_variable(i,'RHO-DIL', rho_dil, 
                pickle_path, 'new_line_data') 

if (probe_flag == True):
    for i in probe_keys:
        print(i) 
    # Plot Scatter 
        probe_radius = probe.x_axis(i, sampling_flag=sub_sampling_flag) 
        # Constant length scale (U-X)
        const_cutoff_k  = probe.const_cutoff_k(i, probe_radius, 
                correlation_lag=probe_correlation_lag)  

        if (scatter_probe_flag == True):
            probe.plot_scatter(i, 'P', 'RHO', probe_radius, const_cutoff_k,  
                               correlation_lag=probe_correlation_lag,
                               saving_path=scatter_save) 

            probe.plot_scatter(i, 'P', 'T', probe_radius, const_cutoff_k,  
                               correlation_lag=probe_correlation_lag,
                               saving_path=scatter_save) 

            probe.plot_scatter(i, 'RHO', 'T', probe_radius, const_cutoff_k, 
                               correlation_lag=probe_correlation_lag,
                               saving_path=scatter_save) 

            probe.plot_scatter(i, 'P', 'DIL', probe_radius, const_cutoff_k, 
                               correlation_lag=probe_correlation_lag,
                               saving_path=scatter_save) 

            probe.plot_scatter(i, 'T', 'DIL', probe_radius, const_cutoff_k, 
                               correlation_lag=probe_correlation_lag,
                               saving_path=scatter_save) 

            probe.plot_scatter(i, 'RHO', 'DIL', probe_radius, const_cutoff_k,  
                               correlation_lag=probe_correlation_lag,
                               saving_path=scatter_save) 

            probe.plot_scatter(i, 'RHO', 'P-DIL', probe_radius, const_cutoff_k,  
                               correlation_lag=probe_correlation_lag,
                               saving_path=scatter_save) 

            probe.plot_scatter(i, 'T', 'P-DIL', probe_radius, const_cutoff_k,  
                               correlation_lag=probe_correlation_lag,
                               saving_path=scatter_save) 

        for j in variables:  
            print(i, j)  
        # General Functions 
            variable             = probe.working_data[i][j] 
            variable             = probe.sub_sampling(variable) 
            variable_fluctuation = probe.reynolds_decomposition(variable) 
            auto_correlation     = probe.auto_correlation(probe_radius,
                                    variable_fluctuation, 
                                    auto_correlation_len=probe_correlation_lag) 
            pwr_spectral         = probe.filter_decay(variable) 
            length_scales        = probe.length_scales(
                             auto_correlation['correlation_radius'], 
                                    auto_correlation['correlation'],  
                                    variable_fluctuation, pwr_spectral) 
            boxcar_dict          = probe.boxcar_filter(probe_radius, 
                                   variable, const_cutoff_k) 
            variable_moments     = probe.raw_stat_moments(variable)
            legendre_dict        = probe.legendre_interpolation(boxcar_dict) 
            moments_str          = probe.statistical_moments_str(boxcar_dict) 

            # Shifting Factors 
            if (j == 'RHO'):
                shifting_factor = 1.6 
            elif (j == 'P'):
                shifting_factor = 1e9 
            elif (j == 'DIL') :
                shifting_factor = 1e12 
            elif (j == 'RHOE'):
                shifting_factor = 1e10 
            elif (j == 'VORTMAG'):
                shifting_factor = 1e12 
            elif (j == 'P-DIL'):
                shifting_factor = 1e21 
            elif (j == 'RHO-DIL'):
                shifting_factor = 1e10
            else: 
                shifting_factor = 1e7 

            # Plots  
            probe.plot_boxcar(i,j, boxcar_dict, moments_str, 
                        saving_path=boxcar_save) 
            probe.plot_legendre(i, j, boxcar_dict, legendre_dict,
                                saving_path=legendre_save) 
            probe.plot_results(i,j, probe_radius, variable, pwr_spectral,
                    auto_correlation, length_scales, moments_str,
                  inertial_shifting_factor=shifting_factor,
                  saving_path=probe_save) 





if (line_flag == True):
    for i in line_keys:
        for j in variables:
            print(i, j)  
            if (j == 'RHO'):
                shifting_factor = 1.5 
            elif (j == 'P'):
                shifting_factor = 1e9 
            else: 
                shifting_factor = 1e6 

            # Plot all results 
            line.plot_results(i,j, 
                sampling_rate=line_sampling_rate,
                correlation_len=line_correlation_lag, 
                inertial_shifting_factor=shifting_factor, 
                saving_path=line_save) 
