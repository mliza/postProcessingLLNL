#!/opt/homebrew/bin/python3.9
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
import line 

# Configuration Parameters 
line_flag             = False    
probe_flag            = True  
new_data_flag         = False  
scatter_probe_flag    = True  
probe_sampling_rate   = 1 
sub_sampling_flag     = False
probe_correlation_lag = 100 
line_correlation_lag  = 40 
time_sub_sampling     = 500
spatial_sub_sampling  = 100 

# Paths 
pickle_path   = '/Users/martin/Documents/Research/UoA/Projects/LLNL/data/data_5/pickle' 
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

# Line and crunched data  
temporal_dict = { }
spatial_dict  = { }
for i in line_keys: 
    temp_dict = line.temporal_data(i, n_points=spatial_sub_sampling, 
                                auto_correlation_len=line_correlation_lag) 
    spat_dict = line.spatial_data(i, n_points=time_sub_sampling, 
                                auto_correlation_len=line_correlation_lag) 
    temporal_dict[i] = line.data_cruncher(temp_dict) 
    spatial_dict[i]  = line.data_cruncher(spat_dict) 
# Calculate Scales 
    line.plot_correlation_spe( temporal_dict[i], 
                           spatial_dict[i], dataset=i, 
                           variable='U-X',
                           time_sub_sampling=time_sub_sampling,
                           spatial_sub_sampling=spatial_sub_sampling,
                           saving_path = line_save) 
IPython.embed(colors='Linux') 
        

'''
return_dict = { }
for i in line_keys: 
    radius_z      = line.x_axis(i) 
    time_cutoff_k = float(time_scales[i]['cutoff_k'])
    return_dict[i] = { }
    for j in variables: 
        return_dict[i][j] = { }
        variable = line.working_data[i][j]  
        solution_dict = { } 
        [time_len, spatial_len] = np.shape(variable) 
        solution_dict['full_radius']       = [ ]
        solution_dict['full_variable']     = [ ]
        solution_dict['filtered_variable'] = [ ]
        solution_dict['sgs_variable']      = [ ]
        solution_dict['full_moments']      = [ ]
        solution_dict['filtered_moments']  = [ ]
        solution_dict['sgs_moments']       = [ ]
        for k in chosen_times:  
            boxcar_dict = line.boxcar_filter(radius_z, variable[k], time_cutoff_k) 
            temp_radius = boxcar_dict['full_radius'] 
            temp_full   = boxcar_dict['full_variable']
            temp_filt   = boxcar_dict['filtered_variable']
            temp_sgs    = boxcar_dict['sgs_variable'] 
            solution_dict['full_radius'].append(temp_radius) 
            solution_dict['full_variable'].append(temp_full) 
            solution_dict['filtered_variable'].append(temp_filt) 
            solution_dict['sgs_variable'].append(temp_sgs_m)

        IPython.embed(colors='Linux') 


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
'''
# Add new data to structures if this is on 
if (new_data_flag is True):
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

if (probe_flag is True):
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
