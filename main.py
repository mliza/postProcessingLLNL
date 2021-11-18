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
line_flag             = True 
probe_flag            = True  
new_data_flag         = True  
scatter_probe_flag    = True
sub_sampling_flag     = False
probe_sampling_rate   = 1 
probe_correlation_lag = 50 
line_correlation_lag  = 50 
time_sub_sampling     = 1
spatial_sub_sampling  = 1 

# Paths 
pickle_path    = '/Users/martin/Documents/Research/UoA/Projects/LLNL/data/data_5/pickle' 
save_path      = '/Users/martin/Desktop/workingFiles/results' 
probe_save     = os.path.join(save_path, 'probe', 'total')
probe_boxcar   = os.path.join(save_path, 'probe', 'boxcar')
probe_scatter  = os.path.join(save_path, 'probe', 'scatter')
probe_legendre = os.path.join(save_path,'probe', 'legendre') 
line_save      = os.path.join(save_path, 'line')
line_correlation_save = os.path.join(line_save, 'correlation') 
line_boxcar    = os.path.join(save_path, 'line', 'boxcar')
line_legendre  = os.path.join(save_path,'line', 'legendre') 
temporal_line_legendre = os.path.join(line_legendre, 'temporal') 
spatial_line_legendre  = os.path.join(line_legendre, 'spatial') 
temporal_line_boxcar   = os.path.join(line_boxcar, 'temporal') 
spatial_line_boxcar    = os.path.join(line_boxcar, 'spatial') 

# Loading class's instances  
base  = baseAnalysis.Base_Analysis(pickle_path)
probe = probe.Probe(pickle_path, sampling_rate=probe_sampling_rate) 
line  = line.Line(pickle_path, sampling_rate=probe_sampling_rate) 

# Keys 
probe_keys = list(probe.location.keys()) 
line_keys  = list(line.location.keys()) 
variables  = ['U-X', 'U-Y', 'U-Z', 'P', 'T', 'RHO', 'RHOE', 'GRADRHOMAG', 'DIL', 'VORTMAG', 'P-DIL', 'RHO-DIL'] 
variables  = ['U-X', 'U-Z', 'U-Y', 'P', 'T', 'RHO', 'DIL', 'P_DIL']  

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
        radius_dict = probe.x_axis(i, sampling_flag=sub_sampling_flag) 
        probe_radius = radius_dict['x_axis'] 
        # Constant length scale (U-X)
        const_cutoff_k  = probe.const_cutoff_k(i, probe_radius, 
                correlation_lag=probe_correlation_lag)  

        if (scatter_probe_flag == True):
            probe.plot_scatter(i, 'P', 'RHO', probe_radius, const_cutoff_k,  
                               correlation_lag=probe_correlation_lag,
                               saving_path=probe_scatter) 

            probe.plot_scatter(i, 'P', 'T', probe_radius, const_cutoff_k,  
                               correlation_lag=probe_correlation_lag,
                               saving_path=probe_scatter) 

            probe.plot_scatter(i, 'RHO', 'T', probe_radius, const_cutoff_k, 
                               correlation_lag=probe_correlation_lag,
                               saving_path=probe_scatter) 

            probe.plot_scatter(i, 'P', 'DIL', probe_radius, const_cutoff_k, 
                               correlation_lag=probe_correlation_lag,
                               saving_path=probe_scatter) 

            probe.plot_scatter(i, 'T', 'DIL', probe_radius, const_cutoff_k, 
                               correlation_lag=probe_correlation_lag,
                               saving_path=probe_scatter) 

            probe.plot_scatter(i, 'RHO', 'DIL', probe_radius, const_cutoff_k,  
                               correlation_lag=probe_correlation_lag,
                               saving_path=probe_scatter) 

            probe.plot_scatter(i, 'RHO', 'P-DIL', probe_radius, const_cutoff_k,  
                               correlation_lag=probe_correlation_lag,
                               saving_path=probe_scatter) 

            probe.plot_scatter(i, 'T', 'P-DIL', probe_radius, const_cutoff_k,  
                               correlation_lag=probe_correlation_lag,
                               saving_path=probe_scatter) 

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
                        saving_path=probe_boxcar) 
            probe.plot_legendre(i, j, boxcar_dict, legendre_dict,
                                saving_path=probe_legendre) 
            probe.plot_results(i,j, probe_radius, variable, pwr_spectral,
                    auto_correlation, length_scales, moments_str,
                  inertial_shifting_factor=shifting_factor,
                  saving_path=probe_save) 


# Line and crunched data  
temporal_dict = { }
spatial_dict  = { }
if (line_flag == True):
    for i in line_keys: 
        temporal_dict[i] = { }
        spatial_dict[i]  = { }  
        for j in variables:
            print(i, j)  
            temp_dict = line.temporal_data(i,j, n_points=spatial_sub_sampling, 
                                    auto_correlation_len=line_correlation_lag) 
            spat_dict = line.spatial_data(i,j, n_points=time_sub_sampling, 
                                    auto_correlation_len=line_correlation_lag) 
            temporal_dict[i][j] = line.data_cruncher(temp_dict) 
            spatial_dict[i][j]  = line.data_cruncher(spat_dict) 
    # Calculate Scales 
            line.plot_correlation_spe(temporal_dict[i][j], 
                               spatial_dict[i][j], dataset=i, variable=j,
                               time_sub_sampling=time_sub_sampling,
                               spatial_sub_sampling=spatial_sub_sampling,
                               saving_path = line_correlation_save) 
    # Testing scales 
        for j in variables:
            # Spatial 
            spatial_cutoff_k = spatial_dict[i]['U-Z']['length_scales']['cutoff_k']
            spatial_boxcar   = line.boxcar_filter(spatial_dict[i][j]['radius'], 
                                              spatial_dict[i][j]['variable'], 
                                              spatial_cutoff_k)
            spatial_moments_str = line.statistical_moments_str(spatial_boxcar) 
            spatial_legendre    = line.legendre_interpolation(spatial_boxcar) 
            spatial_dict[i][j]['boxcar']     = spatial_boxcar
            spatial_dict[i][j]['legendre']    = spatial_legendre
            spatial_dict[i][j]['moments_str'] = spatial_moments_str
            # Plots 
            line.plot_boxcar(i,j, spatial_boxcar, spatial_moments_str, 
                        saving_path=spatial_line_boxcar) 
            line.plot_legendre(i, j, spatial_boxcar, spatial_legendre,
                                saving_path=spatial_line_legendre) 

            # Temporal 
            temporal_cutoff_k = temporal_dict[i]['U-Z']['length_scales']['cutoff_k']
            temporal_boxcar   = line.boxcar_filter(temporal_dict[i][j]['radius'], 
                                              temporal_dict[i][j]['variable'], 
                                              temporal_cutoff_k)
            temporal_moments_str = line.statistical_moments_str(temporal_boxcar) 
            temporal_legendre    = line.legendre_interpolation(temporal_boxcar) 
            temporal_dict[i][j]['boxcar']      = temporal_boxcar
            temporal_dict[i][j]['legendre']    = temporal_legendre
            temporal_dict[i][j]['moments_str'] = temporal_moments_str

            # Plots 
            line.plot_boxcar(i,j, temporal_boxcar, temporal_moments_str, 
                        saving_path=temporal_line_boxcar) 
            line.plot_legendre(i, j, temporal_boxcar, temporal_legendre,
                                saving_path=temporal_line_legendre) 