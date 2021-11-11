#!/usr/local/bin/python3.9
'''
    Date:   07/09/2021
    Author: Martin E. Liza
    File:   baseAnalysis.py
    Def:               

    Author           Date         Revision
    --------------------------------------------------------------- 
    Martin E. Liza   06/17/2021   Initial Version.
    Martin E. Liza   06/21/2021   Added Line and Probe classes. 
    Martin E. Liza   06/24/2021   Cleaned up code (coding standards)
                                  and split the child classes in two
                                  different files. 
    Martin E. Liza   07/09/2021   Added raw_stat_moments, function.
    Martin E. Liza   07/11/2021   Added sub_sampling, function.
    Martin E. Liza   07/14/2021   Fixed base_correlation.
    Martin E. Liza   07/15/2021   Fixed base_correlation, rename sub_sampling to 
                                  base_sub_sampling; and clean up code. 
    Martin E. Liza   07/19/2021   Fixed reynolds_decomposition to take into 
                                  consideration 2d, dataset such as line data.
    Martin E. Liza   07/21/2021   Implement vorticity_calculations. 
    Martin E. Liza   08/29/2021   Added add_new_variable. 
'''
import os 
import sys 
import pickle
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec
from scipy import integrate
from scipy.fft import fft 
from scipy.special import legendre 
import IPython

class Base_Analysis:
    # Global class variables 
    flag_type   = None 
    flag_points = None

    def __init__(self, pickle_path, sampling_rate=1): 
        data               = self.load_data(pickle_path) 
        self.sampling_rate = sampling_rate
        # If None returns all data 
        if (self.flag_type == None):
            self.working_data = data
        else: 
            self.working_data = data[self.flag_type] 
            self.location     = data[self.flag_points] 

    # Loads all the .pickle files from a folder and storage 
    # them in  multi dimensions dictionary 
    def load_data(self, pickle_path): 
        pickle_files = os.listdir(pickle_path) 
        data         = { }
        # Iterates through all pickle files in pickle path 
        for i in range(len(pickle_files)): 
            data_pickle       = os.path.join(pickle_path, pickle_files[i]) 
            pickle_name       = pickle_files[i].replace('.pickle','') 
            data[pickle_name] = pickle.load(open(data_pickle, 'rb')) 
        return data 

# Adds a new variable 
    def add_new_variable(self, dataset_number, new_variable_name, 
                new_variable_data, saving_path, saving_name):
        data_total                    = self.working_data[dataset_number]
        data_total[new_variable_name] = new_variable_data 
        tot_saving_path = os.path.join(saving_path, f'{saving_name}.pickle')
        pickle_out      = open(tot_saving_path, 'wb') 
        pickle.dump(self.working_data, pickle_out) 
        pickle_out.close() 

# Calculates Reynolds decomposition
    def reynolds_decomposition(self, variable):
        fluctuation = variable - np.mean(variable)  
        return fluctuation 

# Calculates the auto correlation 
    def auto_correlation(self, x_axis, fluctuation, auto_correlation_len=50):  
        fluctuation_len     = len(fluctuation) - 1  # there should be 1 less element
        numerator           = np.zeros(auto_correlation_len) 
        denominator         = np.zeros(auto_correlation_len)  
        # Re shape x_axis to be the same length as the correlation 
        max_x_len           = np.max(x_axis)  
        auto_correlation_x  = np.linspace(0, max_x_len, auto_correlation_len) 
        for i in range(fluctuation_len): 
            k = i 
            for j in range(auto_correlation_len): 
                numerator[j]   += (fluctuation[i] * fluctuation[k]) 
                denominator[j] += fluctuation[i]**2  
                k += 1
                if (k > fluctuation_len):
                    break 
        auto_correlation = numerator / denominator
        autoCorr_dict = { 'correlation_radius' : auto_correlation_x, 
                          'correlation'        : auto_correlation }
        return autoCorr_dict  

# Calculate length scales 
    def length_scales(self, corr_radius, correlation, fluctuation, spe): 
        delta_x = corr_radius[1] - corr_radius[0] 
    # Calculates the Taylor microscale
        nd_derivative  = (2*correlation[0] - 5*correlation[1]  
                            + 4*correlation[2] - correlation[3]) / delta_x**2
        taylor_scale   = 1 / np.sqrt(-0.5 * nd_derivative) 
    # Calculates the Integral Scales
        # sometimes the length is negative because of the data 
        integral_scale = np.abs(integrate.simpson(correlation, dx=delta_x)) 

    # Converts them to k_space 
        dist_k = range(len(spe))
    # Mapping from a small data set to a big data set <?> 
        ord_mag       = np.floor(np.log10(dist_k[-1]))
        integral_norm = integral_scale / corr_radius[-1] 
        taylor_norm   = taylor_scale / corr_radius[-1] 
    # Convert them to k-space 
        dist_ratio = len(fluctuation) / dist_k[-1]  
        integral_k = (1 / integral_norm) * ord_mag * dist_ratio 
        taylor_k   = (1 / taylor_norm) * ord_mag * dist_ratio  
        k_delta    = (taylor_k + integral_k) / 2

    # Returns a dictionary 
        scales_dict = { 'taylor_m'   : "{:.3e}".format(taylor_scale),
                        'integral_m' : "{:.3e}".format(integral_scale), 
                        'integral_k' : "{:.3e}".format(integral_k), 
                        'cutoff_k'   :  k_delta,  
                        'taylor_k'   : "{:.3e}".format(taylor_k) }
        return scales_dict 

# Calculates filter decay 
    def filter_decay(self, variable, n_bin_number=2): 
        # Filter setup 1/10 decay
        data_size            = np.shape(variable)[0] 
        variable_hat         = fft(variable, data_size) 
        pwr_spectral_density = ( variable_hat * np.conj(variable_hat) 
                                 / data_size ).real   
        bin_number           = int( np.floor(data_size / n_bin_number) + 1 ) 
        bin_vector           = range(1, bin_number + 1)
        pwf                  = np.zeros(bin_number) 
        log_freq             = np.log10(bin_vector) 

        # Butterworth filter spectrum with 1/10th decade filter 
        for i in range(bin_number - 1):
            count   = 0
            sum_pwr_spectral_density = 0
            log_freq_0 = log_freq[i] - 0.05 #initial 
            log_freq_t = log_freq[i] + 0.05 #final 
            for j in range(bin_number - 1):
                if (log_freq[j] > log_freq_0 and log_freq[j] < log_freq_t):
                    count += 1
                    sum_pwr_spectral_density += pwr_spectral_density[j] 
            if (count != 0):
                pwf[i] = sum_pwr_spectral_density / count 
        # Append a 0 at the beginning of the pwr spectral 
        pwr_spectral = np.append(0, pwf[0:bin_number-2]) 
        return pwr_spectral  

# Sub sampling 
    def sub_sampling(self, variable, sampling_in=[]): 
        if (sampling_in == []): 
            sampling = self.sampling_rate 
        else: 
            sampling = sampling_in 

        if (sampling == 1):
            sampled_var = variable
        else:
            if (variable.ndim == 1):
                total_len      = len(variable) # initial length  
                max_len        = int(np.ceil(len(variable) / sampling))
                sampling_index = range(0, max_len * sampling, sampling) 
                sampled_var    = variable[sampling_index] 
          
            if (variable.ndim == 2):
                [rows, columns]    = np.shape(variable) 
                row_sampling_rate  = sampling[0] 
                col_sampling_rate  = sampling[1] 
                max_row_len        = int(np.ceil(rows / row_sampling_rate)) 
                max_col_len        = int(np.ceil(columns / 
                                     col_sampling_rate)) 
                row_sampling_index = range(0, max_row_len*row_sampling_rate,
                                     row_sampling_rate)     
                col_sampling_index = range(0, max_col_len*col_sampling_rate,
                                     col_sampling_rate)     
                sampled_var = [ ]
                for i in row_sampling_index: 
                    sampled_var.append(variable[i][col_sampling_index]) 
        return sampled_var 

# Calculates the statistic moments 
    def raw_stat_moments(self, variable):
        var2     = variable**2 
        var3     = variable**3 
        var4     = variable**4 
        mu1      = np.mean(variable) 
        mu2      = np.mean(var2)  
        mu3      = np.mean(var3)  
        mu4      = np.mean(var4)  
        
        # Calculates the statistics moments 
        mean     = np.mean(variable) 
        variance = mu2 - mu1**2 
        skewness = (mu3 - 3 * mu2 * mu1 + 2 * mu1**3) / variance**(3/2) 
        kurtosis = ( (mu4 - 4 * mu3 * mu1 + 6 * mu2 * mu1**2 - 
                    3* mu1**4) / variance**2 )
        # Returns a dictionary  
        raw_moments = { 'mean'    : "{:.3e}".format(mean), 
                        'std_dev' : "{:.3e}".format(np.sqrt(variance)),
                        'skewness': "{:.3e}".format(skewness), 
                        'kurtosis': "{:.3e}".format(kurtosis) }
        return raw_moments  

# Boxcar filter (Moving average filter, Finite Impulse Response filter)  
    def boxcar_filter(self, radius, variable, cutoff_scale): 
     # Loading variables 
        window_size  = int(np.floor(cutoff_scale)) 
        # Divide the window size by 2 to pick the value at the center
        delta        = int(np.floor((window_size)/2))
        # Moving Avg. filter, boxcar filter 
        filtered_variable = [ ]
        for i in range(delta,len(variable)-delta):
            temp_x = [ ]
            j = i - delta
            for k in range(window_size):
                temp_x.append(variable[j+k]) 
            filtered_variable.append(np.sum(temp_x)/window_size) 
    # Calculating SGS and cleaning variable in and radius 
        filtered_variable = np.array(filtered_variable)
        # Box variable is raw data within the box 
        box_variable      = variable[delta:len(variable) - delta]
        box_radius        = radius[delta:len(variable) - delta]
        sgs_variable      = filtered_variable - box_variable 
        les_variable      = self.sub_sampling(filtered_variable,
                                               window_size) 
        les_radius        = self.sub_sampling(box_radius, window_size)
        filter_dict   = { 'full_radius'       : box_radius,  
                          'full_variable'     : box_variable,  
                          'filtered_variable' : filtered_variable, 
                          'sgs_variable'      : sgs_variable,  
                          'les_radius'        : les_radius, 
                          'les_variable'      : les_variable,
                          'window_size'       : window_size }
        return filter_dict 

# Legendre Interpolation 
    def legendre_interpolation(self, filter_dict): 
        les_variable = filter_dict['les_variable']
    # Legendre documentation 
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.legendre.html 
    # Define the Legendre polynomials 
        legendre_pol_0  = legendre(0)
        legendre_pol_1  = legendre(1)
        legendre_pol_2  = legendre(2)
    # Create a (3x3) Legendre matrix  
        x_array = np.linspace(-1,1,3)
        legendre_matrix = [ ]
        for i in x_array: 
            pol_0 = legendre_pol_0(i)
            pol_1 = legendre_pol_1(i)
            pol_2 = legendre_pol_2(i)
            legendre_matrix.append([pol_0, pol_1, pol_2]) 
        legendre_matrix = np.matrix(legendre_matrix).transpose()
    # Filter  
        sol_matrix = [ ]
        for i in range(len(les_variable)-3):
            les_variable_temp = les_variable[i:i+3] 
            sol_temp = np.linalg.solve(legendre_matrix,
                                      les_variable_temp) 
            sol_matrix.append(sol_temp) 
    # Calculating delta variables 
        constants_matrix = np.asarray(sol_matrix) 
        a_0 = np.squeeze(constants_matrix[:,0])
        variable_delta  = les_variable[:-3] - a_0  
        variable_2delta = a_0
        radius_delta    = filter_dict['les_radius'][:-3]  
    # Return a dictionary 
        filter_dict = { 'delta_variable'  : variable_delta,  
                        'delta2_variable' : variable_2delta,  
                        'delta_radius'    : radius_delta }
        return filter_dict 

# Moments String 
    def statistical_moments_str(self, boxcar_dict):  
        full_moments     = self.raw_stat_moments(boxcar_dict['full_variable']) 
        resolved_moments = self.raw_stat_moments(boxcar_dict['filtered_variable']) 
        sgs_moments      = self.raw_stat_moments(boxcar_dict['sgs_variable']) 
        # Creating strings 
        moment_labels    = list(full_moments.keys()) 
        full_mom_val     = list(full_moments.values())
        resolved_mom_val = list(resolved_moments.values())
        sgs_mom_val      = list(sgs_moments.values())

        # Create moment strings 
        nl = '\n'
        full_str    = f'{moment_labels[0]} = {full_mom_val[0]}{nl}{moment_labels[1]} = {full_mom_val[1]}{nl}{moment_labels[2]} = {full_mom_val[2]}{nl}{moment_labels[3]} = {full_mom_val[3]}' 
        resolved_str     = f'{moment_labels[0]} = {resolved_mom_val[0]}{nl}{moment_labels[1]} = {resolved_mom_val[1]}{nl}{moment_labels[2]} = {resolved_mom_val[2]}{nl}{moment_labels[3]} = {resolved_mom_val[3]}' 
        sgs_str     = f'{moment_labels[0]} = {sgs_mom_val[0]}{nl}{moment_labels[1]} = {sgs_mom_val[1]}{nl}{moment_labels[2]} = {sgs_mom_val[2]}{nl}{moment_labels[3]} = {sgs_mom_val[3]}' 
        moment_str  =  { 'full_str'      : full_str,  
                         'resolved_str'  : resolved_str, 
                         'sgs_str'       : sgs_str }
        return moment_str

