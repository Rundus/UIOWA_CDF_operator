# ------------------------------------------------------------------------------------------------------
# A Variables File that contains everything relevant in the parent eepaa_low and eepaa_high counts files
# ------------------------------------------------------------------------------------------------------
import numpy as np
from files import EEPAA_file_high,EEPAA_file_low,Flux_aligned_high,Flux_aligned_low,dist_file_high,dist_file_low,J_par_file_low
from files import counts_file_high,counts_file_low,Mag_file_low,Mag_file_high,pitch_file_actual_high,pitch_file_actual_low
from files import J_par_file_high,root
from files import cal_file_counts_high,cal_file_counts_low
from files import cal_diff_file_high,cal_diff_file_low
from files import sun_spike_noise_file_high,sun_spike_noise_file_low
from files import spike_rem_counts_high_file,spike_rem_counts_low_file,spike_rem_file_statistics_low,spike_rem_file_statistics_high
# --------------
# RAW COUNT DATA
# --------------

# The counts data has 1 supportingDocs slice in time and needs to be removed to match EPOCH
counts_high_raw = counts_file_high.varget('eepaa')
counts_low_raw = counts_file_low.varget('eepaa')
EPOCH_high_counts = counts_file_high.varget('Epoch')
EPOCH_low_counts = counts_file_low.varget('Epoch')
counts_high_info = counts_file_high.cdf_info()
counts_low_info = counts_file_low.cdf_info()
zvars_counts_high = counts_high_info['zVariables']
zvars_counts_low = counts_low_info['zVariables']

counts_high = np.zeros(shape=(len(counts_high_raw)-1,len(counts_high_raw[0]),len(counts_high_raw[0][0])), dtype = 'float64')
counts_low = np.zeros(shape=(len(counts_low_raw)-1,len(counts_low_raw[0]),len(counts_low_raw[0][0])), dtype = 'float64')

for tme in range(1,len(counts_high)):
    counts_high[tme] = counts_high_raw[tme]

for tme in range(1,len(counts_low)):
    counts_low[tme] = counts_low_raw[tme]


# ---------------------
# CALIBRATED COUNT DATA
# ---------------------

#using calibrated counts
# cal_counts_high = cal_file_counts_high.varget('Calibrated_Counts')
# cal_counts_low = cal_file_counts_low.varget('Calibrated_Counts')

#Using Rem-spike removed counts
cal_counts_high = cal_file_counts_high.varget('removed_spike_counts_threshed_data')
cal_counts_low = cal_file_counts_low.varget('removed_spike_counts_threshed_data')
cal_counts_high_info = cal_file_counts_high.cdf_info()
cal_counts_low_info = cal_file_counts_low.cdf_info()
zvars_cal_counts_high = cal_counts_high_info['zVariables']
zvars_cal_counts_low = cal_counts_low_info['zVariables']


# ---------------------------------
# CALIBRATED DIFFERENTIAL FLUX DATA
# ---------------------------------
diffEFlux_high = cal_diff_file_high.varget(variable="Differential_Energy_Flux")
diffNFlux_high = cal_diff_file_high.varget(variable="Differential_Number_Flux")
diffEFlux_low = cal_diff_file_low.varget(variable="Differential_Energy_Flux")
diffNFlux_low = cal_diff_file_low.varget(variable="Differential_Number_Flux")
cal_diff_high_info = cal_diff_file_high.cdf_info()
cal_diff_low_info = cal_diff_file_low.cdf_info()
zvars_cal_diff_high = cal_diff_high_info['zVariables']
zvars_cal_diff_low = cal_diff_low_info['zVariables']

# -----------------------
# INITIAL VARIABLES TRICE
# -----------------------

# get some variables
geometric_high = EEPAA_file_high.varget('geometric_factor')
deadtime_high = EEPAA_file_high.varget('dead_time')
geometric_low = EEPAA_file_low.varget('geometric_factor')
deadtime_low = EEPAA_file_low.varget('dead_time')
pitch = EEPAA_file_low.varget('Pitch_Angle')
Energies = EEPAA_file_high.varget('Energy')
count_interval_high = EEPAA_file_high.varget('count_interval')
count_interval_low = EEPAA_file_low.varget('count_interval')



#Get the input file variable information
EEPAA_file_high_info = EEPAA_file_high.cdf_info()
EEPAA_file_low_info = EEPAA_file_low.cdf_info()

# Get the input file zVariables
zvars_high = EEPAA_file_high_info['zVariables']
zvars_low = EEPAA_file_low_info['zVariables']

#DiffFlux (UNCALIBRATED)
DiffEFlux_high = EEPAA_file_high.varget(variable="Differential_Energy_Flux")
DiffNFlux_high = EEPAA_file_high.varget(variable="Differential_Number_Flux")
DiffEFlux_low = EEPAA_file_low.varget(variable="Differential_Energy_Flux")
DiffNFlux_low = EEPAA_file_low.varget(variable="Differential_Number_Flux")

#Get the fillvalues from the Diffflux variables
atts_high = EEPAA_file_high.varattsget(zvars_high[13])
atts_low = EEPAA_file_low.varattsget(zvars_low[13])
Fillhigh = atts_high["FILLVAL"]
Filllow =  atts_low["FILLVAL"]
Energy_High = EEPAA_file_high.varget('Energy')
Energy_low = EEPAA_file_low.varget('Energy')





# EPOCHS
EPOCH_High = EEPAA_file_high.varget('Epoch')
EPOCH_low = EEPAA_file_low.varget('Epoch')
#
# print('------------ Energy & PITCH -------------')
# print(counts_high_info)
# print(counts_file_high.varattsget(zvars_counts_high[16]))
# print(counts_file_high.varinq(zvars_counts_high[16]))
# print('\n')
# print(counts_file_high.varattsget(zvars_counts_high[9]))
# print(counts_file_high.varinq(zvars_counts_high[9]))
#
#
# print('\n')
# print('------------ flux -------------')
# print(atts_high)
# print(EEPAA_file_high.varinq(zvars_high[13]))
# print('\n')
# print(atts_low)
# print(EEPAA_file_low.varinq(zvars_low[13]))
# print('\n')
# print('------------ epoch -------------')
# print(EEPAA_file_high.varattsget(zvars_high[0]))
# print(EEPAA_file_high.varinq(zvars_high[0]))
# print('\n')
# print(EEPAA_file_high.varattsget(zvars_high[0]))
# print(EEPAA_file_high.varinq(zvars_high[0]))



#Pitch Angles)
pitch_angle_actual_high_info = pitch_file_actual_high.cdf_info()
pitch_angle_actual_low_info = pitch_file_actual_low.cdf_info()
zvars_pitch_high = pitch_angle_actual_high_info['zVariables']
zvars_pitch_low = pitch_angle_actual_low_info['zVariables']
pitch_output_high = pitch_file_actual_high.varget(variable= 'Pitch_Angle_Actual')
pitch_output_low = pitch_file_actual_low.varget(variable= 'Pitch_Angle_Actual')

# Aligned Flux
DiffN_aligned_high = Flux_aligned_high.varget('Differential_Number_Flux_Aligned')
DiffN_aligned_low = Flux_aligned_low.varget('Differential_Number_Flux_Aligned')
Flux_aligned_high_info = Flux_aligned_high.cdf_info()
Flux_aligned_low_info = Flux_aligned_low.cdf_info()
zvars_flux_high = Flux_aligned_high_info['zVariables']
zvars_flux_low = Flux_aligned_low_info['zVariables']

# Distribution_Function
dist_high = dist_file_high.varget('Distribution_Function')
dist_low = dist_file_low.varget('Distribution_Function')
dist_high_info = dist_file_high.cdf_info()
dist_low_info = dist_file_low.cdf_info()
zvars_high_dist = dist_high_info['zVariables']
zvars_low_dist = dist_low_info['zVariables']

# Distribution Function for J_parallel_Calculation
dist_jpar_high = dist_file_high.varget('Distribution_Function_jpar')
dist_jpar_low = dist_file_low.varget('Distribution_Function_jpar')


#Get the working information
pitch_actual_high_info = pitch_file_actual_high.cdf_info()
pitch_actual_low_info = pitch_file_actual_low.cdf_info()

# Parallel Current
J_par_high = J_par_file_high.varget('Parallel_Current_Density')
J_par_high_info = J_par_file_high.cdf_info()
zvars_J_par_high = J_par_high_info['zVariables']

J_par_low = J_par_file_low.varget('Parallel_Current_Density')
J_par_low_info = J_par_file_low.cdf_info()
zvars_J_par_low = J_par_low_info['zVariables']




# -----------------
# NOISE INFORMATION
# -----------------
#Ranges for the noise-collection program
Maskvals = [3,4]
noise_threshholds = [0,0]
bin_spacing = 1  #space between the points (in degrees)

t_noise_start1_low = 10089 # corresponds to TRICE_high 08:31:10.014
# t_noise_start1_low = 12889 # corresponds to TRICE_high 08:33:30.01
t_noise_end1_low = 15289 # corresponds to TRICE_high 08:35:30.016
t_noise_start2_low = 22789 # corresponds to TRICE_high 08:41:45.018
t_noise_end2_low = 23589 # corresponds to TRICE_high 08:42:25.18

t_noise_start1_high = 7610# corresponds to TRICE_high 08:28:59.99
# t_noise_start1_high = 13310 # corresponds to TRICE_high 08:33:44.99
t_noise_end1_high = 14090 # corresponds to TRICE_high 08:34:23.99
t_noise_start2_high = 20510 # corresponds to TRICE_high 08:39:44.994
t_noise_end2_high = 22410 # corresponds to TRICE_high 08:41:19.995


roll_reduced_low = sun_spike_noise_file_low.varget('Roll_Reduced')
roll_reduced_high = sun_spike_noise_file_high.varget('Roll_Reduced')

roll_epoch_low = sun_spike_noise_file_low.varget('Roll_Epoch')
roll_epoch_high = sun_spike_noise_file_high.varget('Roll_Epoch')

noise_vals_low = sun_spike_noise_file_low.varget('noise_vals')
noise_vals_high = sun_spike_noise_file_high.varget('noise_vals')

masked_counts_low = sun_spike_noise_file_low.varget('counts_masked')
masked_counts_high = sun_spike_noise_file_high.varget('counts_masked')


#Noise Statistics
spike_statistics_means_low = spike_rem_file_statistics_low.varget('noise_means')
spike_statistics_stddevs_low = spike_rem_file_statistics_low.varget('noise_stds')
spike_statistics_means_high = spike_rem_file_statistics_high.varget('noise_means')
spike_statistics_stddevs_high = spike_rem_file_statistics_high.varget('noise_stds')


# #Noise Reduced COUNTS
# rem_sun_spike_counts_low = spike_rem_counts_low_file.varget('removed_spike_counts_threshed_data')
# rem_sun_spike_counts_high = spike_rem_counts_high_file.varget('removed_spike_counts_threshed_data')


#Ranges for data removal
# LOW_ranges = [[-82,-57],[97,125]]
# HIGH_ranges = [[-93,-59],[92,125]]

#threshold value to determine noise
hist_val_engy = 24
hist_val_threshold = 4

#------------------------
# General EEPAA Variables
# -----------------------

# Fill Values
placeholder = EEPAA_file_high.varattsget(zvars_high[1])
placeholder1 = placeholder['FILLVAL']
fillval_high = placeholder1[0]

placeholder = EEPAA_file_low.varattsget(zvars_low[1])
placeholder1 = placeholder['FILLVAL']
fillval_low = placeholder1[0]


ranges_high = [
    range(len(EPOCH_High)),
    range(21),
    range(49)]

ranges_low = [
    range(len(EPOCH_low)),
    range(21),
    range(49)]

ranges_counts_high = [
    range(len(EPOCH_high_counts)),
    range(21),
    range(49)
    ]

ranges_counts_low = [
    range(len(EPOCH_low_counts)),
    range(21),
    range(49)]

t_start_high = 14200 # corresponds to TRICE_high 08:34:29.493
t_end_high = 20810 +1# corresponds to TRICE_high 08:39:59.995
chi_ranges_high = [
range(t_start_high,t_end_high),
range(15,49)]

t_start_low = 15769 # corresponds to TRICE_low 08:35:54.016
t_end_low = 22729 +1# corresponds to TRICE_low 08:41:42.018
chi_ranges_low = [
range(t_start_low,t_end_low),
range(15,49)]

# Random things
J_high_tag = 'J_parallel_high'
J_low_tag = 'J_parallel_low'

rocket_names = ['LOW','HIGH']




