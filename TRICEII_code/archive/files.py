from cdflib import cdfread
from os import path

# ---------------
# FILE REPOSITORY
# ---------------
# root = "C:/Users/Rundus/Desktop/TRICE/" #My home PC
root = 'C:/Users/cfeltman/Desktop/TRICE/' # My labtop
output_path = 'Output/'
Mag_file_high = cdfread.CDF( root + "mag_52003_cal_v2.cdf")
Mag_file_low = cdfread.CDF(root + "mag_52004_cal_v2.cdf")
EEPAA_file_high = cdfread.CDF(root + "TRICE_52003_l2_eepaa_20181208T082239_v1.1.2.cdf")
EEPAA_file_low = cdfread.CDF(root + "TRICE_52004_l2_eepaa_20181208T082243_v1.1.2.cdf")
Flux_aligned_high = cdfread.CDF(root + "DiffFlux_sorted_high.cdf")
Flux_aligned_low = cdfread.CDF(root + "DiffFlux_sorted_low.cdf")
dist_file_high = cdfread.CDF(root + "dist_high.cdf")
dist_file_low = cdfread.CDF(root + "dist_low.cdf")
pitch_file_actual_high = cdfread.CDF(root + 'pitch_actual_high.cdf')
pitch_file_actual_low = cdfread.CDF(root + 'pitch_actual_low.cdf')

# Raw Data Counts
counts_file_high = cdfread.CDF(root + 'TRICE_52003_l1_eepaa_20181208T082239_v1.1.2_COUNTS.cdf')
counts_file_low = cdfread.CDF(root + 'TRICE_52004_l1_eepaa_20181208T082243_v1.1.2_COUNTS.cdf')

#Sun-Spike Removed Counts
spike_rem_counts_high_file = cdfread.CDF(root + 'spike_rem_counts_high')
spike_rem_counts_low_file = cdfread.CDF(root + 'spike_rem_counts_low')
cal_diff_file_high = cdfread.CDF(root + 'diff_flux_high_cal' + '.cdf')
cal_diff_file_low = cdfread.CDF(root + 'diff_flux_low_cal'+'.cdf')

cal_file_counts_high = cdfread.CDF(root + 'spike_rem_counts_high' + '.cdf')
cal_file_counts_low = cdfread.CDF(root + 'spike_rem_counts_low' + '.cdf')




J_par_file_high = cdfread.CDF(root + 'J_parallel_high.cdf')
J_par_file_low = cdfread.CDF(root + 'J_parallel_low' + '.cdf')

attitude_control_file_high = cdfread.CDF(root + output_path +'TRICE_52003_20181208T082239_attitude_control&_v1.1' + '.cdf')
attitude_control_file_low = cdfread.CDF(root + output_path + 'TRICE_52004_20181208T082243_attitude_control&_v1.1'+ '.cdf')

sun_spike_noise_file_high = cdfread.CDF(root + 'sun_spike_noise_dat_high' + '.cdf')
sun_spike_noise_file_low = cdfread.CDF(root + 'sun_spike_noise_dat_low' + '.cdf')
spike_rem_file_statistics_high = cdfread.CDF(root + 'spike_rem_statistics_high' + '.cdf')
spike_rem_file_statistics_low = cdfread.CDF(root + 'spike_rem_statistics_low' + '.cdf')


