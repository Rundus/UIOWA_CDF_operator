from cdflib import cdfread

# ---------------
# FILE REPOSITORY
# ---------------
user_path = "C:/Users/cfeltman/Desktop/"
root = 'CAPERII/data'
output_root ='CAPERII'
mag = 'CAPERII_52005_l0_mag_20190104_v01'
eepaa1 = 'CAPERII_52005_l1_eepaa_20190104_v1.1.2'
eepaa2 = 'CAPERII_52005_l1_eepaa2_20190104_v1.1.2'
correlator = 'CAPERII_52005_l0_correlator_20190104_v1.2.1'

mag_file = cdfread.CDF(user_path + root + "/" + mag)
EEPAA_file_1 = cdfread.CDF(user_path + root + "/" + eepaa1 + '.cdf')
EEPAA_file_2 = cdfread.CDF(user_path + root + "/" + eepaa2 + '.cdf')
diffFlux_file_1 = cdfread.CDF(user_path  + root + "/" + 'diff_flux_1' + '.cdf')
diffFlux_file_2 = cdfread.CDF(user_path + root + "/" + 'diff_flux_2' + '.cdf')
pitch_actual_file_1 = cdfread.CDF(user_path + root + "/" + 'pitch_actual_1' + '.cdf')
pitch_actual_file_2 = cdfread.CDF(user_path + root + "/" + 'pitch_actual_2' + '.cdf')
flux_aligned_file_1 = cdfread.CDF(user_path + root + "/" + 'diffFlux_aligned_1' + '.cdf')
flux_aligned_file_2 = cdfread.CDF(user_path + root + "/" + 'diffFlux_aligned_2' + '.cdf')
# dist_file_1 = cdfread.CDF(user_path + root + "/" + 'dist_1' + '.cdf')
# dist_file_2 = cdfread.CDF(user_path + root + "/" + 'dist_2' + '.cdf')
correlator_file = cdfread.CDF(user_path + root +'/' + correlator + '.cdf')


# EXTRA FILES FOR COMPARISON: TRICE
Mag_file_high = cdfread.CDF("C:/Users/cfeltman/Desktop/TRICE/mag_52003_cal_v2.cdf")
EEPAA_file_low = cdfread.CDF("C:/Users/cfeltman/Desktop/TRICE/TRICE_52004_l2_eepaa_20181208T082243_v1.1.2.cdf")
EEPAA_file_high = cdfread.CDF("C:/Users/cfeltman/Desktop/TRICE/TRICE_52003_l2_eepaa_20181208T082239_v1.1.2.cdf")
counts_file_high = cdfread.CDF('C:/Users/cfeltman/Desktop/TRICE/TRICE_52003_l1_eepaa_20181208T082239_v1.1.2_COUNTS.cdf')
counts_file_low = cdfread.CDF('C:/Users/cfeltman/Desktop/TRICE/TRICE_52004_l1_eepaa_20181208T082243_v1.1.2_COUNTS.cdf')