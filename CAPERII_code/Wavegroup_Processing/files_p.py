from cdflib import cdfread


#Get the files for EEPAA
root = 'CAPERII'
user_path = "C:/Users/cfeltman/Desktop/"
mag = 'CAPERII_52005_l0_mag_20190104_v01'
eepaa1 = 'CAPERII_52005_l1_eepaa_20190104_v1.1.2'
eepaa2 = 'CAPERII_52005_l1_eepaa2_20190104_v1.1.2'
correlator = 'CAPERII_52005_l0_correlator_20190104_v1.2.1'
mag_file = cdfread.CDF(user_path + root + "/data/" + mag)
EEPAA_file_1 = cdfread.CDF(user_path + root + "/data/" + eepaa1 + '.cdf')
EEPAA_file_2 = cdfread.CDF(user_path + root + "/data/" + eepaa2 + '.cdf')

#Get the files for EEPAA
files = [user_path + root + "/" + eepaa1 + '.cdf',user_path + root + "/" + eepaa2 + '.cdf',user_path + root + "/" + mag]

# Get Spencer's Files
file_freqs = user_path + 'DC_Spectral/FFT_freqs__1024pt__490-577s.txt'
file_LH_mag = user_path + 'DC_Spectral/FFT_LH_magnitudes__1024pt__490-577s.txt'
file_LR_ratio = user_path +' +DC_Spectral/FFT_LR_ratios__1024pt__490-577s.txt'
file_RH_mag = user_path + 'DC_Spectral/FFT_RH_magnitudes__1024pt__490-577s.txt'
file_times = user_path +'DC_Spectral/FFT_times__1024pt__490-577s.txt'
file_integrate = user_path + "/DC_Spectral/Integrated_PSDs__1024pt__490-577s.txt"
correlator_file = cdfread.CDF(user_path + root +'/Data/' + correlator + '.cdf')