# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ------- TRICEII Calculate DiffFlux code------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


print('importing variables: ', end='')

import time
start_time = time.time()

import numpy as np
from cdflib import cdfwrite
import itertools

# Get the files
from files import cal_file_counts_high,cal_file_counts_low,root,EEPAA_file_high,EEPAA_file_low
from Variables import cal_counts_high,cal_counts_low
from Variables import ranges_high,ranges_low,Energy_High,Energy_low,EEPAA_file_high_info,EEPAA_file_low_info
from Variables import geometric_high,geometric_low,count_interval_high,count_interval_low,deadtime_high,deadtime_low,zvars_high,zvars_low,DiffEFlux_high


diff_file_high = cdfwrite.CDF(root + "/" + 'diff_flux_high_cal' + '.cdf',cdf_spec=EEPAA_file_high_info,delete=True)
diff_file_low = cdfwrite.CDF(root + "/" + 'diff_flux_low_cal' + '.cdf',cdf_spec=EEPAA_file_low_info,delete=True)


diffNflux_high = np.ndarray(shape=(len(cal_counts_high),len(cal_counts_high[0]),len(cal_counts_high[0][0])), dtype='float64')
diffEflux_high = np.ndarray(shape=(len(cal_counts_high),len(cal_counts_high[0]),len(cal_counts_high[0][0])), dtype='float64')
diffNflux_low = np.ndarray(shape=(len(cal_counts_low),len(cal_counts_low[0]),len(cal_counts_low[0][0])), dtype='float64')
diffEflux_low = np.ndarray(shape=(len(cal_counts_low),len(cal_counts_low[0]),len(cal_counts_low[0][0])), dtype='float64')


#Find max/min of data
minE_low = 1000000000
maxE_low = 10
minN_low = 1000000000
maxN_low = 10
minE_high = 1000000000
maxE_high = 10
minN_high = 1000000000
maxN_high = 10

print('Done')
print('Calculating diffFluxs for EEPAA HIGH: ', end='')
for tme,ptch,engy in itertools.product(*ranges_high):
    diffEflux_high[tme][ptch][engy] = (Energy_High[engy] * cal_counts_high[tme][ptch][engy]) / (Energy_High[engy] * geometric_high[ptch] * (count_interval_high[0] - (cal_counts_high[tme][ptch][engy] * deadtime_high)))
    diffNflux_high[tme][ptch][engy] = (cal_counts_high[tme][ptch][engy]) / (Energy_High[engy] * geometric_high[ptch] * (count_interval_high[0] - (cal_counts_high[tme][ptch][engy] * deadtime_high)))
    if diffEflux_high[tme][ptch][engy] > maxE_high:
        maxE_high =diffEflux_high[tme][ptch][engy]
    if diffEflux_high[tme][ptch][engy] < minE_high:
        minE_high = diffEflux_high[tme][ptch][engy]
    if diffNflux_high[tme][ptch][engy] > maxN_high:
        maxN_high =diffNflux_high[tme][ptch][engy]
    if diffNflux_high[tme][ptch][engy] < minN_high:
        minN_high = diffNflux_high[tme][ptch][engy]

print('Done')
print('Calculating diffFluxs for EEPAA LOW: ', end='')
for tme,ptch,engy in itertools.product(*ranges_low):
    diffEflux_low[tme][ptch][engy] = (Energy_low[engy] * cal_counts_low[tme][ptch][engy]) / (Energy_low[engy] * geometric_low[ptch] * (count_interval_low[0] - (cal_counts_low[tme][ptch][engy] * deadtime_low)))
    diffNflux_low[tme][ptch][engy] = (cal_counts_low[tme][ptch][engy]) / (Energy_low[engy] * geometric_low[ptch] * (count_interval_low[0] - (cal_counts_low[tme][ptch][engy] * deadtime_low)))
    if diffEflux_low[tme][ptch][engy] > maxE_low:
        maxE_low =diffEflux_low[tme][ptch][engy]
    if diffEflux_low[tme][ptch][engy] < minE_low:
        minE_low = diffEflux_low[tme][ptch][engy]
    if diffNflux_low[tme][ptch][engy] > maxN_low:
        maxN_low =diffNflux_low[tme][ptch][engy]
    if diffNflux_low[tme][ptch][engy] < minN_low:
        minN_low = diffNflux_low[tme][ptch][engy]
print('Done')

# Look at the min/max's
# print(maxE_low,maxN_low,minE_low,minN_low)



#####OUTPUT DISTRIBUTION DATA#####

#HIGH

#DiffNFlux
vardata= diffNflux_high[:-1]

varinfo = EEPAA_file_high.varinq(zvars_high[-4])
varinfo['Variable'] = 'Differential_Number_Flux'

varattrs_1 = EEPAA_file_high.varattsget(zvars_high[-4], expand=True)
varattrs_1['CATDESC'] = 'Differential_Number_Flux'
varattrs_1['FIELDNAM'] = 'Differential_Number_Flux'
varattrs_1['UNITS'] = 'counts!cm!U-2!N str!U-1!N s!U-1!N eV!U-1!N'
varattrs_1['SCALETYP'] = 'log'
varattrs_1['VALIDMIN'] = [maxN_high]
varattrs_1['VALIDMAX'] = [minN_high]
varattrs_1['LABLAXIS'] = 'Differential_Number_Flux'
varattrs_1['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_1['LABL_PTR_2'] = 'LABL_PTR_2'

diff_file_high.write_var(varinfo,var_attrs=varattrs_1,var_data=vardata)



#DiffEFlux
vardata= diffEflux_high[:-1]

varinfo = EEPAA_file_high.varinq(zvars_high[6])
varinfo['Variable'] = 'Differential_Energy_Flux'

varattrs_1 = EEPAA_file_high.varattsget(zvars_high[6], expand=True)
varattrs_1['CATDESC'] = 'Differential_Energy_Flux'
varattrs_1['FIELDNAM'] = 'Differential_Energy_Flux'
varattrs_1['UNITS'] = 'cm!U-2!N str!U-1!N s!U-1!N eV/eV'
varattrs_1['SCALETYP'] = 'log'
varattrs_1['VALIDMIN'] = [maxE_high]
varattrs_1['VALIDMAX'] = [minE_high]
varattrs_1['LABLAXIS'] = 'Differential_Energy_Flux'
varattrs_1['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_1['LABL_PTR_2'] = 'LABL_PTR_2'

diff_file_high.write_var(varinfo,var_attrs=varattrs_1,var_data=vardata)



#LOW
vardata= diffNflux_low[:-1]

varinfo = EEPAA_file_low.varinq(zvars_low[-4])
varinfo['Variable'] = 'Differential_Number_Flux'

varattrs_2 = EEPAA_file_low.varattsget(zvars_low[-4], expand=True)
varattrs_2['CATDESC'] = 'Differential_Number_Flux'
varattrs_2['FIELDNAM'] = 'Differential_Number_Flux'
varattrs_2['UNITS'] = 'N (counts)!cm!U-2!N str!U-1!N s!U-1!N eV!U-1!N'
varattrs_2['SCALETYP'] = 'log'
varattrs_2['VALIDMIN'] = [maxN_low]
varattrs_2['VALIDMAX'] = [minN_low]
varattrs_2['LABLAXIS'] = 'Differential_Number_Flux'
varattrs_2['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_2['LABL_PTR_2'] = 'LABL_PTR_2'

diff_file_low.write_var(varinfo,var_attrs=varattrs_2,var_data=vardata)


#DiffEFlux
vardata= diffEflux_low[:-1]

varinfo = EEPAA_file_low.varinq(zvars_low[6])
varinfo['Variable'] = 'Differential_Energy_Flux'

varattrs_2 = EEPAA_file_low.varattsget(zvars_low[6], expand=True)
varattrs_2['CATDESC'] = 'Differential_Energy_Flux'
varattrs_2['FIELDNAM'] = 'Differential_Energy_Flux'
varattrs_2['UNITS'] = 'cm!U-2!N str!U-1!N s!U-1!N eV/eV'
varattrs_2['SCALETYP'] = 'log'
varattrs_2['VALIDMIN'] = [maxE_low]
varattrs_2['VALIDMAX'] = [minE_low]
varattrs_2['LABLAXIS'] = 'Differential_Energy_Flux'
varattrs_2['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_2['LABL_PTR_2'] = 'LABL_PTR_2'

diff_file_low.write_var(varinfo,var_attrs=varattrs_2,var_data=vardata)









#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cal_file_counts_high.close()
cal_file_counts_low.close()
EEPAA_file_high.close()
EEPAA_file_low.close()
diff_file_high.close()
diff_file_low.close()


print("--- %s seconds for Calc_Differential_Flux---" % (time.time() - start_time) ,'\n')


