# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ------- CAPERII Calculate DiffFluc code------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


from cdflib import cdfwrite
import itertools, numpy as np,time

print('importing variables: ', end='')

start_time = time.time()

# Get the files
from files import EEPAA_file_1,EEPAA_file_2,root,user_path
from Variables import zvars_1,zvars_2,EEPAA_file_1_info as info_1,EEPAA_file_2_info as info_2, counts_1,counts_2,energies_1,energies_2,geometric_1,geometric_2,count_interval_1,count_interval_2,ranges_1,ranges_2,TRICE_deadtime


print('Done')

# Crete output file
diff_flux_1_tag = 'diff_flux_1'
diff_flux_2_tag = 'diff_flux_2'

diff_file_1 = cdfwrite.CDF(user_path + root + "/" + diff_flux_1_tag + '.cdf',cdf_spec=info_1,delete=True)
diff_file_2 = cdfwrite.CDF(user_path + root + "/" + diff_flux_2_tag + '.cdf',cdf_spec=info_2,delete=True)


diffNflux_1 = np.ndarray(shape=(len(counts_1),len(counts_1[0]),len(counts_1[0][0])), dtype='float64')
diffEflux_1 = np.ndarray(shape=(len(counts_1),len(counts_1[0]),len(counts_1[0][0])), dtype='float64')
diffNflux_2 = np.ndarray(shape=(len(counts_2),len(counts_2[0]),len(counts_2[0][0])), dtype='float64')
diffEflux_2 = np.ndarray(shape=(len(counts_2),len(counts_2[0]),len(counts_2[0][0])), dtype='float64')



print('Calculating diffFluxs for EEPAA1: ', end='')

# EEPAA1
for time_index,pitch_index,energy_index in itertools.product(*ranges_1):
    diffEflux_1[time_index][pitch_index][energy_index] = (energies_1[energy_index] * counts_1[time_index][pitch_index][energy_index]) /(energies_1[energy_index]*geometric_1[pitch_index]*(count_interval_1[time_index]/1000 - (counts_1[time_index][pitch_index][energy_index] * TRICE_deadtime)) )
    diffNflux_1[time_index][pitch_index][energy_index] = (counts_1[time_index][pitch_index][energy_index]) /(energies_1[energy_index]*geometric_1[pitch_index]*(count_interval_1[time_index]/1000 - (counts_1[time_index][pitch_index][energy_index] * TRICE_deadtime)) )

print('Done')
print('Calculating diffFluxs for EEPAA2: ', end='')
# EEPAA2
for time_index,pitch_index,energy_index in itertools.product(*ranges_2):
    diffEflux_2[time_index][pitch_index][energy_index] = (energies_2[energy_index] * counts_2[time_index][pitch_index][energy_index]) /(energies_2[energy_index]*geometric_2[pitch_index]*(count_interval_2[time_index]/1000 - (counts_2[time_index][pitch_index][energy_index] * TRICE_deadtime)) )
    diffNflux_2[time_index][pitch_index][energy_index] = (counts_2[time_index][pitch_index][energy_index]) /(energies_2[energy_index]*geometric_2[pitch_index]*(count_interval_2[time_index]/1000 - (counts_2[time_index][pitch_index][energy_index] * TRICE_deadtime)) )

print('Done')

#NEED TO FIND MIN/MAX OF DIFFEFLUX and DIFFNFLUX
min1 = 0
max1 = 1000000000
min2 = 0
max2 = 1000000000


#####OUTPUT DISTRIBUTION DATA#####


#DiffNFlux
vardata= diffNflux_1

varinfo = EEPAA_file_1.varinq(zvars_1[-5])
varinfo['Variable'] = 'Differential_Number_Flux'

varattrs_1 = EEPAA_file_1.varattsget(zvars_1[-5], expand=True)
varattrs_1['CATDESC'] = 'Differential_Number_Flux'
varattrs_1['FIELDNAM'] = 'Differential_Number_Flux'
varattrs_1['UNITS'] = 'counts!cm!U-2!N str!U-1!N s!U-1!N eV!U-1!N'
varattrs_1['SCALETYP'] = 'log'
varattrs_1['VALIDMIN'] = [min1]
varattrs_1['VALIDMAX'] = [max1]
varattrs_1['LABLAXIS'] = 'Differential_Number_Flux'
#Things that need changing
# varattrs_high['LABL_PTR_1'] = Label_1
# varattrs_high['LABL_PTR_2'] = Label_2
varattrs_1['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_1['LABL_PTR_2'] = 'LABL_PTR_2'

diff_file_1.write_var(varinfo,var_attrs=varattrs_1,var_data=vardata)



#DiffEFlux
vardata= diffEflux_1

varinfo = EEPAA_file_1.varinq(zvars_1[-5])
varinfo['Variable'] = 'Differential_Energy_Flux'

varattrs_1 = EEPAA_file_1.varattsget(zvars_1[-5], expand=True)
varattrs_1['CATDESC'] = 'Differential_Energy_Flux'
varattrs_1['FIELDNAM'] = 'Differential_Energy_Flux'
varattrs_1['UNITS'] = 'cm!U-2!N str!U-1!N s!U-1!N eV/eV'
varattrs_1['SCALETYP'] = 'log'
varattrs_1['VALIDMIN'] = [min1]
varattrs_1['VALIDMAX'] = [max1]
varattrs_1['LABLAXIS'] = 'Differential_Energy_Flux'
#Things that need changing
# varattrs_high['LABL_PTR_1'] = Label_1
# varattrs_high['LABL_PTR_2'] = Label_2
varattrs_1['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_1['LABL_PTR_2'] = 'LABL_PTR_2'

diff_file_1.write_var(varinfo,var_attrs=varattrs_1,var_data=vardata)





#eepaa 2
vardata= diffNflux_2

varinfo = EEPAA_file_2.varinq(zvars_2[-5])
varinfo['Variable'] = 'Differential_Number_Flux'

varattrs_2 = EEPAA_file_2.varattsget(zvars_2[-5], expand=True)
varattrs_2['CATDESC'] = 'Differential_Number_Flux'
varattrs_2['FIELDNAM'] = 'Differential_Number_Flux'
varattrs_2['UNITS'] = 'N (counts)!cm!U-2!N str!U-1!N s!U-1!N eV!U-1!N'
varattrs_2['SCALETYP'] = 'log'
varattrs_2['VALIDMIN'] = [min2]
varattrs_2['VALIDMAX'] = [max2]
varattrs_2['LABLAXIS'] = 'Differential_Number_Flux'
#Things that need changing
# varattrs_high['LABL_PTR_1'] = Label_1
# varattrs_high['LABL_PTR_2'] = Label_2
varattrs_2['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_2['LABL_PTR_2'] = 'LABL_PTR_2'

diff_file_2.write_var(varinfo,var_attrs=varattrs_2,var_data=vardata)


#DiffEFlux
vardata= diffEflux_2

varinfo = EEPAA_file_2.varinq(zvars_2[-5])
varinfo['Variable'] = 'Differential_Energy_Flux'

varattrs_2 = EEPAA_file_2.varattsget(zvars_2[-5], expand=True)
varattrs_2['CATDESC'] = 'Differential_Energy_Flux'
varattrs_2['FIELDNAM'] = 'Differential_Energy_Flux'
varattrs_2['UNITS'] = 'cm!U-2!N str!U-1!N s!U-1!N eV/eV'
varattrs_2['SCALETYP'] = 'log'
varattrs_2['VALIDMIN'] = [min1]
varattrs_2['VALIDMAX'] = [max1]
varattrs_2['LABLAXIS'] = 'Differential_Energy_Flux'
#Things that need changing
# varattrs_high['LABL_PTR_1'] = Label_1
# varattrs_high['LABL_PTR_2'] = Label_2
varattrs_2['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_2['LABL_PTR_2'] = 'LABL_PTR_2'

diff_file_2.write_var(varinfo,var_attrs=varattrs_2,var_data=vardata)









#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EEPAA_file_1.close()
EEPAA_file_2.close()
diff_file_1.close()
diff_file_2.close()



print("--- %s seconds for Calc_Differential_Flux---" % (time.time() - start_time) ,'\n')




