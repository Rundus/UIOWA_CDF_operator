# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE Distribution Funct. output/sorting Code--------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np, time
from cdflib import cdfwrite
import itertools

start_time = time.time()
print('importing variables: ', end='')

# ------------
# Import Files
# ------------

from files import Flux_aligned_high,Flux_aligned_low,EEPAA_file_high,EEPAA_file_low
from Variables import EEPAA_file_low_info,EEPAA_file_high_info,zvars_low,zvars_high,ranges_high,ranges_low
from Variables import EPOCH_low,EPOCH_High,fillval_low,fillval_high,DiffN_aligned_low,DiffN_aligned_high,Fillhigh,Filllow,Energy_High,Energy_low


# ------------
# OUTPUT FILES
# ------------
dist_high = cdfwrite.CDF('C:/Users/Rundus/Desktop/TRICE/dist_high',cdf_spec=EEPAA_file_high_info,delete=True)
dist_low = cdfwrite.CDF('C:/Users/Rundus/Desktop/TRICE/dist_low',cdf_spec=EEPAA_file_low_info,delete=True)

print('Done')




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -----PRODUCE DISTRIBUTION FUNCTION DATA-----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Define the constants for the distribution function
m = 9.11*10**(-31)
q = 1.602176565 * 10**(-19)
scale = 100

print('Calculating the Distribution Functions for HIGH: ',end='')




Dist_Funct_high = np.ndarray(shape=(len(EPOCH_High),21,49),dtype='float64')
Dist_Funct_low = np.ndarray(shape=(len(EPOCH_low),21,49),dtype='float64')


max1 = 0
min1 = 10000000000000

dist_time1 = time.time()

# High Distribution Function
for tme in range(len(EPOCH_High)):
    for ptch in range(21):
        Dist_Funct_high[tme][ptch] = [ (((0.5*((scale*m)**2) * DiffN_aligned_high[tme][ptch][engy]) *( (Energy_High[engy])*q*q)**(-1))) for engy in range(49)]
        for engy in range(49):
            if Dist_Funct_high[tme][ptch][engy] != Fillhigh and Dist_Funct_high[tme][ptch][engy] >=0:
                if Dist_Funct_high[tme][ptch][engy] > max1:
                    max1 = Dist_Funct_high[tme][ptch][engy]
                if Dist_Funct_high[tme][ptch][engy] < min1:
                    min1 = Dist_Funct_high[tme][ptch][engy]
            elif Dist_Funct_high[tme][ptch][engy] != Fillhigh and Dist_Funct_high[tme][ptch][engy] < 0:
                Dist_Funct_high[tme][ptch][engy] = fillval_high

print('Done')

print("--- %s seconds for Distribution HIGH---" % (time.time() - dist_time1) ,'\n')


max2 = 0
min2 = 10000000000000

print('Calculating the Distribution Functions for LOW: ',end='')
dist_time2 = time.time()


# Low Distribution Function
for tme in range(len(EPOCH_low)):
    for ptch in range(21):
        Dist_Funct_low[tme][ptch] = [ (((0.5*((scale*m)**2) * DiffN_aligned_low[tme][ptch][engy]) * ((Energy_low[engy])*q*q   )**(-1)) ) for engy in range(49)]
        for engy in range(49):
            if Dist_Funct_low[tme][ptch][engy] != Filllow and Dist_Funct_low[tme][ptch][engy] >=0:
                if Dist_Funct_low[tme][ptch][engy] > max2:
                    max2 = Dist_Funct_low[tme][ptch][engy]
                if Dist_Funct_low[tme][ptch][engy] < min2:
                    min2 = Dist_Funct_low[tme][ptch][engy]
            elif Dist_Funct_low[tme][ptch][engy] != Filllow and Dist_Funct_low[tme][ptch][engy] < 0:
                Dist_Funct_low[tme][ptch][engy] = fillval_low

print('Done')

print("--- %s seconds for Distribution LOW---" % (time.time() - dist_time2),'\n')



print('Creating Distribution function for J_parallel Calculation')

Dist_j_parallel_high = Dist_Funct_high
Dist_j_parallel_low = Dist_Funct_low

for tme, ptch, engy in itertools.product(*ranges_high):
    if Dist_j_parallel_high[tme][ptch][engy] == fillval_high:
        Dist_j_parallel_high[tme][ptch][engy] = 0

for tme, ptch, engy in itertools.product(*ranges_low):
    if Dist_j_parallel_low[tme][ptch][engy] == fillval_low:
        Dist_j_parallel_low[tme][ptch][engy] = 0

print('Done')

#####OUTPUT DISTRIBUTION DATA#####


# High
vardata= Dist_Funct_high

varinfo = EEPAA_file_high.varinq(zvars_high[6])
varinfo['Variable'] = 'Distribution_Function'

varattrs_high = EEPAA_file_high.varattsget(zvars_high[6], expand=True)
varattrs_high['CATDESC'] = 'Distribution_Function'
varattrs_high['FIELDNAM'] = 'Distribution_Function'
varattrs_high['UNITS'] = 'Counts!U 1!N m!U-6!N s!U3!'
varattrs_high['SCALETYP'] = 'log'
varattrs_high['VALIDMIN'] = [min1]
varattrs_high['VALIDMAX'] = [max1]
varattrs_high['LABLAXIS'] = 'Distribution_Function'
#Things that need changing
# varattrs_high['LABL_PTR_1'] = Label_1
# varattrs_high['LABL_PTR_2'] = Label_2
varattrs_high['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_high['LABL_PTR_2'] = 'LABL_PTR_2'

dist_high.write_var(varinfo,var_attrs=varattrs_high,var_data=vardata)


#Write out the J_par_dist
vardata = Dist_j_parallel_high
varinfo['Variable'] = 'Distribution_Function_jpar'
varattrs_high['CATDESC'] = 'Distribution_Function_jpar'
varattrs_high['FIELDNAM'] = 'Distribution_Function_jpar'
dist_high.write_var(varinfo,var_attrs=varattrs_high,var_data=vardata)



#Low
vardata= Dist_Funct_low

varinfo = EEPAA_file_low.varinq(zvars_high[6])
varinfo['Variable'] = 'Distribution_Function'

varattrs_low = EEPAA_file_low.varattsget(zvars_low[6], expand=True)
varattrs_low['CATDESC'] = 'Distribution_Function'
varattrs_low['FIELDNAM'] = 'Distribution_Function'
varattrs_low['UNITS'] = 'Counts!U 1!N m!U-6!N s!U-3!'
varattrs_low['SCALETYP'] = 'log'
varattrs_low['VALIDMIN'] = [min2]
varattrs_low['VALIDMAX'] = [max2]
varattrs_low['LABLAXIS'] = 'Distribution_Function'
varattrs_low['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_low['LABL_PTR_2'] = 'LABL_PTR_2'

dist_low.write_var(varinfo,var_attrs=varattrs_low,var_data=vardata)

#Write out the
vardata = Dist_j_parallel_low
varinfo['Variable'] = 'Distribution_Function_jpar'
varattrs_low['CATDESC'] = 'Distribution_Function_jpar'
varattrs_low['FIELDNAM'] = 'Distribution_Function_jpar'
dist_low.write_var(varinfo,var_attrs=varattrs_low,var_data=vardata)

#
# print("\nDist one max: ",max1,'\nDist one min: ',min1,'\n')
# print("\nDist two max: ",max2,'\nDist two min: ',min2,'\n')






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Flux_aligned_high.close()
Flux_aligned_low.close()
EEPAA_file_high.close()
EEPAA_file_low.close()
dist_high.close()
dist_low.close()


print("--- %s seconds for Dist_Funcs---" % (time.time() - start_time) ,'\n')