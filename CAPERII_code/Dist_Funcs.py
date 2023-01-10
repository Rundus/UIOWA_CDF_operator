# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE Distribution Funct. output/sorting Code--------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



import numpy as np, time
from cdflib import cdfread, cdfwrite

start_time = time.time()
print('importing variables: ', end='')


# ------------
# Import Files
# ------------


from files import flux_aligned_file_1,flux_aligned_file_2,EEPAA_file_1,EEPAA_file_2,root,user_path
from Variables import EEPAA_file_1_info,EEPAA_file_2_info,flux_aligned_1_info,flux_aligned_2_info,zvars_flux_aligned_1,zvars_flux_aligned_2,zvars_1,zvars_2,fillval_1,fillval_2,Nflux_aligned_1,Nflux_aligned_2,energies_1,energies_2,epoch_1,epoch_2,ranges_1,ranges_2


# ------------
# OUTPUT FILES
# ------------
dist_1_tag = 'dist_1'
dist_2_tag = 'dist_2'
dist_1 = cdfwrite.CDF(user_path + root + "/" + dist_1_tag + '.cdf',cdf_spec=EEPAA_file_1_info,delete=True)
dist_2 = cdfwrite.CDF(user_path + root + "/" + dist_2_tag + '.cdf',cdf_spec=EEPAA_file_2_info,delete=True)

print('Done')


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -----PRODUCE DISTRIBUTION FUNCTION DATA-----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Define the constants for the distribution function
m = 9.11*10**(-31)
q = 1.602176565 * 10**(-19)
scale = 100

print('Calculating the Distribution Functions for eepaa1: ',end='')

dist_dat_1 = np.ndarray(shape=(len(epoch_1),21,49),dtype='float64')
dist_dat_2 = np.ndarray(shape=(len(epoch_2),21,49),dtype='float64')

max1 = 0
min1 = 10000000000000

dist_time1 = time.time()

# eepaa1 Distribution Function
for tme in ranges_1[0]:
    for ptch in ranges_1[1]:
        dist_dat_1[tme][ptch] = [ (((0.5*((scale*m)**2) * Nflux_aligned_1[tme][ptch][engy]) *( (energies_1[engy])*q*q)**(-1))) for engy in ranges_1[2]]
        for engy in ranges_1[2]:
            if dist_dat_1[tme][ptch][engy] != fillval_1 and dist_dat_1[tme][ptch][engy] >=0:
                if dist_dat_1[tme][ptch][engy] > max1:
                    max1 = dist_dat_1[tme][ptch][engy]
                if dist_dat_1[tme][ptch][engy] < min1:
                    min1 = dist_dat_1[tme][ptch][engy]
            elif dist_dat_1[tme][ptch][engy] != fillval_1 and dist_dat_1[tme][ptch][engy] < 0:
                dist_dat_1[tme][ptch][engy] = fillval_1

print('Done')

print("--- %s seconds for Distribution HIGH---" % (time.time() - dist_time1) ,'\n')


max2 = 0
min2 = 10000000000000

print('Calculating the Distribution Functions for eepaa2: ',end='')
dist_time2 = time.time()


# eepaa2 Distribution Function
for tme in ranges_2[0]:
    for ptch in ranges_2[1]:
        dist_dat_2[tme][ptch] = [ (((0.5*((scale*m)**2) * Nflux_aligned_2[tme][ptch][engy]) *( (energies_2[engy])*q*q)**(-1))) for engy in ranges_1[2]]
        for engy in ranges_2[2]:
            if dist_dat_2[tme][ptch][engy] != fillval_2 and dist_dat_2[tme][ptch][engy] >=0:
                if dist_dat_2[tme][ptch][engy] > max2:
                    max2 = dist_dat_2[tme][ptch][engy]
                if dist_dat_2[tme][ptch][engy] < min2:
                    min2 = dist_dat_2[tme][ptch][engy]
            elif dist_dat_2[tme][ptch][engy] != fillval_2 and dist_dat_2[tme][ptch][engy] < 0:
                dist_dat_2[tme][ptch][engy] = fillval_2

print('Done')

print("--- %s seconds for Distribution LOW---" % (time.time() - dist_time2) ,'\n')



#---------------------------------
#####OUTPUT DISTRIBUTION DATA#####
#---------------------------------

#eepaa1
tag = 18
vardata= dist_dat_1
varinfo = EEPAA_file_1.varinq(zvars_1[tag])
varinfo['Variable'] = 'Distribution_Function'
varattrs_high = EEPAA_file_1.varattsget(zvars_1[tag], expand=True)
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
dist_1.write_var(varinfo,var_attrs=varattrs_high,var_data=vardata)




#eepaa2
vardata= dist_dat_2
varinfo = EEPAA_file_2.varinq(zvars_2[tag])
varinfo['Variable'] = 'Distribution_Function'
varattrs_low = EEPAA_file_2.varattsget(zvars_2[tag], expand=True)
varattrs_low['CATDESC'] = 'Distribution_Function'
varattrs_low['FIELDNAM'] = 'Distribution_Function'
varattrs_low['UNITS'] = 'Counts!U 1!N m!U-6!N s!U-3!'
varattrs_low['SCALETYP'] = 'log'
varattrs_low['VALIDMIN'] = [min2]
varattrs_low['VALIDMAX'] = [max2]
varattrs_low['LABLAXIS'] = 'Distribution_Function'
varattrs_low['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_low['LABL_PTR_2'] = 'LABL_PTR_2'
dist_2.write_var(varinfo,var_attrs=varattrs_low,var_data=vardata)







#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flux_aligned_file_1.close()
flux_aligned_file_2.close()
EEPAA_file_1.close()
EEPAA_file_2.close()
dist_1.close()
dist_2.close()


print("--- %s seconds for Dist_Funcs---" % (time.time() - start_time) ,'\n')