# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE Sorting the DiffEFlux Variable--------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


print('importing variables: ', end='')


# ------------
# Import Files
# ------------
import numpy as np, itertools,time
start_time = time.time()
from cdflib import cdfwrite
from files import pitch_actual_file_1,pitch_actual_file_2,EEPAA_file_1,EEPAA_file_2,user_path,root
from functions import dataplacer,dataplacer_checker

# ---------
# VARIABLES
# ---------

from Variables import diffEflux_1,diffEflux_2,diffNflux_1,diffNflux_2,pitches_1 as pitch,pitch_angle_actual_1 as pitch_actual_1,pitch_angle_actual_2 as pitch_actual_2,fillval_1,fillval_2,ranges_1,ranges_2,pitch_actual_1_info,pitch_actual_2_info,EEPAA_file_1_info,EEPAA_file_2_info,zvars_pitch_1,zvars_pitch_2,zvars_1,zvars_2

DiffE_list_1 = diffEflux_1.tolist()
DiffN_list_1 = diffNflux_1.tolist()
DiffE_list_2 = diffEflux_2.tolist()
DiffN_list_2 = diffNflux_2.tolist()

ranges_1_ex = [
    ranges_1[0],
    ranges_1[1],
    ranges_1[2]]

ranges_2_ex = [
    ranges_2[0],
    ranges_2[1],
    ranges_2[2]]

# ----------------------
# Output Files/Variables
# ----------------------
sortedflux_1_tag = 'diffFlux_aligned_1'
sortedflux_2_tag ='diffFlux_aligned_2'
Sortedflux_high = cdfwrite.CDF(user_path + root + "/" + sortedflux_1_tag + '.cdf',cdf_spec=EEPAA_file_1_info,delete=True)
Sortedflux_low = cdfwrite.CDF(user_path + root + "/" + sortedflux_2_tag + '.cdf',cdf_spec=EEPAA_file_2_info,delete=True)


print('Done')


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ----------------------Sorting the DiffFlux data ---------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flux_time1 = time.time()
print('Reorganizing the Differential energy fluxes for eepaa1: ', end='')


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ----------- eepaa1 -----------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# If you find a pitch angle actual that's out of place, look through the other pitch angles at the time and energy. Find where it fits between +/-5 and put it there when it does. Special cases are given to check -10 and 190 because these need to be treated like 10 and 170
for time_index,pitch_index,energy_index in itertools.product(*ranges_1_ex):
    # FOR ANGLES 0-180

    if (0 < pitch_index < 20) and ((pitch[pitch_index] - 5) <= pitch_actual_1[time_index][pitch_index][energy_index] <= (pitch[pitch_index] + 5)  ):
        dataplacer(time_index,pitch_index,energy_index,diffEflux_1,diffNflux_1,DiffE_list_1,DiffN_list_1)

    elif (0 < pitch_index < 20) and ( (pitch[pitch_index] - 5) > pitch_actual_1[time_index][pitch_index][energy_index] or  (pitch[pitch_index] + 5) < pitch_actual_1[time_index][pitch_index][energy_index]):
        for checker in range(21):
            if checker == 0 and (    (pitch[2] -5) <= pitch_actual_1[time_index][pitch_index][energy_index] <= (pitch[2] + 5)    ):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_1,diffNflux_1,DiffE_list_1,DiffN_list_1,checker)
            elif checker == 20 and (    (pitch[18] -5) <= pitch_actual_1[time_index][pitch_index][energy_index] <= (pitch[18] + 5)    ):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_1,diffNflux_1,DiffE_list_1,DiffN_list_1,checker)
            elif (0 < checker < 20) and (    (pitch[checker] -5) <= pitch_actual_1[time_index][pitch_index][energy_index] <= (pitch[checker] + 5)    ):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_1,diffNflux_1,DiffE_list_1,DiffN_list_1,checker)


    # FOR ANGLE -10
    if (pitch_index == 0) and ((pitch[2] - 5) <= pitch_actual_1[time_index][pitch_index][energy_index] <= (pitch[2] + 5)  ):
        dataplacer(time_index,pitch_index,energy_index,diffEflux_1,diffNflux_1,DiffE_list_1,DiffN_list_1)  # Place the data into the right spot

    elif (pitch_index == 0) and ( (pitch[2] - 5) > pitch_actual_1[time_index][pitch_index][energy_index] or  (pitch[2] + 5) < pitch_actual_1[time_index][pitch_index][energy_index]):
        for checker in range(21):
            if checker == 0 and ((pitch[2] - 5) <= pitch_actual_1[time_index][pitch_index][energy_index] <= (pitch[2] + 5)):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_1,diffNflux_1,DiffE_list_1,DiffN_list_1,checker)
            elif checker == 20 and ((pitch[18] - 5) <= pitch_actual_1[time_index][pitch_index][energy_index] <= (pitch[18] + 5)):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_1,diffNflux_1,DiffE_list_1,DiffN_list_1,checker)
            elif (0 < checker < 20) and ((pitch[checker] - 5) <= pitch_actual_1[time_index][pitch_index][energy_index] <= (pitch[checker] + 5)):
                # 2) Place the data into the right spot
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_1,diffNflux_1,DiffE_list_1,DiffN_list_1,checker)



    # FOR ANGLE 190
    if (pitch_index == 20) and ((pitch[18] - 5) <= pitch_actual_1[time_index][pitch_index][energy_index] <= (pitch[18] + 5)  ):
        dataplacer(time_index,pitch_index,energy_index,diffEflux_1,diffNflux_1,DiffE_list_1,DiffN_list_1)


    elif (pitch_index == 20) and ( (pitch[18] - 5) > pitch_actual_1[time_index][pitch_index][energy_index] or  (pitch[18] + 5) < pitch_actual_1[time_index][pitch_index][energy_index]):
        for checker in range(21):
            if checker == 0 and ((pitch[2] - 5) <= pitch_actual_1[time_index][pitch_index][energy_index] <= (pitch[2] + 5)):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_1,diffNflux_1,DiffE_list_1,DiffN_list_1,checker)
            elif checker == 20 and ((pitch[18] - 5) <= pitch_actual_1[time_index][pitch_index][energy_index] <= (pitch[18] + 5)):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_1,diffNflux_1,DiffE_list_1,DiffN_list_1,checker)
            elif (0 < checker < 20) and ((pitch[checker] - 5) <= pitch_actual_1[time_index][pitch_index][energy_index] <= (pitch[checker] + 5)):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_1,diffNflux_1,DiffE_list_1,DiffN_list_1,checker)




print('Done')
print('Averaging and filling in differential data for eepaa1: ', end='')

# Removed all data points that are not in lists and replace with fill values and AVERAGE the ones that are lists
for time_index,pitch_index,energy_index in itertools.product(*ranges_1_ex):
    # Remove
    if not isinstance(DiffE_list_1[time_index][pitch_index][energy_index],list):
        DiffE_list_1[time_index][pitch_index][energy_index] = fillval_1
    if not isinstance(DiffN_list_1[time_index][pitch_index][energy_index],list):
        DiffN_list_1[time_index][pitch_index][energy_index] = fillval_1
    # Average
    if isinstance(DiffE_list_1[time_index][pitch_index][energy_index], list):
        DiffE_list_1[time_index][pitch_index][energy_index] = sum(DiffE_list_1[time_index][pitch_index][energy_index])/len(DiffE_list_1[time_index][pitch_index][energy_index])
    if isinstance(DiffN_list_1[time_index][pitch_index][energy_index], list):
        DiffN_list_1[time_index][pitch_index][energy_index] = sum(DiffN_list_1[time_index][pitch_index][energy_index])/len(DiffN_list_1[time_index][pitch_index][energy_index])


print('Done')


print("--- %s seconds for sorting eepaa1---" % (time.time() - flux_time1) ,'\n')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# ----------- LOW -----------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
flux_time2 = time.time()
print('Reorganizing the diffEflux for eepaa2: ', end='')

for time_index,pitch_index,energy_index in itertools.product(*ranges_2_ex):
    # FOR ANGLES 0-180

    if (0 < pitch_index < 20) and ((pitch[pitch_index] - 5) <= pitch_actual_2[time_index][pitch_index][energy_index] <= (pitch[pitch_index] + 5)  ):
        dataplacer(time_index,pitch_index,energy_index,diffEflux_2,diffNflux_2,DiffE_list_2,DiffN_list_2)

    elif (0 < pitch_index < 20) and ( (pitch[pitch_index] - 5) > pitch_actual_2[time_index][pitch_index][energy_index] or  (pitch[pitch_index] + 5) < pitch_actual_2[time_index][pitch_index][energy_index]):
        for checker in range(21):
            if checker == 0 and (    (pitch[2] -5) <= pitch_actual_2[time_index][pitch_index][energy_index] <= (pitch[2] + 5)    ):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_2,diffNflux_2,DiffE_list_2,DiffN_list_2,checker)
            elif checker == 20 and (    (pitch[18] -5) <= pitch_actual_2[time_index][pitch_index][energy_index] <= (pitch[18] + 5)    ):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_2,diffNflux_2,DiffE_list_2,DiffN_list_2,checker)
            elif (0 < checker < 20) and (    (pitch[checker] -5) <= pitch_actual_2[time_index][pitch_index][energy_index] <= (pitch[checker] + 5)    ):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_2,diffNflux_2,DiffE_list_2,DiffN_list_2,checker)


    # FOR ANGLE -10
    if (pitch_index == 0) and ((pitch[2] - 5) <= pitch_actual_2[time_index][pitch_index][energy_index] <= (pitch[2] + 5)  ):
        dataplacer(time_index,pitch_index,energy_index,diffEflux_2,diffNflux_2,DiffE_list_2,DiffN_list_2)  # Place the data into the right spot

    elif (pitch_index == 0) and ( (pitch[2] - 5) > pitch_actual_2[time_index][pitch_index][energy_index] or  (pitch[2] + 5) < pitch_actual_2[time_index][pitch_index][energy_index]):
        for checker in range(21):
            if checker == 0 and ((pitch[2] - 5) <= pitch_actual_2[time_index][pitch_index][energy_index] <= (pitch[2] + 5)):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_2,diffNflux_2,DiffE_list_2,DiffN_list_2,checker)
            elif checker == 20 and ((pitch[18] - 5) <= pitch_actual_2[time_index][pitch_index][energy_index] <= (pitch[18] + 5)):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_2,diffNflux_2,DiffE_list_2,DiffN_list_2,checker)
            elif (0 < checker < 20) and ((pitch[checker] - 5) <= pitch_actual_2[time_index][pitch_index][energy_index] <= (pitch[checker] + 5)):
                # 2) Place the data into the right spot
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_2,diffNflux_2,DiffE_list_2,DiffN_list_2,checker)



    # FOR ANGLE 190
    if (pitch_index == 20) and ((pitch[18] - 5) <= pitch_actual_2[time_index][pitch_index][energy_index] <= (pitch[18] + 5)  ):
        dataplacer(time_index,pitch_index,energy_index,diffEflux_2,diffNflux_2,DiffE_list_2,DiffN_list_2)


    elif (pitch_index == 20) and ( (pitch[18] - 5) > pitch_actual_2[time_index][pitch_index][energy_index] or  (pitch[18] + 5) < pitch_actual_2[time_index][pitch_index][energy_index]):
        for checker in range(21):
            if checker == 0 and ((pitch[2] - 5) <= pitch_actual_2[time_index][pitch_index][energy_index] <= (pitch[2] + 5)):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_2,diffNflux_2,DiffE_list_2,DiffN_list_2,checker)
            elif checker == 20 and ((pitch[18] - 5) <= pitch_actual_2[time_index][pitch_index][energy_index] <= (pitch[18] + 5)):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_2,diffNflux_2,DiffE_list_2,DiffN_list_2,checker)
            elif (0 < checker < 20) and ((pitch[checker] - 5) <= pitch_actual_2[time_index][pitch_index][energy_index] <= (pitch[checker] + 5)):
                dataplacer_checker(time_index,pitch_index,energy_index,diffEflux_2,diffNflux_2,DiffE_list_2,DiffN_list_2,checker)

print('Done')


# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# -------- AVERAGING --------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
print('Averaging and filling in differential data for eepaa2: ', end='')

# Removed all data points that are not in lists and replace with fill values and AVERAGE the ones that are lists
for time_index,pitch_index,energy_index in itertools.product(*ranges_2_ex):
    # Remove
    if not isinstance(DiffE_list_2[time_index][pitch_index][energy_index],list):
        DiffE_list_2[time_index][pitch_index][energy_index] = fillval_2
    if not isinstance(DiffN_list_2[time_index][pitch_index][energy_index],list):
        DiffN_list_2[time_index][pitch_index][energy_index] = fillval_2
    # Average
    if isinstance(DiffE_list_2[time_index][pitch_index][energy_index], list):
        DiffE_list_2[time_index][pitch_index][energy_index] = sum(DiffE_list_2[time_index][pitch_index][energy_index])/len(DiffE_list_2[time_index][pitch_index][energy_index])
    if isinstance(DiffN_list_2[time_index][pitch_index][energy_index], list):
        DiffN_list_2[time_index][pitch_index][energy_index] = sum(DiffN_list_2[time_index][pitch_index][energy_index])/len(DiffN_list_2[time_index][pitch_index][energy_index])

print('Done')


print("--- %s seconds for program---" % (time.time() - flux_time2) ,'\n')




# -----------------------------------------
# Store the DiffE/N flux into it's own file
# -----------------------------------------

tag = 18

#-----DiffEFlux_aligned_1-----

distENAM = 'Differential_Energy_Flux_Aligned'
vardata= np.array(DiffE_list_1,dtype = 'float64')


#Var Info
varinfo = EEPAA_file_1.varinq(zvars_1[tag])
varinfo['Variable'] = distENAM

#Var Attributes
varattrs_high = EEPAA_file_1.varattsget(zvars_1[tag], expand=True)  # Get the variable's attributes
varattrs_high['CATDESC'] = distENAM
varattrs_high['FIELDNAM'] = distENAM
Sortedflux_high.write_var(varinfo,var_attrs=varattrs_high,var_data=vardata)


#-----DiffNFlux_aligned_1-----
distNNAM = 'Differential_Number_Flux_Aligned'
vardata= np.array(DiffN_list_1,dtype = 'float64')

#Var Info
varinfo = EEPAA_file_1.varinq(zvars_1[tag])
varinfo['Variable'] = distNNAM

#Var Attributes
varattrs_high = EEPAA_file_1.varattsget(zvars_1[tag], expand=True)  # Get the variable's attributes
varattrs_high['CATDESC'] = distNNAM
varattrs_high['FIELDNAM'] = distNNAM
Sortedflux_high.write_var(varinfo,var_attrs=varattrs_high,var_data=vardata)



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# -----DiffEFlux_aligned_2-----

distENAM = 'Differential_Energy_Flux_Aligned'
vardata= np.array(DiffE_list_2,dtype = 'float64')


#Var Info
varinfo = EEPAA_file_2.varinq(zvars_2[tag])
varinfo['Variable'] = distENAM

#Var Attributes
varattrs_low = EEPAA_file_2.varattsget(zvars_2[tag], expand=True)  # Get the variable's attributes
varattrs_low['CATDESC'] = distENAM
varattrs_low['FIELDNAM'] = distENAM
Sortedflux_low.write_var(varinfo,var_attrs=varattrs_low,var_data=vardata)


#-----DiffNFlux_Low-----
distNNAM = 'Differential_Number_Flux_Aligned'
vardata= np.array(DiffN_list_2,dtype = 'float64')


#Var Info
varinfo = EEPAA_file_2.varinq(zvars_2[tag])
varinfo['Variable'] = distNNAM

#Var Attributes
varattrs_low = EEPAA_file_2.varattsget(zvars_2[tag], expand=True)  # Get the variable's attributes
varattrs_low['CATDESC'] = distNNAM
varattrs_low['FIELDNAM'] = distNNAM
Sortedflux_low.write_var(varinfo,var_attrs=varattrs_low,var_data=vardata)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pitch_actual_file_1.close()
pitch_actual_file_2.close()
EEPAA_file_1.close()
EEPAA_file_2.close()
Sortedflux_high.close()
Sortedflux_low.close()


print("--- %s seconds for Sorting_DiffFlux---" % (time.time() - start_time) ,'\n')
