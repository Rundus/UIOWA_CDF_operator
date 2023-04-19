# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE Sorting the DiffEFlux Variable--------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


print('importing variables: ', end='')


import numpy as np, itertools,time
from cdflib import cdfwrite


start_time = time.time()




# ------------
# Import Files
# ------------

from files import EEPAA_file_high,EEPAA_file_low,pitch_file_actual_high,pitch_file_actual_low
from functions import dataplacerchecker_low,dataplacerchecker_high,dataplacer_low,dataplacer_high
from Variables import EEPAA_file_high_info,EEPAA_file_low_info,zvars_high,zvars_low,zvars_pitch_high,zvars_pitch_low,pitch_output_high as Pitch_output_high,pitch_output_low as Pitch_output_low
from Variables import DiffNFlux_low,DiffNFlux_high,DiffEFlux_low,DiffEFlux_high,EPOCH_low,EPOCH_High,pitch
from Variables import fillval_high,fillval_low,ranges_high,ranges_low,DiffE_list_low,DiffN_list_low,DiffE_list_high,DiffN_list_high




# ----------------------
# Output Files/Variables
# ----------------------
Sortedflux_high = cdfwrite.CDF('C:/Users/Rundus/Desktop/TRICE/DiffFlux_sorted_high',cdf_spec=EEPAA_file_high_info,delete=True)
Sortedflux_low = cdfwrite.CDF('C:/Users/Rundus/Desktop/TRICE/DiffFlux_sorted_low',cdf_spec=EEPAA_file_low_info,delete=True)


print('Done')


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ----------------------Sorting the DiffEData ---------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flux_time1 = time.time()
print('Reorganizing the Differential energy fluxes for high: ', end='')


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ----------- HIGH -----------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for time_index,pitch_index,energy_index in itertools.product(*ranges_high):
    # FOR ANGLES 0-180
    if (0 < pitch_index < 20) and ((pitch[pitch_index] - 5) <= Pitch_output_high[time_index][pitch_index][energy_index] <= (pitch[pitch_index] + 5)  ):
        dataplacer_high(time_index,pitch_index,energy_index)

    # If you find a pitch angle actual that's out of place, look through the other pitch angles at the time and energy. Find where it fits between +/-5 and put it there when it does. Special cases are given to check -10 and 190 because these need to be treated like 10 and 170
    elif (0 < pitch_index < 20) and ( (pitch[pitch_index] - 5) > Pitch_output_high[time_index][pitch_index][energy_index] or  (pitch[pitch_index] + 5) < Pitch_output_high[time_index][pitch_index][energy_index]):
        for checker in range(21):
            if checker == 0 and (    (pitch[2] -5) <= Pitch_output_high[time_index][pitch_index][energy_index] <= (pitch[2] + 5)    ):
                dataplacerchecker_high(time_index,pitch_index,energy_index,checker)
            elif checker == 20 and (    (pitch[18] -5) <= Pitch_output_high[time_index][pitch_index][energy_index] <= (pitch[18] + 5)    ):
                dataplacerchecker_high(time_index, pitch_index, energy_index, checker)
            elif (0 < checker < 20) and (    (pitch[checker] -5) <= Pitch_output_high[time_index][pitch_index][energy_index] <= (pitch[checker] + 5)    ):
                # 2) Place the data into the right spot
                dataplacerchecker_high(time_index,pitch_index,energy_index,checker)

    # FOR ANGLE -10
    if (pitch_index == 0) and ((pitch[2] - 5) <= Pitch_output_high[time_index][pitch_index][energy_index] <= (pitch[2] + 5)  ):
        dataplacer_high(time_index, pitch_index, energy_index)  # Place the data into the right spot

    elif (pitch_index == 0) and ( (pitch[2] - 5) > Pitch_output_high[time_index][pitch_index][energy_index] or  (pitch[2] + 5) < Pitch_output_high[time_index][pitch_index][energy_index]):
        for checker in range(21):
            if checker == 0 and ((pitch[2] - 5) <= Pitch_output_high[time_index][pitch_index][energy_index] <= (pitch[2] + 5)):
                dataplacerchecker_high(time_index, pitch_index, energy_index, checker)
            elif checker == 20 and ((pitch[18] - 5) <= Pitch_output_high[time_index][pitch_index][energy_index] <= (pitch[18] + 5)):
                dataplacerchecker_high(time_index, pitch_index, energy_index, checker)
            elif (0 < checker < 20) and ((pitch[checker] - 5) <= Pitch_output_high[time_index][pitch_index][energy_index] <= (pitch[checker] + 5)):
                # 2) Place the data into the right spot
                dataplacerchecker_high(time_index, pitch_index, energy_index, checker)

    # FOR ANGLE 190
    if (pitch_index == 20) and ((pitch[18] - 5) <= Pitch_output_high[time_index][pitch_index][energy_index] <= (pitch[18] + 5)  ):
        dataplacer_high(time_index, pitch_index, energy_index)

    elif (pitch_index == 20) and ( (pitch[18] - 5) > Pitch_output_high[time_index][pitch_index][energy_index] or  (pitch[18] + 5) < Pitch_output_high[time_index][pitch_index][energy_index]):
        for checker in range(21):
            if checker == 0 and ((pitch[2] - 5) <= Pitch_output_high[time_index][pitch_index][energy_index] <= (pitch[2] + 5)):
                dataplacerchecker_high(time_index, pitch_index, energy_index, checker)
            elif checker == 20 and ((pitch[18] - 5) <= Pitch_output_high[time_index][pitch_index][energy_index] <= (pitch[18] + 5)):
                dataplacerchecker_high(time_index, pitch_index, energy_index, checker)
            elif (0 < checker < 20) and ((pitch[checker] - 5) <= Pitch_output_high[time_index][pitch_index][energy_index] <= (pitch[checker] + 5)):
                dataplacerchecker_high(time_index, pitch_index, energy_index, checker)


print('Done')
print('Averaging and filling in differential data for high: ', end='')

# Removed all data points that are not in lists and replace with fill values and AVERAGE the ones that are lists
for time_index,pitch_index,energy_index in itertools.product(*ranges_high):
    # Remove
    if not isinstance(DiffE_list_high[time_index][pitch_index][energy_index],list):
        DiffE_list_high[time_index][pitch_index][energy_index] = fillval_high
    if not isinstance(DiffN_list_high[time_index][pitch_index][energy_index],list):
        DiffN_list_high[time_index][pitch_index][energy_index] = fillval_high
    # Average
    if isinstance(DiffE_list_high[time_index][pitch_index][energy_index], list):
        DiffE_list_high[time_index][pitch_index][energy_index] = sum(DiffE_list_high[time_index][pitch_index][energy_index])/len(DiffE_list_high[time_index][pitch_index][energy_index])
    if isinstance(DiffN_list_high[time_index][pitch_index][energy_index], list):
        DiffN_list_high[time_index][pitch_index][energy_index] = sum(DiffN_list_high[time_index][pitch_index][energy_index])/len(DiffN_list_high[time_index][pitch_index][energy_index])


print('Done')


print("--- %s seconds for sorting high---" % (time.time() - flux_time1) ,'\n')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# ----------- LOW -----------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
flux_time2 = time.time()
print('Reorganizing the diffEflux for low: ', end='')

for time_index,pitch_index,energy_index in itertools.product(*ranges_low):
    # FOR ANGLES 0-180
    if (0 < pitch_index < 20) and ((pitch[pitch_index] - 5) <= Pitch_output_low[time_index][pitch_index][energy_index] <= (pitch[pitch_index] + 5)  ):
        dataplacer_low(time_index,pitch_index,energy_index)

    elif (0 < pitch_index < 20) and ( (pitch[pitch_index] - 5) > Pitch_output_low[time_index][pitch_index][energy_index] or  (pitch[pitch_index] + 5) < Pitch_output_low[time_index][pitch_index][energy_index]):
        for checker in range(21):
            if checker == 0 and (    (pitch[2] -5) <= Pitch_output_low[time_index][pitch_index][energy_index] <= (pitch[2] + 5)    ):
                dataplacerchecker_low(time_index,pitch_index,energy_index,checker)
            elif checker == 20 and (    (pitch[18] -5) <= Pitch_output_low[time_index][pitch_index][energy_index] <= (pitch[18] + 5)    ):
                dataplacerchecker_low(time_index, pitch_index, energy_index, checker)
            elif (0 < checker < 20) and (    (pitch[checker] -5) <= Pitch_output_low[time_index][pitch_index][energy_index] <= (pitch[checker] + 5)    ):
                dataplacerchecker_low(time_index,pitch_index,energy_index,checker)


    # FOR ANGLE -10
    if (pitch_index == 0) and ((pitch[2] - 5) <= Pitch_output_low[time_index][pitch_index][energy_index] <= (pitch[2] + 5)  ):
        dataplacer_low(time_index, pitch_index, energy_index)
    elif (pitch_index == 0) and ( (pitch[2] - 5) > Pitch_output_low[time_index][pitch_index][energy_index] or  (pitch[2] + 5) < Pitch_output_low[time_index][pitch_index][energy_index]):
        for checker in range(21):
            if checker == 0 and ((pitch[2] - 5) <= Pitch_output_low[time_index][pitch_index][energy_index] <= (pitch[2] + 5)):
                dataplacerchecker_low(time_index, pitch_index, energy_index, checker)
            elif checker == 20 and ((pitch[18] - 5) <= Pitch_output_low[time_index][pitch_index][energy_index] <= (pitch[18] + 5)):
                dataplacerchecker_low(time_index, pitch_index, energy_index, checker)
            elif (0 < checker < 20) and ((pitch[checker] - 5) <= Pitch_output_low[time_index][pitch_index][energy_index] <= (pitch[checker] + 5)):
                dataplacerchecker_low(time_index, pitch_index, energy_index, checker)

    # FOR ANGLE 190
    if (pitch_index == 20) and ((pitch[18] - 5) <= Pitch_output_low[time_index][pitch_index][energy_index] <= (pitch[18] + 5)  ):
        dataplacer_low(time_index, pitch_index, energy_index)

    elif (pitch_index == 20) and ( (pitch[18] - 5) > Pitch_output_low[time_index][pitch_index][energy_index] or  (pitch[18] + 5) < Pitch_output_low[time_index][pitch_index][energy_index]):
        for checker in range(21):
            if checker == 0 and ((pitch[2] - 5) <= Pitch_output_low[time_index][pitch_index][energy_index] <= (pitch[2] + 5)):
                dataplacerchecker_low(time_index, pitch_index, energy_index, checker)
            elif checker == 20 and ((pitch[18] - 5) <= Pitch_output_low[time_index][pitch_index][energy_index] <= (pitch[18] + 5)):
                dataplacerchecker_low(time_index, pitch_index, energy_index, checker)
            elif (0 < checker < 20) and ((pitch[checker] - 5) <= Pitch_output_low[time_index][pitch_index][energy_index] <= (pitch[checker] + 5)):
                dataplacerchecker_low(time_index, pitch_index, energy_index, checker)

print('Done')



# %%%%%%%%%%%%%%%%%
# --- DIFFNFlUX ---
# %%%%%%%%%%%%%%%%%





# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# -------- AVERAGING --------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
print('Averaging and filling in differential data for low: ', end='')

# Removed all data points that are not in lists and replace with fill values and AVERAGE the ones that are lists
for time_index,pitch_index,energy_index in itertools.product(*ranges_low):
    # Remove
    if not isinstance(DiffE_list_low[time_index][pitch_index][energy_index],list):
        DiffE_list_low[time_index][pitch_index][energy_index] = fillval_low
    if not isinstance(DiffN_list_low[time_index][pitch_index][energy_index],list):
        DiffN_list_low[time_index][pitch_index][energy_index] = fillval_low
    # Average
    if isinstance(DiffE_list_low[time_index][pitch_index][energy_index], list):
        DiffE_list_low[time_index][pitch_index][energy_index] = sum(DiffE_list_low[time_index][pitch_index][energy_index])/len(DiffE_list_low[time_index][pitch_index][energy_index])
    if isinstance(DiffN_list_low[time_index][pitch_index][energy_index], list):
        DiffN_list_low[time_index][pitch_index][energy_index] = sum(DiffN_list_low[time_index][pitch_index][energy_index])/len(DiffN_list_low[time_index][pitch_index][energy_index])

print('Done')


print("--- %s seconds for program---" % (time.time() - flux_time2) ,'\n')




# -----------------------------------------
# Store the DiffE/N flux into it's own file
# -----------------------------------------

#-----DiffEFlux_aligned_High-----

distENAM = 'Differential_Energy_Flux_Aligned'
vardata= np.array(DiffE_list_high,dtype = 'float64')
tag = 6

#Var Info
varinfo = EEPAA_file_high.varinq(zvars_high[tag])
varinfo['Variable'] = distENAM

#Var Attributes
varattrs_high = EEPAA_file_high.varattsget(zvars_high[tag], expand=True)  # Get the variable's attributes
varattrs_high['CATDESC'] = distENAM
varattrs_high['FIELDNAM'] = distENAM
Sortedflux_high.write_var(varinfo,var_attrs=varattrs_high,var_data=vardata)


#-----DiffNFlux_aligned_High-----
distNNAM = 'Differential_Number_Flux_Aligned'
vardata= np.array(DiffN_list_high,dtype = 'float64')
tag = 13

#Var Info
varinfo = EEPAA_file_high.varinq(zvars_high[tag])
varinfo['Variable'] = distNNAM

#Var Attributes
varattrs_high = EEPAA_file_high.varattsget(zvars_high[tag], expand=True)  # Get the variable's attributes
varattrs_high['CATDESC'] = distNNAM
varattrs_high['FIELDNAM'] = distNNAM
Sortedflux_high.write_var(varinfo,var_attrs=varattrs_high,var_data=vardata)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#-----DiffEFlux_low-----

distENAM = 'Differential_Energy_Flux_Aligned'
vardata= np.array(DiffE_list_low,dtype = 'float64')
tag = 6

#Var Info
varinfo = EEPAA_file_low.varinq(zvars_low[tag])
varinfo['Variable'] = distENAM

#Var Attributes
varattrs_low = EEPAA_file_low.varattsget(zvars_low[tag], expand=True)  # Get the variable's attributes
varattrs_low['CATDESC'] = distENAM
varattrs_low['FIELDNAM'] = distENAM
Sortedflux_low.write_var(varinfo,var_attrs=varattrs_low,var_data=vardata)


#-----DiffNFlux_Low-----
distNNAM = 'Differential_Number_Flux_Aligned'
vardata= np.array(DiffN_list_low,dtype = 'float64')
tag = 13

#Var Info
varinfo = EEPAA_file_low.varinq(zvars_low[tag])
varinfo['Variable'] = distNNAM

#Var Attributes
varattrs_low = EEPAA_file_low.varattsget(zvars_low[tag], expand=True)  # Get the variable's attributes
varattrs_low['CATDESC'] = distNNAM
varattrs_low['FIELDNAM'] = distNNAM
Sortedflux_low.write_var(varinfo,var_attrs=varattrs_low,var_data=vardata)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pitch_actual_low_new.close()
pitch_actual_high_new.close()
EEPAA_file_high.close()
EEPAA_file_low.close()
Sortedflux_high.close()
Sortedflux_low.close()


print("--- %s seconds for Sorting_DiffFlux---" % (time.time() - start_time) ,'\n')
