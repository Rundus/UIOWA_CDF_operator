# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE Chi Square correction factors --------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ----
# IDEA
# ----

# Goal : There may be pads that are a bit too high or too low on the TRICE data set and must be corrected via chi-square minimization to fix this

# Method:
# [1] We first look for a pair of detector pads that overlay in their actual pitch angle, namely pitch angle with fall within 5deg of one paid (e.g. pad 40deg is at 47deg true, and pad 50deg at 49 deg true)
#       -> The logic here is that two pads looking at approximately the same true pitch angle should have about the same flux value, if they don't its due to the pad calibration not the physics
# [2] For any one time slice, sweep through in energy the two pads of interest and store necesary data (counts, true pitch angle, and index identifiers).
# [4] The data is now paired; one count value & pitch for one pad also one count value and pitch for the second pad. One of these pads is the principal look angles, the other is the uncal pad.
# [5] Continue process [4] through all time slices and collect all the pairs between the two pads
# [6] You will now preform the chi-square minimization technique given by \chi ~ 1/(N-1) \sum \frac{data_{principal} - \alpha * data_{uncal}   }{ error_{principal} + error_{uncal}}
# [7] Create a method to sweep through possible values for the alpha fitting parameter such that \chi becomes ~1
# [8] Repeat steps [1] - [7] for the next pad pairs
# [9] Likely I will move from low pads to high pads, do the results become teh same if I perform them in reverse?


import itertools
import time
print('importing variables: ', end='')
start_time = time.time()

# -----------
# IMPORT INFO
# -----------
import matplotlib.pyplot as plt
import numpy as np,math
from files import counts_file_high,counts_file_low,root
from cdflib import cdfwrite
from Variables import counts_low,counts_high,pitch_output_high,pitch_output_low,chi_ranges_high,chi_ranges_low,pitch,zvars_counts_high,zvars_counts_low,counts_high_info,counts_low_info
from Variables import EPOCH_High,EPOCH_low,ranges_high,ranges_low,rem_sun_spike_counts_low,rem_sun_spike_counts_high
print('Done')

# ------------
# OUTPUT FILES
# ------------
cal_file_counts_high = cdfwrite.CDF(root + 'cal_counts_high',cdf_spec=counts_high_info,delete=True)
cal_file_counts_low = cdfwrite.CDF(root + 'chi',cdf_spec=counts_low_info,delete=True)

# ----------------
# PRELIMINARY DATA
# ----------------
file_choice = [0,1]
names = ['LOW','HIGH']
checker_ranges = [ranges_low,ranges_high]
Titles = ['TRICE LOW','TRICE HIGH']
output_cal_counts_high = rem_sun_spike_counts_high
output_cal_counts_low = rem_sun_spike_counts_low
count_vars = [counts_low,counts_high]
outputdat = [output_cal_counts_low,output_cal_counts_high]
pitch_vars = [pitch_output_low,pitch_output_high]
chi_vars_ranges = [chi_ranges_low,chi_ranges_high] # Checks all the energies from 2979eV to 60eV (we don't need to check the higher energies because the counts are basically nothing there)
epochs = [range(len(counts_low)),range(len(counts_high))]

#-------=
# TOGGLES
#--------
var_name ='Calibrated_Counts'
degree = 1 # The degree of fitting equation for determining the correction
thresh = 0 # Determines the number that data values must exceed to be counted in the analysis

for select in file_choice:
    print('Collecting pairs for file ' + names[select] +':', end='')
    print('\n')

    # ------------------------------
    # GET INFO FOR THE SPECIFIC FILE
    # ------------------------------

    Which_roc = Titles[select]
    countvar = count_vars[select]
    pitchvar = pitch_vars[select]
    chi_range= chi_vars_ranges[select]
    output = outputdat[select]
    output_data = output

    # --------------------------
    # AQUIRE AND SORT DATA PAIRS
    # --------------------------

    dat_var = [[] for i in range(21)]
    corrections = [[] for i in range(21)]
    averagess = []
    uplift = [['-10deg',[]],['0deg',[]],['10deg',[]],['20deg',[]],['30deg',[]],['40deg',[]],['50deg',[]],['60deg',[]],['70deg',[]],['80deg',[]],['90deg',[]],['100deg',[]],['110deg',[]],['120deg',[]],['130deg',[]],['140deg',[]],['150deg',[]],['160deg',[]],['170deg',[]],['180deg',[]],['190deg',[]]]
    pitch_choice = [[2],range(len(pitch))]


    for prin in pitch_choice[1]:

        pair_no = [0 for k in range(21)]
        pair_dat = [[] for j in range(21)]
        out_dat = [[] for i in range(21)]

        for uncal in range(len(pitch)):

            if prin != uncal:

                for tme in chi_range[0]:

                    for engy in chi_range[1]:

                        if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                            check1 = countvar[tme][prin][engy] <= thresh
                            check2 = countvar[tme][uncal][engy] <= thresh

                            if not check1 and not check2: #Makes sure neither of the values are zero

                                dat_var[prin].append(   [countvar[tme][prin][engy],countvar[tme][uncal][engy]]   )
                                pair_no[uncal] += 1
                                pair_dat[uncal].append([countvar[tme][prin][engy],countvar[tme][uncal][engy]])

                                if uncal not in uplift[prin][1]:
                                    uplift[prin][1].append(uncal)

        uplift[prin][1].append(pair_no)
        uplift[prin].append(pair_dat)

        # Store the polyfitted parameters in out_dat and store out_dat in the appropiate uplift bin
        counter = -1
        for dataset in pair_dat:
            counter += 1
            if len(dataset) != 0:
                uncalsdata = np.array([pair[1] for pair in dataset])
                prindata = np.array([pair[0] for pair in dataset])
                slope, intercept = np.polyfit(uncalsdata, prindata, degree)
                out_dat[counter].append([slope,intercept])

        uplift.append(out_dat)
        temp_data = []
        for i in range(len(out_dat)):
            if out_dat[i] != []:
                temp_data.append([pitch[i],pair_no[i],out_dat[i]])

        if temp_data != []:
            print("Principal Pad:", str(pitch[prin]))
        for correction in temp_data:
            print(correction)
        if temp_data != []:
            print('\n')

    print('Done')


    # ---------------------------------
    # ####OUTPUT UPDATED COUNT DATA####
    # ---------------------------------
    max1 = 1000000000000
    min1 = 0
    max2 = 1000000000000
    min2 = 0


    if select == 0:
        # low
        vardata = output_data
        name = var_name

        varinfo = counts_file_low.varinq(zvars_counts_low[18])
        varinfo['Variable'] = name

        varattrs_low = counts_file_low.varattsget(zvars_counts_low[18], expand=True)
        varattrs_low['CATDESC'] = name
        varattrs_low['FIELDNAM'] = name
        varattrs_low['UNITS'] = 'Counts'
        varattrs_low['SCALETYP'] = 'log'
        varattrs_low['VALIDMIN'] = [min2]
        varattrs_low['VALIDMAX'] = [max2]
        varattrs_low['LABLAXIS'] = name
        varattrs_low['LABL_PTR_1'] = 'LABL_PTR_1'
        varattrs_low['LABL_PTR_2'] = 'LABL_PTR_2'

        cal_file_counts_low.write_var(varinfo, var_attrs=varattrs_low, var_data=vardata)


    elif select == 1:

        # High
        vardata = output_data
        name = var_name

        varinfo = counts_file_high.varinq(zvars_counts_high[18])
        varinfo['Variable'] = name

        varattrs_high = counts_file_high.varattsget(zvars_counts_high[18], expand=True)
        varattrs_high['CATDESC'] = name
        varattrs_high['FIELDNAM'] = name
        varattrs_high['UNITS'] = 'Counts'
        varattrs_high['SCALETYP'] = 'log'
        varattrs_high['VALIDMIN'] = [min1]
        varattrs_high['VALIDMAX'] = [max1]
        varattrs_high['LABLAXIS'] = name
        varattrs_high['LABL_PTR_1'] = 'LABL_PTR_1'
        varattrs_high['LABL_PTR_2'] = 'LABL_PTR_2'

        cal_file_counts_high.write_var(varinfo, var_attrs=varattrs_high, var_data=vardata)




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



print("--- %s seconds for Chi_square_corrections---" % (time.time() - start_time) ,'\n')