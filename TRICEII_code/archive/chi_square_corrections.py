# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE Chi Square correction factors --------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ----
# IDEA
# ----

# Goal : There may be pads that are a bit too high or two low on the TRICE data set and must be corrected via chi-square minimization to fix this

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
from Variables import EPOCH_high_counts,EPOCH_low_counts,Energies
from functions import write_var_to_file
print('Done')

# ------------
# OUTPUT FILES
# ------------
cal_file_counts_high = cdfwrite.CDF(root + 'cal_counts_high',cdf_spec=counts_high_info,delete=True)
cal_file_counts_low = cdfwrite.CDF(root + 'cal_counts_low',cdf_spec=counts_low_info,delete=True)

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
thresh = 5
max1 = 1000000000000
min1 = 0
max2 = 1000000000000
min2 = 0


cor_pairs_high = [[2,1],[2,0],[0,3],[17,20],[20,16],[20,18],[20,19]]
cor_pairs_low = [[2,0],[0,1],[0,3],[0,4],[17,20],[20,16],[20,18],[18,19]]

for select in file_choice:
    print('\n\n-------------------------------------')
    print('--- Collecting pairs for file ' + names[select] +' ---\n', end='')
    print('-------------------------------------')

    # ------------------------------
    # GET INFO FOR THE SPECIFIC FILE
    # ------------------------------

    Which_roc = Titles[select]
    countvar = count_vars[select]
    pitchvar = pitch_vars[select]
    chi_range= chi_vars_ranges[select]
    rawdata = outputdat[select]
    output_data = rawdata
    pitch_choice = [[2], range(len(pitch))]
    degree = 1
    data = []


    # ---------------------
    # CORRECT THE DATA PADS
    # ---------------------

    # ---
    # LOW
    # ---
    if select == 0:
        # ---------------
        # -10deg --> 10deg
        # ---------------
        prin = cor_pairs_low[0][0]
        uncal = cor_pairs_low[0][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = countvar[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([countvar[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)
        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # ---------------
        # 0deg --> 10deg
        # ---------------
        prin = cor_pairs_low[1][0]
        uncal = cor_pairs_low[1][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = countvar[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([countvar[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # ---------------
        # 20deg --> -10deg
        # ---------------
        prin = cor_pairs_low[2][0]
        uncal = cor_pairs_low[2][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = output_data[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([output_data[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # -----------------
        # 30deg --> -10deg
        # -----------------
        prin = cor_pairs_low[3][0]
        uncal = cor_pairs_low[3][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = output_data[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([output_data[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # ---------------
        # 190deg --> 160deg
        # ---------------
        prin = cor_pairs_low[4][0]
        uncal = cor_pairs_low[4][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = countvar[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([countvar[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # -----------------
        # 150deg --> 190deg
        # -----------------
        prin = cor_pairs_low[5][0]
        uncal = cor_pairs_low[5][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = output_data[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([output_data[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # -----------------
        # 170deg --> 190deg
        # -----------------
        prin = cor_pairs_low[6][0]
        uncal = cor_pairs_low[6][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = output_data[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([output_data[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # -----------------
        # 180deg --> 170deg
        # -----------------
        prin = cor_pairs_low[7][0]
        uncal = cor_pairs_low[7][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = output_data[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([output_data[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

            # low
        vardata = output_data
        name = 'Calibrated_Counts'

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

        # -------------------------------------
        # Write out EPOCH,pitch,energy
        # -------------------------------------
        vardata = EPOCH_low_counts
        attributes = ['Epoch', 'ns', 'linear', EPOCH_low_counts.min(), EPOCH_low_counts.max(),counts_file_low.varattsget(zvars_counts_low[0], expand=True)]
        varinfo = counts_file_low.varinq(zvars_counts_low[0])
        write_var_to_file(cal_file_counts_low, varinfo, vardata, attributes)

        vardata = pitch
        attributes = ['Pitch_Angle', 'deg', 'linear', pitch.min(), pitch.max(),counts_file_low.varattsget(zvars_counts_low[9], expand=True)]
        varinfo = counts_file_low.varinq(zvars_counts_low[9])
        write_var_to_file(cal_file_counts_low, varinfo, vardata, attributes)

        vardata = Energies
        attributes = ['Energy', 'eV', 'log', Energies.min(), Energies.max(),counts_file_low.varattsget(zvars_counts_low[16], expand=True)]
        varinfo = counts_file_low.varinq(zvars_counts_low[16])
        write_var_to_file(cal_file_counts_low, varinfo, vardata, attributes)



    # ----
    # HIGH
    # ----
    elif select == 1:

        # ---------------
        # 0deg --> 10deg
        # ---------------
        prin = cor_pairs_high[0][0]
        uncal = cor_pairs_high[0][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = countvar[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([countvar[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # ---------------
        # -10deg --> 10deg
        # ---------------
        prin = cor_pairs_high[1][0]
        uncal = cor_pairs_high[1][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = countvar[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([countvar[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # ---------------
        # 20deg --> -10deg
        # ---------------
        prin = cor_pairs_high[2][0]
        uncal = cor_pairs_high[2][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = output_data[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([output_data[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # -----------------
        # 190deg --> 160deg
        # -----------------
        prin = cor_pairs_high[3][0]
        uncal = cor_pairs_high[3][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = countvar[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([countvar[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # ---------------
        # 150deg --> 190deg
        # ---------------
        prin = cor_pairs_high[4][0]
        uncal = cor_pairs_high[4][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = output_data[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([output_data[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # -----------------
        # 170deg --> 190deg
        # -----------------
        prin = cor_pairs_high[5][0]
        uncal = cor_pairs_high[5][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = output_data[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([output_data[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # -----------------
        # 180deg --> 190deg
        # -----------------
        prin = cor_pairs_high[6][0]
        uncal = cor_pairs_high[6][1]
        for tme in chi_range[0]:
            for engy in chi_range[1]:
                if abs(pitchvar[tme][prin][engy] - pitchvar[tme][uncal][engy]) <= 5:
                    check1 = output_data[tme][prin][engy] <= thresh
                    check2 = countvar[tme][uncal][engy] <= thresh
                    if not check1 and not check2:
                        data.append([output_data[tme][prin][engy], countvar[tme][uncal][engy]])

        uncalsdata = np.array([pair[1] for pair in data])
        prindata = np.array([pair[0] for pair in data])
        data = []
        slope, intercept = np.polyfit(uncalsdata, prindata, degree)

        print(pitch[prin],pitch[uncal],slope, intercept)

        # Apply Correction
        for tme in range(len(output_data)):
            for engy in range(49):
                tempdat = round(slope * output_data[tme][uncal][engy] + intercept)
                if tempdat < 0:
                    output_data[tme][uncal][engy] = 0
                else:
                    output_data[tme][uncal][engy] = tempdat

        # -------------------
        # WRITE OUT High Data
        # -------------------

        vardata = output_data
        name = 'Calibrated_Counts'
        # name = 'Count_diff'

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

        # -------------------------------------
        # Write out EPOCH,pitch,energy
        # -------------------------------------
        vardata = EPOCH_high_counts
        attributes = ['Epoch', 'ns', 'linear', EPOCH_high_counts.min(), EPOCH_high_counts.max(),counts_file_high.varattsget(zvars_counts_high[0], expand=True)]
        varinfo = counts_file_high.varinq(zvars_counts_high[0])
        write_var_to_file(cal_file_counts_high, varinfo, vardata, attributes)

        vardata = pitch
        attributes = ['Pitch_Angle', 'deg', 'linear', pitch.min(), pitch.max(),counts_file_high.varattsget(zvars_counts_high[9], expand=True)]
        varinfo = counts_file_high.varinq(zvars_counts_high[9])
        write_var_to_file(cal_file_counts_high, varinfo, vardata, attributes)

        vardata = Energies
        attributes = ['Energy', 'eV', 'log', Energies.min(), Energies.max(),counts_file_high.varattsget(zvars_counts_high[16], expand=True)]
        varinfo = counts_file_high.varinq(zvars_counts_high[16])
        write_var_to_file(cal_file_counts_high, varinfo, vardata, attributes)

print('Done')



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



print("--- %s seconds for Chi_square_corrections---" % (time.time() - start_time) ,'\n')















# --------------
    # PLOT THE DATA
    # --------------
    # plt.scatter(uncalsdata,prindata,color= 'b')
    # plt.scatter(polyfit_dat,prindata,color='r')
    # plt.scatter(der_data,prindata,color='g')
    # plt.xlabel('uncalibrated data [Counts]')
    # plt.ylabel('Principal Data [Counts]')
    # names = ['Raw Data','Polyfit Method: '+ 'm='  + str(slope) + ', b='+ str(intercept)  , 'Derivative Method: ' + r'$\alpha$=' + str(alpha_der)]
    # plt.legend( names )
    # plt.title(Which_roc + '\n' +'Principal Pad: ' + str(pitch[choice_pad[0]]) +'deg')
    # plt.show()