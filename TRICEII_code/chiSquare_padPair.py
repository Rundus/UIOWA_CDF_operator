# --- chiSquare_padPair.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: A code to collect data and perform post-flight calibration of TRICE II data.

# Method:
# [1] look for a pair of detector pads that overlay in their actual pitch angle, namely pitch angle with fall within 5deg of one paid (e.g. pad 40deg is at 47deg true, and pad 50deg at 49 deg true)
#       -> two pads looking at approximately the same true pitch angle should have about the same flux value, if they don't its due to the pad calibration not the physics
# [2] For any one time slice, sweep through in energy the two pads of interest and store necesary data (counts, true pitch angle, and index identifiers).
# [4] The data is now paired; one count value & pitch for one pad also one count value and pitch for the second pad. One of these pads is the principal look angles, the other is the uncal pad.
# [5] Continue process [4] through all time slices and collect all the pairs between the two pads
# [6] You will now preform the chi-square minimization technique given by \chi ~ 1/(N-1) \sum \frac{data_{principal} - \alpha * data_{uncal}   }{ error_{principal} + error_{uncal}}
# [7] Create a method to sweep through possible values for the alpha fitting parameter such that \chi becomes ~1
# [8] Repeat steps [1] - [7] for the next pad pairs
# [9] Likely I will move from low pads to high pads, do the results become the same if I perform them in reverse?



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
wRocket = 0

# 0 --> High Flier
# 1 --> Low Flier

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
wRocket = []
# [] --> Both Fliers
# [0] --> High Flier
# [1] --> Low Flier

value_threshold = 10 # BOTH values must be larger than this to be stored
outlier_threshold = 10 # both values must have a ratio less than this and greater than 1/this to be stored





# --- --- --- --- ---
import time
from class_var_func import Done, setupPYCDF,prgMsg
start_time = time.time()
# --- --- --- --- ---





# --- --- --- ---
# --- import ---
# --- --- --- ---
import numpy as np
from copy import deepcopy
from tqdm import tqdm
from itertools import product
from os import remove,path
from data_paths import fliers,TRICE_data_folder, COUNTS_p1_files,CALC_pitch_angle_p3_files,COUNTS_p2_files
from class_var_func import color,L1_TRICE
from missionAttributes import TRICE_mission_dicts
setupPYCDF()
from spacepy import pycdf

print(color.BOLD +color.CYAN +  'chiSquare_padPair.py'  +  color.END + color.END)
def chiSquare_padPair(wRocket):
    print(color.UNDERLINE + f'{fliers[wRocket]} flier' + color.END)
    prgMsg('Acquiring data')

    # --- rocket attribute data ---
    rocketAttrs,missionAttributes,data_dict_temp = TRICE_mission_dicts()
    pycdf.lib.set_backward(False)

    # --- get the Calc_Pitch/Mag data ---
    data_dict_pitch = {}
    inputFile = CALC_pitch_angle_p3_files
    with pycdf.CDF(inputFile[wRocket]) as cdfFile:
        for key, val in cdfFile.items():
            data_dict_pitch = {**data_dict_pitch, **{key: [np.array(cdfFile[key][...]), {key: val for key, val in cdfFile[key].attrs.items()}]}}

    # --- get the masked eepaa data ---
    data_dict_counts = {}
    inputFile = COUNTS_p2_files
    with pycdf.CDF(inputFile[wRocket]) as cdfFile:
        for key, val in cdfFile.items():
            data_dict_counts = {**data_dict_counts, **{key: [np.array(cdfFile[key][...]), {key: val for key, val in cdfFile[key].attrs.items()}]}}

    eepaaEpoch_datetime = data_dict_counts['Epoch'][0]
    eepaaEpoch = [pycdf.lib.datetime_to_tt2000(eepaaEpoch_datetime[i]) for i in range(len(eepaaEpoch_datetime))]
    Pitch_Angle = data_dict_counts['Pitch_Angle'][0]
    Energy = data_dict_counts['Energy'][0]
    counts_p2 = data_dict_counts['counts_p2'][0]

    # --- --- --- --- --- --- --- --- ---
    # --- Collect calibration pairs ---
    # --- --- --- --- --- --- --- --- ---
    padPairs = []
    # append format:
    # [0] pad_principle,
    # [1] counts_principle,
    # [2] pad_uncal,
    # [3] counts_uncal,
    # [4] time_index
    # [5] energy_index


    # Look through every pad and find any close pitch pairs
    spread = 6
    for pad in tqdm(range(len(Pitch_Angle))):
        padPitch = Pitch_Angle[pad]
        checkPtches = [num for num in range(pad - spread, pad + spread + 1) if (num >= 0 and num <= 20)]
        ranges = [range(len(eepaaEpoch)), range(min(checkPtches),max(checkPtches)+1), range(len(Energy))]


        for tme, ptch, engy in product(*ranges):

            if ptch != pad: # make sure pad doesn't check itself

                if pad != 0 and pad != 20: # Special case to make sure -10,190 don't check the 10deg and 170deg bins
                    diff = abs(padPitch - data_dict_pitch['Calculated_Pitch_Angle'][0][tme][ptch][engy])

                    if diff <= 5 and data_dict_pitch['Calculated_Pitch_Angle'][0][tme][ptch][engy] != -1:
                        if counts_p2[tme][pad][engy] > value_threshold and counts_p2[tme][ptch][engy] > value_threshold:
                            if counts_p2[tme][pad][engy]/counts_p2[tme][ptch][engy] <= outlier_threshold and counts_p2[tme][pad][engy]/counts_p2[tme][ptch][engy] >= 1/outlier_threshold:
                                padPairs.append([pad, counts_p2[tme][pad][engy], ptch, counts_p2[tme][ptch][engy], tme, engy])

                elif pad == 0 and ptch != 2:
                    diff = abs(padPitch - data_dict_pitch['Calculated_Pitch_Angle'][0][tme][ptch][engy])

                    if diff <= 5 and data_dict_pitch['Calculated_Pitch_Angle'][0][tme][ptch][engy] != -1:
                        if counts_p2[tme][pad][engy] > value_threshold and counts_p2[tme][ptch][engy] > value_threshold:
                            if counts_p2[tme][pad][engy] / counts_p2[tme][ptch][engy] <= outlier_threshold and counts_p2[tme][pad][engy] / counts_p2[tme][ptch][engy] >= 1/outlier_threshold:
                                padPairs.append([pad, counts_p2[tme][pad][engy], ptch, counts_p2[tme][ptch][engy], tme, engy])

                elif pad == 20 and ptch != 18:
                    diff = abs(padPitch - data_dict_pitch['Calculated_Pitch_Angle'][0][tme][ptch][engy])

                    if diff <= 5 and data_dict_pitch['Calculated_Pitch_Angle'][0][tme][ptch][engy] != -1:
                        if counts_p2[tme][pad][engy] > value_threshold and counts_p2[tme][ptch][engy] > value_threshold:
                            if counts_p2[tme][pad][engy] / counts_p2[tme][ptch][engy] <= outlier_threshold and counts_p2[tme][pad][engy] / counts_p2[tme][ptch][engy] >= 1/outlier_threshold:
                                padPairs.append([pad, counts_p2[tme][pad][engy], ptch, counts_p2[tme][ptch][engy], tme, engy])

    # --- --- --- --- --- --- --- ---
    # --- Prepare Data for Output ---
    # --- --- --- --- --- --- --- ---
    data_dict = {}
    data_dict = {**data_dict, **{'chiSquare_padPairs': [np.array(padPairs), deepcopy(data_dict_temp['data'][1])]}}

    # --- --- --- --- ---
    # --- OUTPUT Data ---
    # --- --- --- --- ---

    prgMsg('Writing out Data')
    fileoutName = f'TRICE_{rocketAttrs.rocketID[wRocket]}_padPairs'
    outputPath = f'{TRICE_data_folder}\\science\\{fliers[wRocket]}\\{fileoutName}_p4.cdf'

    # --- delete output file if it already exists ---
    if path.exists(outputPath):
        remove(outputPath)
    pycdf.lib.set_backward(False)

    # --- open the output file ---
    with pycdf.CDF(outputPath, '') as cdfFile:
        cdfFile.readonly(False)

        # --- WRITE OUT GLOBAL ATTRS ---
        globalAttrsMod = {'nothing': None}
        ModelData = L1_TRICE(wRocket)[0]
        inputGlobDic = ModelData.cdfFile.globalattsget()

        for key, val in inputGlobDic.items():
            if key in globalAttrsMod:
                cdfFile.attrs[key] = globalAttrsMod[key]
            else:
                cdfFile.attrs[key] = val

        # --- WRITE OUT DATA ---
        for varKey, varVal in data_dict.items():
            if len(varVal[0]) != 0:
                if varKey == 'Epoch':
                    cdfFile.new(varKey, data=varVal[0], type=33)
                else:
                    cdfFile.new(varKey, data=varVal[0])

                # --- Write out the attributes and variable info ---
                for attrKey, attrVal in data_dict[varKey][1].items():

                    if attrKey == 'VALIDMIN':
                        cdfFile[varKey].attrs[attrKey] = varVal[0].min()
                    elif attrKey == 'VALIDMAX':
                        cdfFile[varKey].attrs[attrKey] = varVal[0].max()
                    elif attrVal != None:
                        cdfFile[varKey].attrs[attrKey] = attrVal
    Done(start_time)
    print('\n')

# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

if wRocket == []:
    for i in range(2):
        if path.exists(CALC_pitch_angle_p3_files[i]):
            chiSquare_padPair(i)
        else:
            print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if path.exists(CALC_pitch_angle_p3_files[wRocket[0]]):
        chiSquare_padPair(wRocket[0])
    else:
        print(color.RED + 'There are no .cdf files in the specified directory' + color.END)