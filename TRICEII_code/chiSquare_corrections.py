# --- chiSquare_corrections.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: get data from chiSquare_padPair.py and apply the
# correction factor obtained from chiSquare_analysis.py to calibrate TRICEII data post-flight


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

print(color.BOLD +color.CYAN +  'chiSquare_corrections.py'  +  color.END + color.END)
def chiSquare_corrections(wRocket):
    print(color.UNDERLINE + f'{fliers[wRocket]} flier' + color.END)
    prgMsg('Acquiring data')

    # --- rocket attribute data ---
    rocketAttrs,missionAttrs, data_dict_temp = TRICE_mission_dicts()
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
    ranges = [range(len(eepaaEpoch)), (range(len(Pitch_Angle))), (range(len(Energy)))]


    # --- --- --- --- --- --- --- --- ---
    # --- Apply ChiSquare calibration ---
    # --- --- --- --- --- --- --- --- ---







    # --- --- --- --- --- --- --- ---
    # --- Prepare Data for Output ---
    # --- --- --- --- --- --- --- ---
    data_dict = {}
    data_dict = {**data_dict, **{'chiSquare_padPairs': [np.array(padPairs), deepcopy(data_dict_temp['data'][1])]}}


    # --- --- --- --- ---
    # --- OUTPUT Data ---
    # --- --- --- --- ---

    prgMsg('Writing out Data')
    fileoutName = f'TRICE_{missionAttrs.rocketID[wRocket]}_padPairs'
    outputPath = f'{TRICE_data_folder}\\science\\{fileoutName}_p4.cdf'

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
        if path.exists(COUNTS_p1_files[i]):
            chiSquare_corrections(i)
        else:
            print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if path.exists(COUNTS_p1_files[wRocket[0]]):
        chiSquare_corrections(wRocket[0])
    else:
        print(color.RED + 'There are no .cdf files in the specified directory' + color.END)