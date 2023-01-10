# --- sun_glint_collection.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:
# [1] Reads in cdf file of attitude data and downsamples the roll angle to match the eepaa data.
# [2] Also produces a new variable counts_p1, where a simple subtraction mask is applied to all the data.
# [3] count data is collected on two regions (before and after) the auroral event that will be used to perform noise reduction statistics in the next code

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

MaskVal = [3,4] # Apply a mask over all the count values. Just a simple subtraction for every count

# --- --- --- --- ---
import time
from class_var_func import Done, setupPYCDF,prgMsg
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- import ---
# --- --- --- ---
import numpy as np
from tqdm import tqdm
from itertools import product
from os import remove,path
from data_paths import fliers,TRICE_data_folder,TRICE_L1_files
from class_var_func import cdf_attitude_TRICE,color,L1_TRICE,noise_region_before_event_in_counts,noise_region_after_event_in_counts
from missionAttributes import TRICE_mission_dicts
from copy import deepcopy
setupPYCDF()
from spacepy import pycdf

print(color.BOLD + color.CYAN + 'sun_glint_collection.py' + color.END + color.END)
def sun_glint_processing(wRocket,MaskVal):
    print(color.UNDERLINE +  f'{fliers[wRocket]} flier' + color.END)
    prgMsg('Acquiring data')

    # --- get the attitude data ---
    pycdf.lib.set_backward(False)

    inputFile = cdf_attitude_TRICE(wRocket)
    with pycdf.CDF(inputFile[0] ) as cdfFile:
        attitudeEpoch = np.array( [pycdf.lib.datetime_to_tt2000(cdfFile['Epoch'][i]) for i in range(len(cdfFile['Epoch'][...]))])
        rollAngle = np.array(cdfFile['roll'][...])

    inputFile =TRICE_L1_files
    with pycdf.CDF(inputFile[wRocket][0] ) as cdfFile:
        Epoch = np.array([ pycdf.lib.datetime_to_tt2000(cdfFile['Epoch'][i]) for i in range(len(cdfFile['Epoch'][...]))])
        counts = np.array(cdfFile['eepaa'][...])
        Energy = np.array(cdfFile['Energy'][...])
        Pitch_Angle = np.array(cdfFile['Pitch_Angle'][...])

    Done(start_time)

    rocketAttrs,missionDict,data_dict_temp = TRICE_mission_dicts()

    # --- Downsample Attitude Data ---
    prgMsg('Matching Epoch indices')

    att_Epoch_Downsampled = []
    roll_Downsampled = []

    nans = []
    for i in range(len(attitudeEpoch)):
        if np.isnan(rollAngle[i]):
            nans.append(i)


    for i in range(len(Epoch)):
        index = np.abs(attitudeEpoch - Epoch[i]).argmin()
        att_Epoch_Downsampled.append(attitudeEpoch[index])
        roll_Downsampled.append(float(rollAngle[index]))

    Done(start_time)


    # --- apply simple subraction mask over the counts data ---
    prgMsg('Subtracting count mask value')
    countsMasked = counts
    ranges = [range(len(counts)), range(len(counts[0])), range(len(counts[0][0]))]
    for tme, ptch, engy in tqdm(product(*ranges)):
        if counts[tme][ptch][engy] >= 65535:
            countsMasked[tme][ptch][engy] = rocketAttrs.fillVal
        elif (counts[tme][ptch][engy] - MaskVal) < 0:
            countsMasked[tme][ptch][engy] = 0
        else:
            countsMasked[tme][ptch][engy] = counts[tme][ptch][engy] - MaskVal

    Done(start_time)

    # collect information on the noise before and after the auroral event for later processing
    prgMsg('Collecting Noise Data')
    noise_data = []
    count_noise = [countsMasked[noise_region_before_event_in_counts[wRocket][0]:noise_region_before_event_in_counts[wRocket][1]],
                   countsMasked[noise_region_after_event_in_counts[wRocket][0]:noise_region_after_event_in_counts[wRocket][1]]]

    roll_noise = [roll_Downsampled[noise_region_before_event_in_counts[wRocket][0]:noise_region_before_event_in_counts[wRocket][1]],
                  roll_Downsampled[noise_region_after_event_in_counts[wRocket][0]:noise_region_after_event_in_counts[wRocket][1]]]

    for i in range(len(count_noise)):
        data = count_noise[i]
        roll = roll_noise[i]
        ranges = [range(len(data)), range(len(data[0])), range(len(data[0][0]))]
        for tme, ptch, engy in product(*ranges):
            noise_data.append([data[tme][ptch][engy], roll[tme], engy, ptch])

    Done(start_time)

    # --- --- --- --- --- --- --- ---
    # --- Prepare data for output ---
    # --- --- --- --- --- --- --- ---
    prgMsg('Writing out Data')

    # --- copy all data from model file to data dict ---
    data_dict = {}
    with pycdf.CDF(TRICE_L1_files[wRocket][0]) as L0DataFile:
        for key, val in L0DataFile.items():
            if key == 'Epoch':
                data_dict = {**data_dict, **{key: [Epoch, {key: val for key, val in L0DataFile[key].attrs.items()}]}}
            else:
                data_dict = {**data_dict, **{key: [ np.array(L0DataFile[key][...]), {key: val for key, val in L0DataFile[key].attrs.items()}]}}

    data_dict = {**data_dict, **{'counts_p1':    [np.array(countsMasked), deepcopy(data_dict_temp['data'][1])]}}
    data_dict = {**data_dict, **{'roll_downsampled': [np.array(roll_Downsampled), deepcopy(data_dict_temp['data'][1])]}}
    data_dict = {**data_dict, **{'Epoch_downsampled':[np.array(att_Epoch_Downsampled), deepcopy(data_dict_temp['Epoch'][1])]}}
    data_dict = {**data_dict, **{'noise_data': [np.array(noise_data), deepcopy(data_dict_temp['data'][1])]}}

    # --- Modify some attributes ---
    data_dict['counts_p1'][1]['VAR_TYPE'] = 'data'
    data_dict['counts_p1'][1]['DEPEND_1'] = 'Pitch_Angle'
    data_dict['counts_p1'][1]['DEPEND_2'] = 'Energy'



    # --- --- --- --- ---
    # --- OUTPUT Data ---
    # --- --- --- --- ---
    fileoutName = TRICE_L1_files[wRocket][0].replace('.cdf','').replace(TRICE_data_folder,'').replace(f'L1\{fliers[wRocket]}','')
    outputPath = f'{TRICE_data_folder}\\science\\{fliers[wRocket]}\\{fileoutName}_p1.cdf'

    # --- delete output file if it already exists ---
    if path.exists(outputPath):
        remove(outputPath)
    pycdf.lib.set_backward(False)

    # --- open the output file ---
    with pycdf.CDF(outputPath,'') as cdfFile:
        cdfFile.readonly(False)

        # --- WRITE OUT GLOBAL ATTRS ---
        globalAttrsMod = {'nothing':None}
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
                if varKey == 'Epoch' or varKey == 'Epoch_downsampled':
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
        if cdf_attitude_TRICE(i)[0] == []:
            print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
        else:
            sun_glint_processing(i,MaskVal[i])
else:
    if cdf_attitude_TRICE(wRocket[0])[0] == []:
        print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
    else:
        sun_glint_processing(wRocket[0],MaskVal[wRocket[0]])