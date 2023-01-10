# --- glint_noise_map.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:
# use data from sun_glint_collection.py's COUNT files to create a map that characterizes the sun glint noise
# For each pitch slice, a 2D map is created with Xaxis -> Roll angle and Yaxis -> energy
# assign the noise data obtained from sun_glint_collection to these slices and average the values
# Create a new counts variable where this mapping subtraction has been applied

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

# resolution of roll angles ( 1 --> 1 deg, 2 --> 2 deg, etc)
bin_spacing = 1




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
from data_paths import fliers,TRICE_data_folder, COUNTS_p1_files,TRICE_L1_files
from class_var_func import color,L1_TRICE
from missionAttributes import TRICE_mission_dicts
setupPYCDF()
from spacepy import pycdf

print(color.BOLD + color.CYAN + 'glint_noise_map.py' + color.END + color.END)
def glint_noise_map(wRocket):

    print(color.UNDERLINE + f'{fliers[wRocket]} flier' + color.END)
    prgMsg('Acquiring data')

    # --- get the COUNTS_p1 data ---
    pycdf.lib.set_backward(False)

    rocketAttrs,missionDict,data_dict_temp = TRICE_mission_dicts()
    data_dict = {}

    with pycdf.CDF(COUNTS_p1_files[wRocket][0]) as cdfFile:
        for key, val in cdfFile.items():
            data_dict = {**data_dict, **{key: [np.array(cdfFile[key][...]), {key: val for key, val in cdfFile[key].attrs.items()}]}}

    # --- --- --- --- --- --- --- ---
    # --- CREATE NOISE MAPPING ---
    # --- --- --- --- --- --- --- ---

    # Map will be slices in pitch which contain 2D arrays of X: roll angle, Y: Energy

    # Create bins for data to be sorted to. Xaxis is roll angle, yaxis is energy:
    roll_bins = np.linspace(-180, 180, int(360 / (bin_spacing)) + 1)
    noise_map = [ [ [ [] for k in range(len(data_dict['Energy'][0]))     ] for j in range(len(roll_bins))   ]  for i in range(len(data_dict['Pitch_Angle'][0]))]
    noise_map_stdDev = [ [ [ [] for k in range(len(data_dict['Energy'][0]))     ] for j in range(len(roll_bins))   ]  for i in range(len(data_dict['Pitch_Angle'][0]))]
    Done(start_time)

    prgMsg('Creating glint noise map')
    print('\n')

    # data format: countMasked, roll_Downsampled, Energy, Pitch_Angle

    for data in tqdm(data_dict['noise_data'][0]):
        ptch = int(data[3])
        engy = int(data[2])
        roll_index = np.abs(roll_bins - data[1]).argmin()
        noise_map[ptch][roll_index][engy].append(data[0])

    # --- Average/Adjust all the data ---

    prgMsg('Averaging and adjusting the mapping')
    ranges = [range(len(noise_map)), range(len(noise_map[0])), range(len(noise_map[0][0]))]

    for tme, ptch, engy in product(*ranges):
        if len(noise_map[tme][ptch][engy]) == 0:
            noise_map[tme][ptch][engy] = 0
        else:
            if wRocket == 0:
                if (ptch == 10 or ptch in range(14, 21)) and engy < 25:
                    mean_multiplier = 2
                    std_dev_multiplier = 2
                elif (ptch == 10 or ptch in range(14, 21)) and engy >= 25:
                    mean_multiplier = 1.6
                    std_dev_multiplier = 1.5
                elif ptch in range(0, 10):
                    mean_multiplier = 1.7
                    std_dev_multiplier = 1.7
                elif (ptch == 11) and engy < 25:
                    mean_multiplier = 2
                    std_dev_multiplier = 3
                elif (ptch == 12 or ptch == 13) and engy < 25:
                    mean_multiplier = 2.8
                    std_dev_multiplier = 3
                elif (ptch == 11 or ptch == 12 or ptch == 13) and engy >= 25:
                    mean_multiplier = 1.9
                    std_dev_multiplier = 1.5
                else:
                    mean_multiplier = 1
                    std_dev_multiplier = 1

            elif wRocket == 1:
                if ptch == 9:
                    mean_multiplier = 1.5
                    std_dev_multiplier = 2
                elif ptch == 11 and engy < 25:
                    mean_multiplier = 4.3
                    std_dev_multiplier = 2
                elif ptch == 11 and engy >= 25:
                    mean_multiplier = 3
                    std_dev_multiplier = 3
                elif ptch == 10 and engy < 25:
                    mean_multiplier = 4
                    std_dev_multiplier = 2
                elif ptch == 10 and engy >= 25:
                    mean_multiplier = 2.5
                    std_dev_multiplier = 2.5
                elif ptch == 12 and engy < 25:
                    mean_multiplier = 3.3
                    std_dev_multiplier = 3
                elif ptch == 12 and engy >= 25:
                    mean_multiplier = 3
                    std_dev_multiplier = 3
                else:
                    mean_multiplier = 1
                    std_dev_multiplier = 1

            noise_map[tme][ptch][engy] = int(mean_multiplier * np.mean(noise_map[tme][ptch][engy]))
            noise_map_stdDev[tme][ptch][engy] = std_dev_multiplier* np.std(noise_map[tme][ptch][engy])

    Done(start_time)

    # --- --- --- --- --- --- ---
    # --- APPLY MAP TO DATA ---
    # --- --- --- --- --- --- ---
    prgMsg('Applying noise map to data')
    print('\n')

    ranges = [range(len(data_dict['counts_p1'][0])),range(len(data_dict['counts_p1'][0][0])),range(len(data_dict['counts_p1'][0][0][0]))]
    counts_mapMasked = np.zeros(shape=(len(data_dict['counts_p1'][0]),len(data_dict['counts_p1'][0][0]),len(data_dict['counts_p1'][0][0][0])))

    for ptch in tqdm(ranges[1]):
        for tme in (ranges[0]):
            roll_index = np.abs(roll_bins - data_dict['roll_downsampled'][0][tme]).argmin()
            counts_mapMasked[tme][ptch] = [data_dict['counts_p1'][0][tme][ptch][engy] - (noise_map[ptch][roll_index][engy] + noise_map_stdDev[ptch][roll_index][engy]) if (data_dict['counts_p1'][0][tme][ptch][engy] - (noise_map[ptch][roll_index][engy] + noise_map_stdDev[ptch][roll_index][engy])  ) >= 0 else 0 for engy in (ranges[2])]


    # --- overwrite counts_p1 data ---
    data_dict['counts_p2'] = data_dict.pop('counts_p1')
    data_dict['counts_p2'][0] = counts_mapMasked

    # --- --- --- --- ---
    # --- OUTPUT Data ---
    # --- --- --- --- ---
    prgMsg('Writing out Data')
    fileoutName = TRICE_L1_files[wRocket][0].replace('.cdf', '').replace(TRICE_data_folder, '').replace(f'L1\{fliers[wRocket]}', '')
    outputPath = f'{TRICE_data_folder}\\science\\{fliers[wRocket]}\\{fileoutName}_p2.cdf'

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
        if path.exists(COUNTS_p1_files[i][0]):
            glint_noise_map(i)
        else:
            print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if path.exists(COUNTS_p1_files[wRocket[0]]):
        glint_noise_map(wRocket[0])
    else:
        print(color.RED + 'There are no .cdf files in the specified directory' + color.END)