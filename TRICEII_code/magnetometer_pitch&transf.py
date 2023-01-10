# --- magnetometer_pitch&transf.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Use the TRICE layout configuration to transform the rocket magnetometer data into EEPAA frame
# and calculate the pitch angle from there

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
from copy import deepcopy
from math import cos,sin,pi,acos
from datetime import datetime as dt
from tqdm import tqdm
from itertools import product
from os import remove,path
from data_paths import fliers,TRICE_data_folder,TRICE_mag_files,TRICE_L1_files
from class_var_func import color,Mag_TRICE,take_closest
from missionAttributes import TRICE_mission_dicts
setupPYCDF()
from spacepy import pycdf

print(color.BOLD + color.CYAN + 'magnetometer_pitch&transf.py' + color.END + color.END)
def magnetometer_transformations(wRocket):

    print(color.UNDERLINE + f'{fliers[wRocket]} flier' + color.END)
    prgMsg('Acquiring data')

    # --- get the magnetometer data ---
    pycdf.lib.set_backward(False)

    rocketAttrs,missionDict,data_dict_temp = TRICE_mission_dicts()
    data_dict, data_dict_counts = {}, {}

    inputFile = TRICE_mag_files
    with pycdf.CDF(inputFile[wRocket][0]) as cdfFile:
        for key, val in cdfFile.items():
            data_dict = {**data_dict, **{key: [np.array(cdfFile[key][...]), {key: val for key, val in cdfFile[key].attrs.items()}]}}

    # --- get the eepaa data ---
    with pycdf.CDF(TRICE_L1_files[wRocket][0]) as cdfFile:
        for key, val in cdfFile.items():
            data_dict_counts = {**data_dict_counts, **{key: [np.array(cdfFile[key][...]), {key: val for key, val in cdfFile[key].attrs.items()}]}}

    eepaaEpoch_datetime = data_dict_counts['Epoch'][0]
    count_interval = data_dict_counts['Count_Interval'][0]
    Pitch_Angle = data_dict_counts['Pitch_Angle'][0]
    Energy = data_dict_counts['Energy'][0]
    eepaaEpoch = [pycdf.lib.datetime_to_tt2000(eepaaEpoch_datetime[i]) for i in range(len(eepaaEpoch_datetime))]
    ranges = [(range(len(eepaaEpoch))),(range(len(Pitch_Angle))),(range(len(Energy)))]


    # --- reorganize the mag data ---
    # NOTE: Mag time data is a timestamp from EPOCH (normally t=0 at 1/01/1970 00:00 UTC). I convert it to a TT2000 value.

    magEpoch = np.array([
                [pycdf.lib.datetime_to_tt2000(dt.utcfromtimestamp(data_dict['MagX'][0][i][0])) for i in range(len(data_dict['MagX'][0]))],
                [pycdf.lib.datetime_to_tt2000(dt.utcfromtimestamp(data_dict['MagY'][0][i][0])) for i in range(len(data_dict['MagY'][0]))],
                [pycdf.lib.datetime_to_tt2000(dt.utcfromtimestamp(data_dict['MagZ'][0][i][0])) for i in range(len(data_dict['MagZ'][0]))]
                ])

    magData = np.array([
                [data_dict['MagX'][0][i][1] for i in range(len(data_dict['MagX'][0]))],
                [data_dict['MagY'][0][i][1] for i in range(len(data_dict['MagY'][0]))],
                [data_dict['MagZ'][0][i][1] for i in range(len(data_dict['MagZ'][0]))]
                ])

    magDataTranspose = np.transpose(magData)
    Done(start_time)

    # --- downsample mag data to match EEPAA data at every epoch and energy value ---
    prgMsg('Downsampling Mag data')
    print('\n')
    magIndexMatrix = np.zeros(shape=(len(eepaaEpoch), 49))
    for i in tqdm(ranges[0]):
        magIndexMatrix[i] = [take_closest(magEpoch[0], (eepaaEpoch[i] + k * count_interval[i]) ) for k in ranges[2] ]


    # --- transform EEPAA frame into Mag Frame ---

    # X,Y,Z unit vectors for the EEPAA frame of reference
    unit_vect = np.zeros(shape=(21, 3))

    # create unit vectors representing EEPAA pitch sectors.
    # Rotate those vectors counterclockwise 90deg so that the EEPAA coordinates align with mag coordinates.
    for i in ranges[1]:
        unit_vect[i][0] = cos((Pitch_Angle[i] - 90) * pi / 180)
        unit_vect[i][1] = sin((Pitch_Angle[i] - 90) * pi / 180)
        unit_vect[i][2] = 0

    # --- Calculate the Pitch Angle from the magnetic field data using dot product ---
    pitch_angle_calc = np.ndarray(shape=(len(eepaaEpoch), len(Pitch_Angle), len(Energy)), dtype='float64')

    for tme, ptch in tqdm(product(ranges[0], ranges[1])):
        pitch_angle_calc[tme][ptch] = [(180 / pi) * acos((np.dot(unit_vect[ptch], magDataTranspose[int(magIndexMatrix[tme][engy])])) / (np.linalg.norm(magDataTranspose[int(magIndexMatrix[tme][engy])]))) for engy in range(len(Energy))]



    # --- --- --- --- --- --- --- ---
    # --- prepare data for output ---
    # --- --- --- --- --- --- --- ---
    del data_dict
    data_dict = {}

    data_dict = {**data_dict, **{'MagX': [np.array(magData[0]), deepcopy(data_dict_temp['data'][1])]}}
    data_dict = {**data_dict, **{'MagX_Epoch': [np.array(magEpoch[0]), deepcopy(data_dict_temp['Epoch'][1])]}}
    data_dict['MagX'][1]['DEPEND_0'] = 'MagX_Epoch'
    data_dict['MagX'][1]['VAR_TYPE'] = 'data'
    data_dict['MagX'][1]['UNITS'] = 'nanotesla'
    data_dict['MagX'][1]['LABLAXIS'] = 'MagX'

    data_dict = {**data_dict, **{'MagY': [np.array(magData[1]), deepcopy(data_dict_temp['data'][1])]}}
    data_dict = {**data_dict, **{'MagY_Epoch': [np.array(magEpoch[1]), deepcopy(data_dict_temp['Epoch'][1])]}}
    data_dict['MagY'][1]['DEPEND_0'] = 'MagY_Epoch'
    data_dict['MagY'][1]['VAR_TYPE'] = 'data'
    data_dict['MagY'][1]['UNITS'] = 'nanotesla'
    data_dict['MagY'][1]['LABLAXIS'] = 'MagY'

    data_dict = {**data_dict, **{'MagZ': [np.array(magData[2]), deepcopy(data_dict_temp['data'][1])]}}
    data_dict = {**data_dict, **{'MagZ_Epoch': [np.array(magEpoch[2]), deepcopy(data_dict_temp['Epoch'][1])]}}
    data_dict['MagZ'][1]['DEPEND_0'] = 'MagZ_Epoch'
    data_dict['MagZ'][1]['VAR_TYPE'] = 'data'
    data_dict['MagZ'][1]['UNITS'] = 'nanotesla'
    data_dict['MagZ'][1]['LABLAXIS'] = 'MagZ'

    data_dict = {**data_dict, **{'Calculated_Pitch_Angle': [np.array(pitch_angle_calc), deepcopy(data_dict_temp['data'][1])]}}
    data_dict['Calculated_Pitch_Angle'][1]['DEPEND_0'] = 'Epoch_EEPAA'
    data_dict['Calculated_Pitch_Angle'][1]['DEPEND_1'] = 'Pitch_Angle'
    data_dict['Calculated_Pitch_Angle'][1]['DEPEND_2'] = 'Energy'
    data_dict['Calculated_Pitch_Angle'][1]['VAR_TYPE'] = 'data'

    data_dict = {**data_dict, **{'Epoch_EEPAA': [np.array(eepaaEpoch), data_dict_counts['Epoch'][1]]}}
    data_dict = {**data_dict, **{'Pitch_Angle': [np.array(Pitch_Angle), data_dict_counts['Pitch_Angle'][1] ] }}
    data_dict = {**data_dict, **{'Energy': [np.array(Energy), data_dict_counts['Energy'][1]]}}

    # --- --- --- --- ---
    # --- OUTPUT Data ---
    # --- --- --- --- ---
    prgMsg('Writing out Data')
    fileoutName = f'TRICE_{rocketAttrs.rocketID[wRocket]}_pitch&mag'
    outputPath = f'{TRICE_data_folder}\\science\\{fliers[wRocket]}\\{fileoutName}_p3.cdf'

    # --- delete output file if it already exists ---
    if path.exists(outputPath):
        remove(outputPath)
    pycdf.lib.set_backward(False)

    # --- open the output file ---
    with pycdf.CDF(outputPath, '') as cdfFile:
        cdfFile.readonly(False)

        # --- WRITE OUT GLOBAL ATTRS ---
        globalAttrsMod = {'nothing': None}
        ModelData = Mag_TRICE(wRocket)[0]
        inputGlobDic = ModelData.cdfFile.globalattsget()

        for key, val in inputGlobDic.items():
            if key in globalAttrsMod:
                cdfFile.attrs[key] = globalAttrsMod[key]
            else:
                cdfFile.attrs[key] = val

        # --- WRITE OUT DATA ---
        for varKey, varVal in data_dict.items():
            if len(varVal[0]) != 0:
                if varKey in ['MagX_Epoch','MagY_Epoch','MagZ_Epoch','Epoch_EEPAA']:
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
        if path.exists(TRICE_mag_files[i][0]):
            magnetometer_transformations(i)
        else:
            print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if path.exists(TRICE_mag_files[wRocket[0]][0]):
        magnetometer_transformations(wRocket[0])
    else:
        print(color.RED + 'There are no .cdf files in the specified directory' + color.END)