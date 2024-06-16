# --- deSpin_practice.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: play around with the TRICE mag&ACS data to learn how to despin


# There are multiple steps to process the mag data:
# [1] Transform mag into payload frame coordinates
# [2] Determine the gain and offsets of the magnetometer
# [3] Apply a linear interpolation to up-sample the ACS data to the mag data
# [4] Apply the attitude solution matrix to convert to ENU coordinates --> scipy SLURP
# [5]



import matplotlib.pyplot as plt

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
wRocket = 0
wComponent = 0

# 0 --> High Flier
# 1 --> Low Flier

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
wRocket = [0]
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
from scipy.spatial.transform import Slerp
from scipy.spatial.transform import Rotation as R
from copy import deepcopy
from datetime import datetime as dt
from os import remove, path
from data_paths import fliers, TRICE_data_folder, TRICE_mag_files, TRICE_ACS_files
from class_var_func import color, Mag_TRICE
from missionAttributes import mission_dicts
setupPYCDF()
from spacepy import pycdf


def euler_to_quaternion(yaw, pitch, roll):
    qx = np.sin(roll / 2) * np.cos(pitch / 2) * np.cos(yaw / 2) - np.cos(roll / 2) * np.sin(pitch / 2) * np.sin(yaw / 2)
    qy = np.cos(roll / 2) * np.sin(pitch / 2) * np.cos(yaw / 2) + np.sin(roll / 2) * np.cos(pitch / 2) * np.sin(yaw / 2)
    qz = np.cos(roll / 2) * np.cos(pitch / 2) * np.sin(yaw / 2) - np.sin(roll / 2) * np.sin(pitch / 2) * np.cos(yaw / 2)
    qw = np.cos(roll / 2) * np.cos(pitch / 2) * np.cos(yaw / 2) + np.sin(roll / 2) * np.sin(pitch / 2) * np.sin(yaw / 2)

    return [qx, qy, qz, qw]

def quaternion_to_euler(x, y, z, w):
    import math
    t0 = +2.0 * (w * x + y * z)
    t1 = +1.0 - 2.0 * (x * x + y * y)
    X = math.degrees(math.atan2(t0, t1))

    t2 = +2.0 * (w * y - z * x)
    t2 = +1.0 if t2 > +1.0 else t2
    t2 = -1.0 if t2 < -1.0 else t2
    Y = math.degrees(math.asin(t2))

    t3 = +2.0 * (w * z + x * y)
    t4 = +1.0 - 2.0 * (y * y + z * z)
    Z = math.degrees(math.atan2(t3, t4))

    return X, Y, Z



print(color.BOLD + color.CYAN + 'deSpin.py' + color.END + color.END)

def deSpin(wRocket,wComponent):
    print(color.UNDERLINE + f'{fliers[wRocket]} flier' + color.END)

    prgMsg('Acquiring data')
    pycdf.lib.set_backward(False)

    missionAttrs, data_dict_temp = mission_dicts()
    data_dict, data_dict_counts = {}, {}

    # --- --- --- --- --- --- --- --- --- --- --- ---
    # --- Step 0: get the magnetometer/ACS data ---
    # --- --- --- --- --- --- --- --- --- --- --- ---

    # --- Mag data ---
    with pycdf.CDF(TRICE_mag_files[wRocket][0]) as cdfFile:
        for key, val in cdfFile.items():
            data_dict = {**data_dict, **{key: [np.array(cdfFile[key][...]), {key: val for key, val in cdfFile[key].attrs.items()}]}}

    # --- ACS data ---
    with pycdf.CDF(TRICE_ACS_files[wRocket][0]) as cdfFile:
        for key, val in cdfFile.items():
            data_dict = {**data_dict, **{key: [np.array(cdfFile[key][...]), {key: val for key, val in cdfFile[key].attrs.items()}]}}

    magEpoch = np.array([
                [pycdf.lib.datetime_to_tt2000(dt.utcfromtimestamp(data_dict['MagX'][0][i][0])) for i in range(len(data_dict['MagX'][0]))],
                [pycdf.lib.datetime_to_tt2000(dt.utcfromtimestamp(data_dict['MagY'][0][i][0])) for i in range(len(data_dict['MagY'][0]))],
                [pycdf.lib.datetime_to_tt2000(dt.utcfromtimestamp(data_dict['MagZ'][0][i][0])) for i in range(len(data_dict['MagZ'][0]))]
                ])

    magData = np.array([
        [data_dict['MagX'][0][i][1],
         data_dict['MagY'][0][i][1],
         data_dict['MagZ'][0][i][1]]
        for i in range(len(data_dict['MagZ'][0]))
    ])

    magDataT = np.transpose(magData)

    # Find the first point that ACS and mag data start roughly together
    ACStime = [pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][i]) for i in range(len(data_dict['Epoch'][0]))]
    nanIndex = np.where(np.isnan(data_dict['a11'][0]))  # Find where data_dict == nan
    startIndex = nanIndex[0][-1] + 1
    startTime = ACStime[startIndex]

    # Get ACS variables
    ACSroll = data_dict['roll'][0]
    ACSpitch = data_dict['pitch'][0]
    ACSyaw = data_dict['yaw'][0]




    ENUmatrix = [
        [
            [data_dict['a11'][0][i], data_dict['a12'][0][i], data_dict['a13'][0][i]],
            [data_dict['a21'][0][i], data_dict['a22'][0][i], data_dict['a23'][0][i]],
            [data_dict['a31'][0][i], data_dict['a32'][0][i], data_dict['a33'][0][i]]
        ]
        for i in range(len(data_dict['a11'][0]))
    ]

    # --- --- --- --- --- --- --- --- --- --- ---
    # --- STEP 1: Rotating into rocket frame ---
    # --- --- --- --- --- --- --- --- --- --- ---

    # Mag requires the following to get to rocket frame:
    # (1) 45deg clockwise around y-axis
    r = R.from_rotvec([0, - np.pi / 4, 0])
    magDataYm45 = r.apply(magData)
    # (2) 90deg clockwise around x-axis
    r = R.from_rotvec([- np.pi / 2, 0, 0])
    magDataXm90 = r.apply(magDataYm45)
    # (3) 90deg clockwise around y-axis
    r = R.from_rotvec([0, - np.pi / 2, 0])
    magDataRF =np.transpose( r.apply(magDataXm90) )


    # --- --- --- --- --- --- --- --- --- --- ---
    # --- STEP 2: Determine MAG Gain/Offset ---
    # --- --- --- --- --- --- --- --- --- --- ---

    # --- --- --- --- --- --- --- --- --- --- --- ---
    # --- STEP 3: Convert ACS data to Quaternions ---
    # --- --- --- --- --- --- --- --- --- --- --- ---

    ACSQuaternions = [euler_to_quaternion((np.pi/180)*ACSyaw[index], (np.pi/180)*ACSroll[index], (np.pi/180)*ACSpitch[index]) for index in range(len(ACStime[startIndex:]))]


    # --- --- --- --- --- --- --- --- --- ---
    # --- STEP 4: interpolate ACS Data ---
    # --- --- --- --- --- --- --- --- --- ---




    # Adjust the mag data startpoint to where ACS starts
    startX = np.abs(magEpoch[0] - startTime).argmin()
    startY = np.abs(magEpoch[1] - startTime).argmin()
    startZ = np.abs(magEpoch[2] - startTime).argmin()

    magDataACSstart =[
        magDataRF[0][startX:],
        magDataRF[1][startY:],
        magDataRF[2][startZ:]
    ]

    magEpochACSstart = [
        magEpoch[0][startX:],
        magEpoch[1][startY:],
        magEpoch[2][startZ:]
    ]


    # --- --- --- --- --- --- --- --- ---
    # --- PROCESSING AND INTERPOLATE ---
    # --- --- --- --- --- --- --- --- ---













    # --- --- --- --- --- --- --- ---
    # --- prepare data for output ---
    # --- --- --- --- --- --- --- ---
    del data_dict
    data_dict = {}

    data_dict = {**data_dict, **{'MagX': [magDataT[0], deepcopy(data_dict_temp['data'][1])]}}
    data_dict = {**data_dict, **{'MagX_Epoch': [magEpoch[0], deepcopy(data_dict_temp['Epoch'][1])]}}
    data_dict['MagX'][1]['DEPEND_0'] = 'MagX_Epoch'
    data_dict['MagX'][1]['VAR_TYPE'] = 'data'
    data_dict['MagX'][1]['UNITS'] = 'nanotesla'
    data_dict['MagX'][1]['LABLAXIS'] = 'MagX'

    data_dict = {**data_dict, **{'MagY': [magDataT[1], deepcopy(data_dict_temp['data'][1])]}}
    data_dict = {**data_dict, **{'MagY_Epoch': [magEpoch[1], deepcopy(data_dict_temp['Epoch'][1])]}}
    data_dict['MagY'][1]['DEPEND_0'] = 'MagY_Epoch'
    data_dict['MagY'][1]['VAR_TYPE'] = 'data'
    data_dict['MagY'][1]['UNITS'] = 'nanotesla'
    data_dict['MagY'][1]['LABLAXIS'] = 'MagY'

    data_dict = {**data_dict, **{'MagZ': [magDataT[2], deepcopy(data_dict_temp['data'][1])]}}
    data_dict = {**data_dict, **{'MagZ_Epoch': [magEpoch[2], deepcopy(data_dict_temp['Epoch'][1])]}}
    data_dict['MagZ'][1]['DEPEND_0'] = 'MagZ_Epoch'
    data_dict['MagZ'][1]['VAR_TYPE'] = 'data'
    data_dict['MagZ'][1]['UNITS'] = 'nanotesla'
    data_dict['MagZ'][1]['LABLAXIS'] = 'MagZ'




    # --- --- --- --- ---
    # --- OUTPUT Data ---
    # --- --- --- --- ---
    prgMsg('Writing out Data')
    fileoutName = f'TRICE_{missionAttrs.rocketID[wRocket]}_deSpin'
    outputPath = f'{TRICE_data_folder}\\science\\{fileoutName}.cdf'

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
            deSpin(i,wComponent)
        else:
            print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if path.exists(TRICE_mag_files[wRocket[0]][0]):
        deSpin(wRocket[0],wComponent)
    else:
        print(color.RED + 'There are no .cdf files in the specified directory' + color.END)