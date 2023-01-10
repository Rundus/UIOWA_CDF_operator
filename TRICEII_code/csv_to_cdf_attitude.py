# --- csv_to_cdf_attitude.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Turn the .cdf files of the TRICE attitude data into cdf files


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
import cdflib
from os import remove,path
from csv import reader
from data_paths import TRICE_attitude_folder,fliers
from class_var_func import csv_attitude_TRICE,color,L0_TRICE
from missionAttributes import TRICE_mission_dicts
from copy import deepcopy
setupPYCDF()
from spacepy import pycdf

print(color.BOLD + color.CYAN + 'csv_to_cdf_attitude.py' + color.END + color.END)

def csv_to_cdf(wRocket):
    inputFile =csv_attitude_TRICE()
    with pycdf.CDF(inputFile(wRocket)[0]) as csvFile:
        csvAllData = [row for row in reader(csvFile)]

    varNames = csvAllData[4]
    csvData = np.array(csvAllData[5:],dtype='float32').transpose()

    # --- Create data Dict ---
    rocketAttrs,missionDict,data_dict_temp = TRICE_mission_dicts()
    data_dict = {varNames[i].replace(' ','') : [ csvData[i] ,  deepcopy(data_dict_temp['data'][1])  ] for i in range(len(varNames))}
    prgMsg(f'Converting csv Files for {rocketAttrs.rocketID[wRocket]}')


    # --- Modify Variable Attributes (the important ones) ---

    # --- time since launch ---
    data_dict['time'][1]['UNITS'] = 'seconds'

    # --- Epoch ---
    startval = cdflib.epochs.CDFepoch.compute_tt2000(rocketAttrs.epoch_starts[wRocket])
    vardata = np.array([(float(value) * (10 ** 9) + startval) for value in data_dict['time'][0]])
    data_dict = {**data_dict,**{  'Epoch' : [vardata,  data_dict_temp['Epoch'][1] ]}   }

    # --- Roll ---
    data_dict['roll'][1]['UNITS'] = 'Degrees'
    data_dict['roll'][0][0] = 0 #remove pesky "Nan" value from the start of the data.


    # --- --- --- --- ---
    # --- Output Data ---
    # --- --- --- --- ---
    fileoutName = f'TRICE_{rocketAttrs.rocketID[wRocket]}_attitude_data_p0.cdf'
    outputPath = f'{TRICE_attitude_folder}\\{fliers[wRocket]}\\{fileoutName}'

    # --- delete output file if it already exists ---
    if path.exists(outputPath):
        remove(outputPath)
    pycdf.lib.set_backward(False)


    # --- open the output file ---
    with pycdf.CDF(outputPath, '') as cdfFile:
        cdfFile.readonly(False)

        # --- WRITE OUT GLOBAL ATTRS ---
        globalAttrsMod = rocketAttrs.globalAttributes[wRocket]
        ModelData = L0_TRICE(wRocket)[0]
        inputGlobDic = ModelData.cdfFile.globalattsget()

        for key, val in inputGlobDic.items():
            if key in globalAttrsMod:
                cdfFile.attrs[key] = globalAttrsMod[key]
            else:
                cdfFile.attrs[key] = val

        # --- WRITE OUT DATA ---
        for varKey, varVal in data_dict.items():
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



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == []:
    for i in range(2):
        if csv_attitude_TRICE(i)[0] == []:
            print(color.RED + 'There are no .csv files in the specified directory' + color.END)
        else:
            csv_to_cdf(i)
else:
    if csv_attitude_TRICE(wRocket[0])[0] == []:
        print(color.RED + 'There are no .csv files in the specified directory' + color.END)
    else:
        csv_to_cdf(wRocket[0])