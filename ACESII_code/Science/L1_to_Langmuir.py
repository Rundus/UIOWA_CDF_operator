# --- L1_to_Langmuir.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert the engineering Langmuir data to scientifically useful units


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools
# --- --- --- --- ---

import time
from ACESII_code.class_var_func import Done, setupPYCDF

start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, L2_ACES_Quick,L2_TRICE_Quick, prgMsg
from glob import glob
from os.path import getsize
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)

from collections import Counter
import matplotlib.pyplot as plt


def L1_to_Langmuir(wRocket, rocketFolderPath, justPrintFileNames, wflyer):

    scienceFolderPath = rocketFolderPath + r'science\\'

    if wRocket in [0,1,4,5]:

        # --- ACES II Flight/Integration Data ---
        rocketAttrs,b,c = ACES_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
        globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'Langmuir'
        L2ModelData = L2_ACES_Quick(wflyer)

    # --- TRICE II ---
    elif wRocket in [2, 3]:
        globalAttrsMod = {}
        rocketAttrs,b,c = TRICE_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        L0ModelData = L2_TRICE_Quick(wflyer)


    # Set the paths for the file names
    L1Files = glob(f'{rocketFolderPath}L1\{fliers[wflyer]}\*_lp_*')
    LangmuirFiles = glob(f'{scienceFolderPath}\{fliers[wflyer]}\*_langmuir_*')

    L1_names = [ifile.replace(f'{rocketFolderPath}L1\{fliers[wflyer]}\\', '') for ifile in L1Files]
    Langmuir_names = [ofile.replace(f'{scienceFolderPath}\{fliers[wflyer]}\\', '') for ofile in LangmuirFiles]

    L1_names_searchable = [ifile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace('l1_', '').replace('_v00', '') for ifile in L1_names]
    Langmuir_names_searchable = [ofile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace('_langmuir_', '').replace('_v00', '').replace('__', '_') for ofile in Langmuir_names]

    dataFile_name = L1_names_searchable[0].replace(f'{rocketFolderPath}L1\{fliers[wflyer]}\\', '').replace('lp_','').replace('ni_','').replace('ne_swept_','').replace('step_','').replace('ni_swept_','').replace('deltaNdivN_','')
    fileoutName = rf'ACESII_{rocketAttrs.rocketID[wRocket-4]}_langmuir_{dataFile_name}'






    if justPrintFileNames:
            for i, file in enumerate(L1Files):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, L1_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Creating langmuir science data for {fliers[wRocket-4]} flyer' + color.END)

        # --- get the data from the tmCDF file ---
        prgMsg('Loading data from L1Files')
        data_dict = {}

        # Collect the LP data except deltaNdivN into a data dict
        for file in L1Files:
            if 'deltaNdivN' not in file:
                with pycdf.CDF(file) as L1DataFile:
                    for key, val in L1DataFile.items():
                        if key not in data_dict:
                            data_dict = {**data_dict, **{key : [L1DataFile[key][...] , {key:val for key,val in L1DataFile[key].attrs.items()  }  ]  }  }


        Done(start_time)



        # --- --- --- --- --- --- --- --- --- --- --- --- ---
        # --- Determine the DAC to Voltage Conversions ---
        # --- --- --- --- --- --- --- --- --- --- --- --- ---

        ### Fixed Probe ###
        fixedCalResistances = rocketAttrs.LPFixed_calResistances[wRocket - 4]
        probeBias = rocketAttrs.LPFixedProbeBias[wRocket - 4]

        ### Swept Probe ###
        voltageRange = rocketAttrs.LPswept_voltage_range[wRocket-4]

        # calculate the epoch index of beginning/end of sample range
        Epoch_step = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_step'][0][i]) for i in range(len(data_dict['Epoch_step'][0]))])
        sampleStart = np.abs(np.array(Epoch_step - pycdf.lib.datetime_to_tt2000(rocketAttrs.Epoch_range_to_determine_stepDAC[wRocket-4][0]))).argmin()
        sampleEnd = np.abs(np.array(Epoch_step - pycdf.lib.datetime_to_tt2000(rocketAttrs.Epoch_range_to_determine_stepDAC[wRocket-4][1]))).argmin()
        adjustments = [[2,2],[-1,-1]]
        stepSample = Counter(data_dict['step'][0][sampleStart:sampleEnd])

        # determine the values of step for each step
        sampleData = data_dict['step'][0][sampleStart+adjustments[wRocket-4][0]:sampleEnd+adjustments[wRocket-4][1]]

        stepsDigitalVals = [round(sum(sampleData[i*10:(i+1)*10])/10) for i in range(round(len(sampleData)/10))]

        counted = Counter(stepsDigitalVals)
        counted_sort = dict(sorted(counted.items()))
        counted_sort_reduced = dict(sorted(counted.items()))
        keys = []
        aside = []

        for key,val in counted_sort_reduced.items():
            keys.append(key)

        for i in range(len(keys)-1):
            if np.abs(keys[i+1] - keys[i]) < 10:
                if counted_sort[keys[i+1]] > counted_sort[keys[i]]:
                    counted_sort_reduced.pop(keys[i], None)
                elif counted_sort[keys[i]] < counted_sort[keys[i]]:
                    counted_sort_reduced.pop(keys[i], None)
                elif counted_sort[keys[i]] == counted_sort[keys[i]]:
                    counted_sort_reduced.pop(keys[i+1], None)
            if counted_sort[keys[i]] < 5:
                counted_sort_reduced.pop(keys[i], None)

        # 4
        # 15
        # 31
        # 46
        # 62
        # 78
        # 94
        # 110
        # 126

        # for key,val in counted_sort.items():
        #     if val <=5:
        #         counted_sort_reduced.pop(key, None)
        #     if val == 10:
        #         aside.append(key)

        for key, val in counted_sort_reduced.items():
            print(key,val)

        print(len(counted_sort_reduced))
        # # --- --- --- --- --- --- --- --- ---
        # # --- Calculate Instrument Data ---
        # # --- --- --- --- --- --- --- --- ---
        #
        # prgMsg('Creating Langmuir science data')
        #





        #
        # # --- --- --- --- --- --- ---
        # # --- WRITE OUT THE DATA ---
        # # --- --- --- --- --- --- ---
        # prgMsg('Creating output file')
        #
        # outputPath = f'{rocketFolderPath}L2\{fliers[wflyer]}\\{fileoutName}'
        #
        #
        # # --- delete output file if it already exists ---
        # if os.path.exists(outputPath):
        #     os.remove(outputPath)
        #
        # # --- open the output file ---
        # with pycdf.CDF(outputPath, '') as LangmuirFile:
        #     LangmuirFile.readonly(False)
        #
        #     # --- write out global attributes ---
        #     inputGlobDic = L2ModelData.cdfFile.globalattsget()
        #     for key, val in inputGlobDic.items():
        #         if key == 'Descriptor':
        #             globalAttrsMod[key] = rocketAttrs.InstrNames_Full[wInstr[0]]
        #         if key in globalAttrsMod:
        #             LangmuirFile.attrs[key] = globalAttrsMod[key]
        #         else:
        #             LangmuirFile.attrs[key] = val
        #
        #     # --- WRITE OUT DATA ---
        #     for varKey, varVal in data_dict.items():
        #         if 'Epoch' in varKey: # epoch data
        #             LangmuirFile.new(varKey, data=varVal[0], type=33)
        #         else: # other data
        #             LangmuirFile.new(varKey, data=varVal[0],type=pycdf.const.CDF_REAL8)
        #
        #         # --- Write out the attributes and variable info ---
        #         for attrKey, attrVal in data_dict[varKey][1].items():
        #             if attrKey == 'VALIDMIN':
        #                 LangmuirFile[varKey].attrs[attrKey] = varVal[0].min()
        #             elif attrKey == 'VALIDMAX':
        #                 LangmuirFile[varKey].attrs[attrKey] = varVal[0].max()
        #             elif attrVal != None:
        #                 LangmuirFile[varKey].attrs[attrKey] = attrVal
        #
        #     Done(start_time)










# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 0:  # ACES II Integration High
    rocketFolderPath = Integration_data_folder
    wflyer = 0
elif wRocket == 1: # ACES II Integration Low
    rocketFolderPath = Integration_data_folder
    wflyer = 1
elif wRocket == 2:  # TRICE II High
    rocketFolderPath = TRICE_data_folder
    wflyer = 0
elif wRocket == 3: # TRICE II Low
    rocketFolderPath = TRICE_data_folder
    wflyer = 1
elif wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}L1\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        L1_to_Langmuir(wRocket, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        L1_to_Langmuir(wRocket, rocketFolderPath, justPrintFileNames,wflyer)