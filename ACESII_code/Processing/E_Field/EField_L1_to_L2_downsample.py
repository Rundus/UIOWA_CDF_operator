# --- EField_L1_to_L2_downsample.py.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: The E_Field is too highly sampled. Downsample it to match its respective MAG counterpart timebase



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import time

import numpy as np

from ACESII_code.myImports import *

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
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [2]


modifier = ''
inputPath_modifier_E_ENU = r'\l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_B_ENU =r'l2'
outputPath_modifier = r'\l2' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import math
def EField_L1_to_L2_downsample(wRocket, wFile, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    globalAttrsMod = rocketAttrs.globalAttributes[wRocket-4]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wRocket-4)

    inputFiles_EENU = glob(f'{rocketFolderPath}{inputPath_modifier_E_ENU}\{fliers[wRocket-4]}{modifier}\*.cdf')
    inputFiles_BENU = glob(f'{rocketFolderPath}{inputPath_modifier_B_ENU}\{fliers[wRocket-4]}{modifier}\*RingCore_ENU.cdf*')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier_E_ENU}\{fliers[wRocket-4]}{modifier}\\', '') for ifile in inputFiles_EENU]
    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier_E_ENU.lower() +'_', '').replace('_v00', '') for ifile in input_names]

    fileoutName = input_names[wFile].replace('l1_','l2_').replace('_flight','_downsampled')


    if justPrintFileNames:
        for i, file in enumerate(inputFiles_EENU):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    # --- get the data from the attitude file ---
    prgMsg(f'Loading data from {inputPath_modifier_B_ENU} Files')
    data_dict_mag = loadDictFromFile(inputFiles_BENU[0])
    Done(start_time)


    # --- get the data from the E-Field or B-Field file ---
    prgMsg(f'Loading data from {inputPath_modifier_E_ENU} Files')
    data_dict_EField,globalAttrs = loadDictFromFile(inputFiles_EENU[wFile],getGlobalAttrs=True, targetVar=[[data_dict_mag['Epoch'][0][0],data_dict_mag['Epoch'][0][-1]],'Epoch'])
    Done(start_time)



    #################################
    # --- Downsample E-Field Data ---
    #################################
    prgMsg('DownSampling E-Field Data')

    E_Epoch =dateTimetoTT2000(InputEpoch= data_dict_EField['Epoch'][0],inverse=False)
    B_Epoch = dateTimetoTT2000(InputEpoch=data_dict_mag['Epoch'][0],inverse=False)

    def find_nearest(array, value):
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
            return idx - 1
        else:
            return idx

    indices = np.zeros(shape=(len(B_Epoch)),dtype='int64')
    for i in tqdm(range(len(B_Epoch))):
        indices[i] = int(find_nearest(E_Epoch,B_Epoch[i]))

    for key in data_dict_EField.keys():
        data_dict_EField[key][0] = deepcopy(data_dict_EField[key][0][indices])
    Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        prgMsg('Writing Out Data')

        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wRocket-4]}\\{fileoutName}'

        outputCDFdata(outputPath, data_dict_EField, globalAttrsMod=globalAttrs,instrNam='E-Field')

        Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder

if len(glob(f'{rocketFolderPath}{inputPath_modifier_E_ENU}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        EField_L1_to_L2_downsample(wRocket, 0, rocketFolderPath, justPrintFileNames)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier_E_ENU}\{fliers[wRocket-4]}\*.cdf')))):
            EField_L1_to_L2_downsample(wRocket, fileNo, rocketFolderPath, justPrintFileNames)
    else:
        for filesNo in wFiles:
            EField_L1_to_L2_downsample(wRocket, filesNo, rocketFolderPath, justPrintFileNames)