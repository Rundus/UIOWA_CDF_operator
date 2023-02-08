# --- L2_Langmuir_to_Temp_&_Density.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Takes in the L2 Langmuir Data and performs the Chi Square analysis to get the Temperature and
# Density from the characteristic curves



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

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'L2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\Langmuir' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = False


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os
import matplotlib.pyplot as plt
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg, L2_ACES_Quick
from glob import glob
from os.path import getsize

setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


def main(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'LangmuirData'
    outputModelData = L2_ACES_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*langmuir*')
    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*Temp&Density*')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')
    fileoutName = dataFile_name.replace('langmuir_','langmuir_Temp&Density_')

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made L2: {:3s} '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict = {}
        with pycdf.CDF(inputFiles[wFile]) as inputDataFile:
            for key, val in inputDataFile.items():
                data_dict = {**data_dict, **{key : [inputDataFile[key][...] , {key:val for key,val in inputDataFile[key].attrs.items()  }  ]  }  }

        data_dict['Epoch_swept_Current'][0] =  np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_swept_Current'][0][i]) for i in (range(len(data_dict['Epoch_swept_Current'][0])  )   )])

        Done(start_time)

        ######################################
        # --- FIND THE SWEPT CURVE INDICES ---
        ######################################

        prgMsg('Finding the swept curve indicies')

        # Find the start and end point for each sweep
        indvEpochThresh = 250000000  # Value, in tt2000, that determines the time diff needed between epoch points to identify particular sweeps
        counter = 0
        sweepIndices = []

        for i in (range(len(data_dict['Epoch_swept_Current'][0]) - 1)):

            if np.abs((data_dict['Epoch_swept_Current'][0][i + 1]/100000 - data_dict['Epoch_swept_Current'][0][i]/100000)) >= indvEpochThresh/100000:

                if counter == 0:
                    sweepIndices.append([0,i])
                    counter += 1
                else:
                    sweepIndices.append([sweepIndices[-1][1]+1, i])

        # include the last sweep
        sweepIndices.append([sweepIndices[-1][1] + 1,len(data_dict['Epoch_swept_Current'][0])])



        # Find the index corresponding to the break in upleg/downleg voltage for the sweeps
        breakIndices = []
        for i in range(len(sweepIndices)):
            start = sweepIndices[i][0]
            end = sweepIndices[i][1]
            sweeps = np.array(data_dict['swept_Current'][0][start:end])
            sweepMAX, sweepMIN = sweeps.max(), sweeps.min()
            threshold = (sweepMAX - sweepMIN)*0.50

            for j in range(len(sweeps)-1):
                if np.abs(sweeps[j+1] - sweeps[j]) >= threshold:
                    breakIndices.append(j+sweepIndices[i][0])

        Done(start_time)

        sweepIndices = np.array(sweepIndices)
        breakIndices = np.array(breakIndices)


        ###############################
        # --- Perform ChiSquare Fit ---
        ###############################



        if outputData:
            # --- --- --- --- --- --- ---
            # --- WRITE OUT THE DATA ---
            # --- --- --- --- --- --- ---
            prgMsg('Creating output file')

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            # --- delete output file if it already exists ---
            if os.path.exists(outputPath):
                os.remove(outputPath)

            # --- open the output file ---
            with pycdf.CDF(outputPath, '') as L2File:
                L2File.readonly(False)

                # --- write out global attributes ---
                inputGlobDic = outputModelData.cdfFile.globalattsget()
                for key, val in inputGlobDic.items():
                    if key in globalAttrsMod:
                        L2File.attrs[key] = globalAttrsMod[key]
                    else:
                        L2File.attrs[key] = val

                # --- WRITE OUT DATA ---
                for varKey, varVal in data_dict.items():
                    if 'Epoch' in varKey: # epoch data
                        L2File.new(varKey, data=varVal[0], type=33)
                    else: # other data
                        L2File.new(varKey, data=varVal[0],type=pycdf.const.CDF_REAL8)

                    # --- Write out the attributes and variable info ---
                    for attrKey, attrVal in data_dict[varKey][1].items():
                        if attrKey == 'VALIDMIN':
                            L2File[varKey].attrs[attrKey] = varVal[0].min()
                        elif attrKey == 'VALIDMAX':
                            L2File[varKey].attrs[attrKey] = varVal[0].max()
                        elif attrVal != None:
                            L2File[varKey].attrs[attrKey] = attrVal

            Done(start_time)







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

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        main(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        main(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)