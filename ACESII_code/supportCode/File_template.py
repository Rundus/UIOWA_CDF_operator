# --- template.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:



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
wFiles = []

modifier = ''
inputPath_modifier = '' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = '' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg
from glob import glob
from os.path import getsize

setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


def main(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    if wRocket in [0, 1, 4, 5]:
        # --- ACES II Flight/Integration Data ---
        rocketAttrs, b, c = ACES_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
        globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
        outputModelData = outputModelData(wflyer)

    # --- TRICE II ---
    elif wRocket in [2, 3]:
        globalAttrsMod = {}
        rocketAttrs, b, c = TRICE_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        outputModelData = outputModelData(wflyer)




    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')
    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')
    fileoutName = dataFile_name.replace(inputPath_modifier.lower(), outputPath_modifier.lower())



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

        data_dict['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][i]) for i in (range(data_dict['Epoch'][0]))])

        Done(start_time)

        #######################
        # --- DO STUFF HERE ---
        #######################








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
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            main(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            main(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)