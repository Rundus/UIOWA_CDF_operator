# --- L0_to_L1.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert electrostatic analyzer data from counts to differential energy flux

# TODO: Need to determine: (1) geometric factor (2) Count Interval (3) Deadtime




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
wFiles = [0,1,7]

useMagCalData = True
if useMagCalData:
    inputPath_modifier = 'L1_mag_cal'  # e.g. 'L1' or 'L1'. It's the name of the broader input folder
    outputPath_modifier = 'L2_mag_cal'  # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
else:
    inputPath_modifier = 'L1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
    outputPath_modifier = 'L2' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

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
from ACESII_code.class_var_func import color, L2_TRICE_Quick, prgMsg
from glob import glob
from os.path import getsize
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


def L1_to_L2(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    if wRocket in [0,1,4,5]:
        # --- ACES II Flight/Integration Data ---
        rocketAttrs,b,c = ACES_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
        globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
        L2ModelData = L2_TRICE_Quick(wflyer)

    # --- TRICE II ---
    elif wRocket in [2, 3]:
        globalAttrsMod = {}
        rocketAttrs,b,c = TRICE_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        L0ModelData = L2_TRICE_Quick(wflyer)

    L1Files = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')
    L2Files = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    L1_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\\', '') for ifile in L1Files]
    L2_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in L2Files]

    L1_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace('l1_', '').replace('_v00', '') for ifile in L1_names]
    L2_names_searchable = [ofile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace('l2_', '').replace('_v00', '').replace('__', '_') for ofile in L2_names]

    dataFile_name = L1Files[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\\', '')
    fileoutName = dataFile_name.replace('l1', 'l2')

    # determine which instrument the file corresponds to:
    for index,instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            wInstr = [index,instr]


    if justPrintFileNames:
            for i, file in enumerate(L1Files):
                anws = ["yes" if L1_names_searchable[i].replace('.cdf', "") in L2_names_searchable else "no"]
                print('[{:.0f}] {:80s}{:5.1f} MB   Made L2: {:3s} '.format(i, L1_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to L2 data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(L1Files[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg('Loading data from L1Files')
        data_dict = {}
        with pycdf.CDF(L1Files[wFile]) as L0DataFile:
            for key, val in L0DataFile.items():
                data_dict = {**data_dict, **{key : [L0DataFile[key][...] , {key:val for key,val in L0DataFile[key].attrs.items()  }  ]  }  }


        Done(start_time)

        # --- --- --- --- --- --- --- --- ---
        # --- Calculate Instrument Data ---
        # --- --- --- --- --- --- --- --- ---

        prgMsg('Creating L2 instrument data')

        # --- PROCESS EEPAA/LEESA DATA ---
        if wInstr[0] in [0, 1, 2]:

            sizes = [len(data_dict[rocketAttrs.InstrNames_LC[wInstr[0]]][0]),len(data_dict['Pitch_Angle'][0]), len(data_dict['Energy'][0])]
            ranges = [range(sizes[0]), range(sizes[1]), range(sizes[2])]
            diffNFlux = np.zeros(shape=(sizes[0], sizes[1], sizes[2]))
            diffEFlux = np.zeros(shape=(sizes[0], sizes[1], sizes[2]))

            Energies = data_dict['Energy'][0]
            counts = data_dict[rocketAttrs.InstrNames_LC[wInstr[0]]][0]
            count_interval = data_dict['Count_Interval'][0]
            count_interval = [917 for i in range(len(data_dict[rocketAttrs.InstrNames_LC[wInstr[0]]][0]))]
            geo_factor = rocketAttrs.geometric_factor[wInstr[0]]


            # --- PROCESS ESA DATA ---
            for tme, ptch, engy in tqdm(itertools.product(*ranges)):
                deltaT = (count_interval[tme] * 10 ** (-6)) - (counts[tme][ptch][engy] * rocketAttrs.deadtime[wflyer])
                diffNFlux[tme][ptch][engy] = int((counts[tme][ptch][engy]) / (Energies[engy] * geo_factor[ptch] * deltaT))
                diffEFlux[tme][ptch][engy] = int((counts[tme][ptch][engy]) / (geo_factor[ptch] * deltaT))

            del data_dict[rocketAttrs.InstrNames_LC[wInstr[0]]]
            outputData = True

        # --- PROCESS LP DATA ---
        elif wInstr[0] in [3]:
            print('potato')
            outputData = False

        Done(start_time)

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        prgMsg('Creating output file')

        if outputData:

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            data_dict = {**data_dict, **{'Differential_Number_Flux':
                                             [diffNFlux, {'LABLAXIS': 'Differential_Number_Flux',
                                                       'DEPEND_0': 'Epoch_esa', 'DEPEND_1': 'Pitch_Angle',
                                                       'DEPEND_2': 'Energy',
                                                       'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                       'UNITS': 'cm!U-2!N str!U-1!N s!U-1!N eV!U-1!N',
                                                       'VALIDMIN': diffNFlux.min(), 'VALIDMAX': diffNFlux.max(),
                                                       'VAR_TYPE': 'data', 'SCALETYP': 'log'}]}}

            data_dict = {**data_dict, **{'Differential_Energy_Flux':
                                             [diffEFlux, {'LABLAXIS': 'Differential_Energy_Flux',
                                                          'DEPEND_0': 'Epoch_esa', 'DEPEND_1': 'Pitch_Angle',
                                                          'DEPEND_2': 'Energy',
                                                          'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                          'UNITS': 'cm!U-2!N str!U-1!N s!U-1!N eV/eV',
                                                          'VALIDMIN': diffEFlux.min(), 'VALIDMAX': diffEFlux.max(),
                                                          'VAR_TYPE': 'data', 'SCALETYP': 'log'}]}}

            # --- delete output file if it already exists ---
            if os.path.exists(outputPath):
                os.remove(outputPath)

            # --- open the output file ---
            with pycdf.CDF(outputPath, '') as L2File:
                L2File.readonly(False)

                # --- write out global attributes ---
                inputGlobDic = L2ModelData.cdfFile.globalattsget()
                for key, val in inputGlobDic.items():
                    if key == 'Descriptor':
                        globalAttrsMod[key] = rocketAttrs.InstrNames_Full[wInstr[0]]
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

if len(glob(f'{rocketFolderPath}L1\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        L1_to_L2(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}L1\{fliers[wflyer]}\*.cdf')))):
            L1_to_L2(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            L1_to_L2(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)