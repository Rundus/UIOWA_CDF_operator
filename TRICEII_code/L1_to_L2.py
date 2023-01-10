# --- L1_to_L2.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert TRICE electrostatic analyzer counts data from counts to differential energy flux

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools
# --- --- --- --- ---

import time
from class_var_func import Done, setupPYCDF

start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 0 -> TRICE II High Flier
# 1 -> TRICE II Low Flier
wRocket = 1

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [1]

# Special case to use the p2 data instead of fully processed data
usep2 = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os

from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from tqdm import tqdm
from missionAttributes import TRICE_mission_dicts
from data_paths import TRICE_data_folder, TRICE_data_folder, fliers
from class_var_func import color, prgMsg, L1_TRICE,L2_TRICE
from glob import glob
from os.path import getsize
setupPYCDF()
from spacepy import pycdf


def L1_to_L2(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    if wRocket in [0,1]:

        # --- TRICE II Flight/Integration Data ---
        rocketAttrs,missionDict,data_dict_template = TRICE_mission_dicts()
        globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
        globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
        L2ModelData = L2_TRICE(wflyer)[0]


    # Set the paths for the file names
    L1Files = glob(f'{rocketFolderPath}L1\{fliers[wflyer]}\*.cdf')
    L2Files = glob(f'{rocketFolderPath}L2\{fliers[wflyer]}\*.cdf')

    L1_names = [ifile.replace(f'{rocketFolderPath}L1\{fliers[wflyer]}\\', '') for ifile in L1Files]
    L2_names = [ofile.replace(f'{rocketFolderPath}L2\{fliers[wflyer]}\\', '') for ofile in L2Files]

    L1_names_searchable = [ifile.replace('TRICE_', '').replace('36359_', '').replace('36364_', '').replace('l1_', '').replace('_v00', '') for ifile in L1_names]
    L2_names_searchable = [ofile.replace('TRICE_', '').replace('36359_', '').replace('36364_', '').replace('l2_', '').replace('_v00', '').replace('__', '_') for ofile in L2_names]

    dataFile_name = L1Files[wFile].replace(f'{rocketFolderPath}L1\{fliers[wflyer]}\\', '')
    fileoutName = dataFile_name.replace('l1', 'l2')


    if justPrintFileNames:
            for i, file in enumerate(L1Files):
                anws = ["yes" if L1_names_searchable[i].replace('.cdf', "") in L2_names_searchable else "no"]
                print('[{:.0f}] {:55s}{:5.1f} MB   Made L2: {:3s} '.format(i, L1_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
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

        if usep2:
            sizes = [len(data_dict[rocketAttrs.InstrNames_LC][0]), len(data_dict['Pitch_Angle'][0]),
                     len(data_dict['Energy'][0])]
        else:
            sizes = [len(data_dict[rocketAttrs.InstrNames_LC][0]), len(data_dict['Pitch_Angle'][0]), len(data_dict['Energy'][0])]

        ranges = [range(sizes[0]), range(sizes[1]), range(sizes[2])]
        diffNFlux = np.zeros(shape=(sizes[0], sizes[1], sizes[2]))
        diffEFlux = np.zeros(shape=(sizes[0], sizes[1], sizes[2]))

        Energies = data_dict['Energy'][0]

        if usep2:
            counts = data_dict['counts_p2'][0]
        else:
            counts = data_dict[rocketAttrs.InstrNames_LC][0]


        count_interval = data_dict['Count_Interval'][0]
        mean_count_interval = sum(count_interval)/len(count_interval)
        for i in range(len(count_interval)):
            if count_interval[i] <= 0:
                count_interval[i] = mean_count_interval


        geo_factor = rocketAttrs.geometric_factor[wRocket]

        # --- PROCESS ESA DATA ---
        for tme, ptch, engy in tqdm(itertools.product(*ranges)):
            try:
                diffNFlux[tme][ptch][engy] = int((counts[tme][ptch][engy])/(Energies[engy] * geo_factor[ptch] * ( (count_interval[tme]*10**(-6)) - (counts[tme][ptch][engy] * rocketAttrs.deadtime[wRocket]))))
                diffEFlux[tme][ptch][engy] = int((Energies[engy] * counts[tme][ptch][engy]) / (Energies[engy] * geo_factor[ptch] * ( (count_interval[tme]*10**(-6)) - (counts[tme][ptch][engy] * rocketAttrs.deadtime[wRocket]))))
            except:
                print(counts[tme][ptch][engy])
                print(Energies[engy])
                print(geo_factor[ptch])
                print(count_interval[tme])
                print(rocketAttrs.deadtime[wRocket])


        del data_dict[rocketAttrs.InstrNames_LC]

        if usep2:
            del data_dict['counts_p2']
        else:
            del data_dict['minor_frame_counter']
            del data_dict['major_frame_counter']


        outputData = True


        Done(start_time)

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        prgMsg('Creating output file')

        if outputData:

            outputPath = f'{rocketFolderPath}L2\{fliers[wflyer]}\\{fileoutName}'

            data_dict = {**data_dict, **{'Differential_Number_Flux':
                                             [diffNFlux, {'LABLAXIS': 'Differential_Number_Flux',
                                                       'DEPEND_0': 'Epoch', 'DEPEND_1': 'Pitch_Angle',
                                                       'DEPEND_2': 'Energy',
                                                       'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                       'UNITS': 'cm!U-2!N str!U-1!N s!U-1!N eV!U-1!N',
                                                       'VALIDMIN': diffNFlux.min(), 'VALIDMAX': diffNFlux.max(),
                                                       'VAR_TYPE': 'data', 'SCALETYP': 'log'}]}}

            data_dict = {**data_dict, **{'Differential_Energy_Flux':
                                             [diffNFlux, {'LABLAXIS': 'Differential_Energy_Flux',
                                                          'DEPEND_0': 'Epoch', 'DEPEND_1': 'Pitch_Angle',
                                                          'DEPEND_2': 'Energy',
                                                          'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                          'UNITS': 'cm!U-2!N str!U-1!N s!U-1!N eV/eV',
                                                          'VALIDMIN': diffEFlux.min(), 'VALIDMAX': diffEFlux.max(),
                                                          'VAR_TYPE': 'data', 'SCALETYP': 'log'}]}}

            # --- delete output file if it already exists ---
            if os.path.exists(outputPath):
                os.remove(outputPath)

            pycdf.lib.set_backward(False)

            # --- open the output file ---
            with pycdf.CDF(outputPath, '') as L2File:
                L2File.readonly(False)

                # --- write out global attributes ---
                inputGlobDic = L2ModelData.cdfFile.globalattsget()
                for key, val in inputGlobDic.items():
                    if key == 'Descriptor':
                        globalAttrsMod[key] = rocketAttrs.InstrNames_Full
                    if key in globalAttrsMod:
                        L2File.attrs[key] = globalAttrsMod[key]
                    else:
                        L2File.attrs[key] = val

                # --- WRITE OUT DATA ---
                for varKey, varVal in data_dict.items():
                    if varKey in ['Epoch', 'Epoch_monitors', 'Epoch_esa']: # epoch data
                        L2File.new(varKey, data=varVal[0], type=33)
                    # elif varKey in ['Differential_Number_Flux', 'Differential_Energy_Flux']: #instrument data
                    #     L2File.new(varKey, data=varVal[0], type=pycdf.const.CDF_DOUBLE)
                    else: # support data
                        L2File.new(varKey, data=varVal[0])

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
if wRocket == 0:  # TRICE II High
    rocketFolderPath = TRICE_data_folder
    wflyer = 0
elif wRocket == 1: # TRICE II Low
    rocketFolderPath = TRICE_data_folder
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