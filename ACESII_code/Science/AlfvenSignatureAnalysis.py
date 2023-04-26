# --- AlfvenSignatureAnalysis.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Pull the L1 eepaa data and look at the specific locations of 17:24:55 to 17:25:10
# and perform an analysis on this dataset to study them further



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import copy
import itertools
# --- --- --- --- ---

import time

import matplotlib.pyplot as plt

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
wFiles = [7]

modifier = ''
inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\AlfvenSignatureAnalysis' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# Subtracting Inverted V
removeInvertedV = True # Tries to remove InvertedV from the data
maskVal = 10
strengthModifier = 1



# PLOTTING THE DATA
plotAlfvenRegion = True
wPitch = 1
epochRange = [[17,24,55],[17,25,12]] # UTC Ranges for plot Epoch. Format: [hour, minute,second]
vExtremes = [2, 20] #colorbar limits
yLimits = [8,1050] # limits of y-axis in eV

# NEXT THING
outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import datetime as dt
import os
import matplotlib.ticker as ticker
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg,outputCDFdata,L1_TRICE_Quick,L2_TRICE_Quick
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
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)


    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')
    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')
    fileoutName = dataFile_name.replace(inputPath_modifier.lower(), outputPath_modifier.lower())

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            descriptiorNam = ['EEPAA', 'LEESA', 'IEPAA', 'Langmuir_Probe']
            wInstr = [index, instr, descriptiorNam[index]]


    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
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

        data_dict['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_esa'][0][i]) for i in (range(len(data_dict['Epoch_esa'][0])))])

        Done(start_time)

        ############################
        # --- REDUCE THE DATASET ---
        ############################
        # Reduce Epoch to only the interested alfven signature range
        targetTimes = [pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, epochRange[0][0], epochRange[0][1], epochRange[0][2], 000000)),
                       pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, epochRange[1][0], epochRange[1][1], epochRange[1][2], 000000))]

        Epoch = data_dict['Epoch_esa'][0]
        limitIndexes = [np.abs(Epoch - targetTimes[0]).argmin(), np.abs(Epoch - targetTimes[1]).argmin()]
        Epoch_reduced = np.array([pycdf.lib.tt2000_to_datetime(Epoch[i]) for i in range(limitIndexes[0], limitIndexes[1])])  # give reduced epoch in datetimes

        # Reduce esaData to only one pitch angle
        esaData = copy.deepcopy(data_dict[wInstr[1]][0][limitIndexes[0]:limitIndexes[1]])

        # pull out other variables
        Energy = data_dict['Energy'][0]
        Pitch = data_dict['Pitch_Angle'][0]

        ############################
        # --- REMOVED INVERTED V ---
        ############################
        if removeInvertedV:

            # get an average value for the ESA
            prgMsg('Removing Inverted V')

            # get the initial esaData we're interested in
            esaData_raw = np.array([esaData[i][wPitch][:] for i in range(len(esaData))])
            esaData = np.array([esaData[i][wPitch][:] for i in range(len(esaData))])


            # Apply a counts mask
            for tme in range(len(esaData)):
                for engy in range(len(esaData[0])):
                    val = esaData[tme][engy] - maskVal

                    if val <= 0:
                        esaData[tme][engy] = 0
                    else:
                        esaData[tme][engy] = val


            # look for values that have no counts +/- 1 tme step away from them. Remove these points
            esaData_temp = esaData

            for tme in range(2,len(esaData) -2):
                for engy in range(len(esaData[0])):
                    if esaData_temp[tme - 1][engy] == 0 and esaData_temp[tme+1][engy] == 0:
                        esaData[tme][engy] = 0

                    if (esaData_temp[tme - 1][engy] == 0 and esaData_temp[tme - 2][engy] == 0) or (esaData_temp[tme + 1][engy] == 0 and esaData_temp[tme + 2][engy] == 0):
                        esaData[tme][engy] = 0


            esaData_masked = esaData_raw

            # remove the mask we created
            for tme in range(len(esaData)):
                for engy in range(len(esaData[0])):
                    val = esaData_raw[tme][engy] - strengthModifier*esaData[tme][engy]
                    if val <= 0:
                        esaData_masked[tme][engy] = 0
                    else:
                        esaData_masked[tme][engy] = val

            Done(start_time)

        #######################
        # --- PLOT THE DATA ---
        #######################
        if plotAlfvenRegion:

            prgMsg('Plotting AlfvenSignature Region')
            if removeInvertedV:
                fig, ax = plt.subplots(2, sharex=True)

                cmap = ax[0].pcolormesh(Epoch_reduced, Energy, esaData.transpose(), vmin=vExtremes[0], vmax=vExtremes[1], cmap="turbo")
                cmap = ax[1].pcolormesh(Epoch_reduced, Energy, esaData_masked.transpose(), vmin=vExtremes[0], vmax=vExtremes[1], cmap="turbo")

            else:
                fig, ax = plt.subplots()
                # colormap
                cmap = ax.pcolormesh(Epoch_reduced, Energy, esaData.transpose(), vmin=vExtremes[0], vmax=vExtremes[1], cmap="turbo")

            # colorbar

            # limits
            if removeInvertedV:
                for i in range(2):
                    ax[i].set_ylim(yLimits[0], yLimits[1])

                    # labels
                    ax[i].set_ylabel('Energy [eV]')

                    cbar = plt.colorbar(cmap, ax = ax[i])
                    cbar.minorticks_on()
                    cbar.set_label('Counts')

                    if i == 0:
                        ax[i].set_title(f'ACESII {rocketAttrs.rocketID[wRocket - 4]} \n {wInstr[1].upper()} \n Pitch {Pitch[wPitch]}$^\circ$')
                    elif i == 1:
                        ax[i].set_xlabel('Epoch')
            else:

                ax.set_ylim(yLimits[0], yLimits[1])

                cbar = plt.colorbar(cmap)
                cbar.minorticks_on()
                cbar.set_label('Counts')

                # labels
                ax.set_ylabel('Energy [eV]')
                ax.set_xlabel('Epoch')
                ax.set_title(f'ACESII {rocketAttrs.rocketID[wRocket - 4]} \n {wInstr[1].upper()} \n Pitch {Pitch[wPitch]}$^\circ$')



            plt.show()
            Done(start_time)










        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        if outputData:
            prgMsg('Creating output file')

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, wInstr[2])

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