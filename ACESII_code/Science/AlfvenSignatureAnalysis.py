# --- AlfvenSignatureAnalysis.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Pull the L1 eepaa data and look at the specific locations of the alfven signature
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
wFiles = [0]

modifier = ''
inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_trajectory = 'trajectories'
outputPath_modifier = 'science\AlfvenSignatureAnalysis' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# Subtracting Inverted V
removeInvertedV = True # Tries to remove InvertedV from the data



# PLOTTING THE DATA
plotAlfvenRegion = True # plots the specific wPitch esa data
wPitch = 2
epochRange = [[17, 24, 56], [17, 25, 9]] # UTC Ranges for plot Epoch. Format: [hour, minute,second]
vExtremes = [1, 50] #colorbar limits
yLimits = [20, 1050] # limits of y-axis in eV


# NEXT THING
outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import datetime as dt
import os
import matplotlib.ticker as ticker
from matplotlib import figure
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg,outputCDFdata,L1_TRICE_Quick,L2_TRICE_Quick,loadDictFromFile
from glob import glob
from os.path import getsize


setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


def AlfvenSignatureAnalysis(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')
    inputFiles_trajectory = glob(f'{rocketFolderPath}{inputPath_modifier_trajectory}\{fliers[wflyer]}{modifier}\*_ILat_ILong*')
    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    input_names_trajectory = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles_trajectory]
    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')
    dataFile_name_trajectory = inputFiles_trajectory[0].replace(f'{rocketFolderPath}{inputPath_modifier_trajectory}\{fliers[wflyer]}{modifier}\\', '')
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

        # --- get the data from the ESA file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict = loadDictFromFile(inputFiles[wFile],{})
        data_dict['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_esa'][0][i]) for i in (range(len(data_dict['Epoch_esa'][0])))])
        Done(start_time)

        # --- get data from Trajectory Files ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict_traj = loadDictFromFile(inputFiles_trajectory[0], {})
        data_dict_traj['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch_esa'][0][i]) for i in (range(len(data_dict_traj['Epoch_esa'][0])))])
        Done(start_time)

        ############################
        # --- REDUCE THE DATASET ---
        ############################
        # Reduce Epoch to only the interested alfven signature range
        targetTimes = [pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, epochRange[0][0], epochRange[0][1], epochRange[0][2], 000000)),
                       pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, epochRange[1][0], epochRange[1][1], epochRange[1][2], 000000))]
        limitIndexes = [np.abs(data_dict['Epoch_esa'][0] - targetTimes[0]).argmin(), np.abs(data_dict['Epoch_esa'][0] - targetTimes[1]).argmin()]
        Epoch_reduced = np.array([data_dict['Epoch_esa'][0][i] for i in range(limitIndexes[0], limitIndexes[1])])  # give reduced epoch in datetimes

        # Reduce esaData to only one pitch angle
        esaData = copy.deepcopy(data_dict[wInstr[1]][0][limitIndexes[0]:limitIndexes[1]])

        # pull out other variables
        Energy = data_dict['Energy'][0]
        Pitch = data_dict['Pitch_Angle'][0]

        ############################
        # --- REMOVED INVERTED V ---
        ############################
        # desciption: the goal is to get an isolated dataset of just the alfven signatures, so we must remove the inverted V.
        # To remove the inverted V singatures from the alfven region, we need to follow this process
        # [1] Apply a counts mask of <10 counts to all the data in the region. Helps remove jitter/noise
        # [2] Create another mask which weights the contribution of higher pitch angles more than
        # lower pitch angles. Only allow energies >85.91eV into this new mask
        # [3] find some way to remove the stray contributions of alfven signatures. Perhaps zero them out


        if removeInvertedV:

            prgMsg('Removing Inverted V')

            # get the data
            esaDataCopy = copy.deepcopy(np.array([esaData[i] for i in range(len(esaData))]))
            esaDataMask = copy.deepcopy(np.array([esaData[i] for i in range(len(esaData))]))

            # --- apply a simple counts mask ---
            maskVal = 8

            # create the mask
            ranges = [range(len(esaDataCopy)), range(len(esaDataCopy[0]))]
            for tme, engy in itertools.product(*ranges):
                if esaDataMask[tme][engy] - maskVal < 0:
                    esaDataMask[tme][engy] = 0
                else:
                    esaDataMask[tme][engy] = (esaDataMask[tme][engy] - maskVal)


            # apply a nearest neighbor filter
            # description: looks at each datapoint and finds points that have 0's on either side of a value and sets that value to 0
            for tme in range(1,len(esaDataCopy)-1):
                for engy in range(len(esaDataCopy[0])):
                    if esaDataMask[tme+1][engy] == 0 and esaDataMask[tme-1][engy] == 0:
                        esaDataMask[tme][engy] = 0

            esaDataMask = np.array(esaDataMask)

            # plot the mask
            fig, ax = plt.subplots()
            ax.set_ylim(yLimits[0], yLimits[1])
            cmap = ax.pcolormesh(Epoch_reduced, Energy, esaDataMask.transpose(), cmap='turbo',vmin=vExtremes[0],vmax=vExtremes[1])
            plt.colorbar(cmap)
            plt.show()

            # apply the mask
            for tme, ptch in itertools.product(*ranges):
                if (esaDataCopy[tme][ptch] - esaDataMask[tme][ptch]) < 0:
                    esaDataCopy[tme][ptch] = 0
                else:
                    esaDataCopy[tme][ptch] = esaDataCopy[tme][ptch] - esaDataMask[tme][ptch]



            esaData = esaDataCopy
            Done(start_time)


        else:
            esaData = np.array([esaData[i][wPitch][:] for i in range(len(esaData))])







        #######################
        # --- PLOT THE DATA ---
        #######################
        if plotAlfvenRegion:

            prgMsg('Plotting Alfven Signature Region')

            fig, ax = plt.subplots()

            # colormap
            cmap = ax.pcolormesh(Epoch_reduced, Energy, esaData.transpose(), vmin=vExtremes[0], vmax=vExtremes[1], cmap="turbo")

            ax.set_ylim(yLimits[0], yLimits[1])

            cbar = plt.colorbar(cmap)
            cbar.minorticks_on()
            cbar.set_label('Counts')

            # labels
            ax.set_ylabel('Energy [eV]')
            ax.set_xlabel('Epoch')
            ax.set_title(f'ACESII {rocketAttrs.rocketID[wRocket - 4]} \n {wInstr[1].upper()} \n Pitch {Pitch[wPitch]}$^\circ$')

            # yticks
            spacing = 20
            Epoch_tick_locations = Epoch_reduced[::spacing]
            Epoch_datetime = np.array([pycdf.lib.tt2000_to_datetime(Epoch_tick_locations[i]).strftime("%M:%S") for i in range(len(Epoch_tick_locations))])
            ax.set_xticks(Epoch_tick_locations) #WHERE the ticks are place
            ax.set_xticklabels(Epoch_datetime)

            #xticks
            spaceing = 100
            ax.set_yticks([i for i in range(0,yLimits[1],spaceing)])










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
if wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        AlfvenSignatureAnalysis(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            AlfvenSignatureAnalysis(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            AlfvenSignatureAnalysis(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)