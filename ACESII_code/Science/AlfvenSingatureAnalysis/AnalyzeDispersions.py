# --- AnalyzeDispersions.py ---
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
justPrintFileNames = False

wRocket = 4

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_trajectory = 'trajectories'
outputPath_modifier = 'science\AlfvenSignatureAnalysis' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder


wPitch = 2

# Plot Single Dispersion Feature
plotKeyDispersions = True
wDispersions = [] # [] -> plot them all, [#,#,#,...] plot specific ones
plotAdditionalDispersions = False


# Subtracting Inverted V
removeInvertedV = False # Tries to remove InvertedV from the data
plotFilteredMask = False # plot the weighted mask used on all the pitch angles. requires removeInvertedV = True


# PLOTTING THE DATA
plotAlfvenRegion = False # plots the specific wPitch esa data
epochRange = [[17, 24, 56], [17, 25, 9]] # UTC Ranges for plot Epoch. Format: [hour, minute,second]
# epochRange = [[17, 25, 1], [17, 25, 9]] # UTC Ranges for plot Epoch. Format: [hour, minute,second]
vExtremes = [0, 40] #colorbar limits
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
from scipy.ndimage import gaussian_filter
from matplotlib import figure
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg,outputCDFdata,L1_TRICE_Quick,L2_TRICE_Quick,loadDictFromFile, q0,m_e
from glob import glob
from os.path import getsize
from scipy.optimize import curve_fit


setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


def AlfvenSignatureAnalysis(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer,wDispersions):

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

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            descriptiorNam = ['EEPAA', 'LEESA', 'IEPAA', 'Langmuir_Probe']
            wInstr = [index, instr, descriptiorNam[index]]

    fileoutName = f'ACESII_{rocketID}_{wInstr[1]}_InvertedV_removed.cdf'

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

        # Reduce esaData to only the data within the alfven region
        esaData = copy.deepcopy(data_dict[wInstr[1]][0][limitIndexes[0]:limitIndexes[1]])
        esaDataRawForPlotting = copy.deepcopy(data_dict[wInstr[1]][0][limitIndexes[0]:limitIndexes[1]])

        # pull out other variables
        Energy = data_dict['Energy'][0]
        Pitch = data_dict['Pitch_Angle'][0]

        if plotKeyDispersions:
            from ACESII_code.Science.AlfvenSingatureAnalysis.dispersionAttributes import dispersionAttributes

            if wDispersions == []:
                wDispersions = [i for i in range(len(dispersionAttributes.keyDispersionTimes))]

            for i in range(len(wDispersions)):

                # isolate a single dispersion trace
                targetTimes = [ pycdf.lib.datetime_to_tt2000(dispersionAttributes.keyDispersionTimes[i][0]),pycdf.lib.datetime_to_tt2000(dispersionAttributes.keyDispersionTimes[i][0])]
                limitIndexes = [np.abs(data_dict['Epoch_esa'][0] - targetTimes[0]).argmin(), np.abs(data_dict['Epoch_esa'][0] - targetTimes[1]).argmin()]
                EpochOneDis = np.array([data_dict['Epoch_esa'][0][i] for i in range(limitIndexes[0], limitIndexes[1])])
                esaDataOneDis = np.array(data_dict[wInstr[1]][0][limitIndexes[0]:limitIndexes[1]])
                esaDataOneDis = np.array([esaDataOneDis[tme][wPitch][:] for tme in range(len(esaDataOneDis))])

                print(len(esaDataOneDis),len(EpochOneDis),targetTimes,limitIndexes)

                # mask the data by maskval
                maskval = 1
                for i in range(len(esaDataOneDis)):
                    for j in range(len(esaDataOneDis[i])):
                        if esaDataOneDis[i][j] - maskval > 0:
                            esaDataOneDis[i][j] = esaDataOneDis[i][j] - maskval
                        else:
                            esaDataOneDis[i][j] = 0

                # create an X-Y dataset for times and energies
                fig, ax = plt.subplots()
                xData = [] # time
                yData = [] # energy
                for tme in range(len(esaDataOneDis)):
                    for engy in range(len(esaDataOneDis[0])):
                        for j in range(int(esaDataOneDis[tme][engy])):
                            xData.append(EpochOneDis[tme])
                            yData.append(Energy[engy])

                # sort the data based on the Epoch (xdata)
                xData, yData = zip(*sorted(zip(xData, yData)))
                xData, yData = np.array(xData), np.array(yData)

                # convert x data to "time since launch, in seconds"
                launchTimeTT2000 = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 000))
                xData = (xData - launchTimeTT2000)/1E9
                yData = np.array([np.sqrt(2*y*q0/(m_e)) for y in yData])

                ax.scatter(xData, yData)
                spacing = 40
                ax.set_ylabel('Energy [eV]')
                ax.set_xlabel('Epoch [seconds]')
                ax.set_title(f'Dispersion No: {wDispersions[i]}')
                # play around with the fit function
                Emin = yData.min()
                Vmin = np.sqrt(2*Emin*q0/m_e)
                tmax = xData.min()
                #
                #
                # # def fitFunc(x, d):
                # #     gamma = ((tmax - x)) * np.cos(np.radians(Pitch[wPitch]) / d)
                # #     y = Emin / ((gamma * np.sqrt(2 * Emin * q0 / m_e) - 1) ** 2)
                # #     return y
                #
                # def fitFunc(x, d):
                #     y = (d/np.cos(np.radians(Pitch[wPitch]))) * (1/tmax - 1/x) + Vmin
                #     return y
                #
                # xDataFit = np.linspace(xData.min(),xData.max(),1000)
                # params,cov = curve_fit(fitFunc,xData,yData)
                # yDataFit = [fitFunc(x, params[0]) for x in xDataFit]
                # plt.plot(xDataFit, yDataFit,color='red')


                plt.show()






        ############################
        # --- REMOVED INVERTED V ---
        ############################
        # desciption: the goal is to get an isolated dataset of just the alfven signatures, so we must remove the inverted V.
        # To remove the inverted V singatures from the alfven region, we need to follow this process
        # [1] Search for all points below a maskval, store this data to be used later
        # [2] Subtract the same maskval from the data
        # [3] take the result and do a gaussian filter
        # [4] remove the low energy contribution from the filtered data
        # [5] subtract another tuning val form the filtered data to ensure alfven contribution isn't strong
        # [6] multiply the mask by some value (maybe do a minimization technqiue?)
        # [7] subtract mask from raw data
        # [8] "re-build" data by looking where counts still exist and go back to raw data and replace new values wherever there are some


        if removeInvertedV:

            # get the data
            esaDataCopy = copy.deepcopy(np.array(esaData))
            esaDataReBuild = copy.deepcopy(np.array(esaData))
            esaDataMask = copy.deepcopy(np.array(esaData))

            # --- --- --- --- --- --- --- --- --- --
            # --- [1] apply a simple counts mask ---
            # --- --- --- --- --- --- --- --- --- --
            prgMsg('Applying Simple Counts Mask')
            # Weights of the various pitch pads
            maskVal1 = 5  # for the simple mask value, remove this many counts from the mask

            # Apply a simple overall counts mask using maskVal
            ranges = [range(len(Epoch_reduced)), range(len(Pitch)), range(len(Energy))]
            for tme, ptch, engy in itertools.product(*ranges):
                if (esaDataMask[tme][ptch][engy] - maskVal1 > 0):
                    esaDataMask[tme][ptch][engy] = (esaDataMask[tme][ptch][engy] - maskVal1)
                else:
                    esaDataMask[tme][ptch][engy] = 0
            Done(start_time)

            # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
            # --- [2] Do a guassian filter on the simple masked data ---
            # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
            # NOTE: the data format changes here. Now filteredData is [ptchSlice0, ptchSlice1, ... ]
            Gsigma = 0.5  # std dev to use on the gaussian filter
            filteredData = []
            NoOfGFilters = 2
            for i in range(NoOfGFilters):
                if i == 0:
                    filteredData =np.array([gaussian_filter(np.array([esaDataMask[tme][ptch][:] for tme in ranges[0]]),sigma=Gsigma) for ptch in ranges[1]])
                else:
                    filteredData = np.array([gaussian_filter(filteredData[ptch], sigma=Gsigma) for ptch in ranges[1]])

            # --- --- --- --- --- --- --- --- --- --- --
            # --- [3] Remove low energy contribution ---
            # --- --- --- --- --- --- --- --- --- --- --
            removeEngy = 34  # corresponds to energy 100.44eV. This energy should be INCLUSIVE. Nominally 31
            for ptch in ranges[1]:
                for tme in ranges[0]:
                    for engy in ranges[2]:
                        if engy > removeEngy:
                            filteredData[ptch][tme][engy] = 0

            # --- --- --- --- --- --- --- --- --- --- --- --- --- --
            # --- [4] Subtract another maskVal from filteredData ---
            # --- --- --- --- --- --- --- --- --- --- --- --- --- --
            maskVal2 = 0
            for ptch in ranges[1]:
                for tme in ranges[0]:
                    for engy in ranges[2]:
                        if filteredData[ptch][tme][engy] - maskVal2 > 0:
                            filteredData[ptch][tme][engy] = filteredData[ptch][tme][engy] - maskVal2
                        else:
                            filteredData[ptch][tme][engy] = 0

            # --- --- --- --- --- --- --- --- --- --- --
            # --- [5] Add multiplier to filteredMask ---
            # --- --- --- --- --- --- --- --- --- --- --
            multipler = 1
            for ptch in ranges[1]:
                for tme in ranges[0]:
                    for engy in ranges[2]:
                        filteredData[ptch][tme][engy] = int(filteredData[ptch][tme][engy]*multipler)


            # plot the weighted mask
            if plotFilteredMask:
                spacing = 20

                fig, ax = plt.subplots()
                plt.title(f'Weighted Mask\n'
                          f'maskVal1: {maskVal1},  Energy Cutoff {Energy[removeEngy]} eV\n'
                          f'Guassian filter Sigma: {Gsigma},  maskval2: {maskVal2}'
                          )
                Epoch_tick_locations = Epoch_reduced[::spacing]
                Epoch_datetime = np.array([pycdf.lib.tt2000_to_datetime(Epoch_tick_locations[i]).strftime("%M:%S") for i in range(len(Epoch_tick_locations))])
                ax.set_ylim(yLimits[0], yLimits[1])
                ax.set_xticks(Epoch_tick_locations)  # WHERE the ticks are place
                ax.set_xticklabels(Epoch_datetime)
                cmap = ax.pcolormesh(Epoch_reduced, Energy, filteredData[wPitch].T, cmap='turbo', vmin=vExtremes[0], vmax=vExtremes[1])
                plt.colorbar(cmap)
                plt.show()

            # --- --- --- --- --- --- --- --- --- ---
            # --- [6] Apply the filteredData Mask ---
            # --- --- --- --- --- --- --- --- --- ---
            for ptch in ranges[1]:
                for tme in ranges[0]:
                    for engy in ranges[2]:
                        dif = esaDataCopy[tme][ptch][engy] - filteredData[ptch][tme][engy]
                        if dif > 0:
                            esaDataCopy[tme][ptch][engy] = dif
                        else:
                            esaDataCopy[tme][ptch][engy] = 0

            # --- --- --- --- --- --- --- --
            # --- [7] "re-build" dataset ---
            # --- --- --- --- --- --- --- --
            rebuildThreshold = 1000
            for tme, ptch, engy in itertools.product(*ranges):
                if np.abs(esaDataCopy[tme][ptch][engy] - esaDataReBuild[tme][ptch][engy]) > 0 and esaDataCopy[tme][ptch][engy] > rebuildThreshold:
                    esaDataCopy[tme][ptch][engy] = int(esaDataReBuild[tme][ptch][engy])

            esaData_removedV = esaDataCopy
            Done(start_time)

        else:
            esaData_removedV = np.array([esaData[i][wPitch][:] for i in range(len(esaData))])



        #######################
        # --- PLOT THE DATA ---
        #######################
        if plotAlfvenRegion:

            prgMsg('Plotting Alfven Signature Region')

            fig, ax = plt.subplots(2)

            fig.suptitle(f'ACESII {rocketAttrs.rocketID[wRocket - 4]} \n {wInstr[1].upper()} \n Pitch {Pitch[wPitch]}$^\circ$')

            for i in range(2):
                if i == 0: # masked data
                    plotThisData = np.array([esaData_removedV[i][wPitch][:] for i in range(len(Epoch_reduced))])

                elif i == 1: # raw data
                    plotThisData = np.array([esaDataRawForPlotting[i][wPitch][:] for i in range(len(Epoch_reduced))])

                cmap = ax[i].pcolormesh(Epoch_reduced, Energy, plotThisData.transpose(), vmin=vExtremes[0], vmax=vExtremes[1], cmap="turbo")
                cbar = plt.colorbar(cmap)
                cbar.minorticks_on()
                cbar.set_label('Counts')
                ax[i].set_ylim(yLimits[0], yLimits[1])

                # labels
                ax[i].set_ylabel('Energy [eV]')
                ax[i].set_xlabel('Epoch')

                # xticks
                spacing = 20
                Epoch_tick_locations = Epoch_reduced[::spacing]
                Epoch_datetime = np.array([pycdf.lib.tt2000_to_datetime(Epoch_tick_locations[i]).strftime("%M:%S") for i in range(len(Epoch_tick_locations))])
                ax[i].set_xticks(Epoch_tick_locations) #WHERE the ticks are place
                ax[i].set_xticklabels(Epoch_datetime)

                # yticks
                spaceing = 100
                ax[i].set_yticks([i for i in range(0, yLimits[1], spaceing)])

            plt.show()
            Done(start_time)






        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        if outputData:
            prgMsg('Creating output file')

            data_dict = {**data_dict, **{f'Epoch_reduced': [Epoch_reduced, {'LABLAXIS': 'Epoch_reduced',
                                                          'DEPEND_0': 'Epoch_reduced',
                                                          'DEPEND_1': None,
                                                          'DEPEND_2': None,
                                                          'FILLVAL': rocketAttrs.epoch_fillVal,
                                                          'FORMAT': 'E12.2', 'UNITS': 'ns',
                                                          'VALIDMIN': Epoch_reduced.min(),
                                                          'VALIDMAX': Epoch_reduced.max(),
                                                          'VAR_TYPE': 'support_data',
                                                          'MONOTON': 'INCREASE',
                                                          'TIME_BASE': 'J2000',
                                                          'TIME_SCALE': 'Terrestrial Time',
                                                          'REFERENCE_POSITION': 'Rotating Earth Geoid',
                                                          'SCALETYP': 'linear'}]}}

            data_dict = {**data_dict, **{f'{wInstr[1]}_no_invertedV': [esaData, {'LABLAXIS': f'{wInstr[1]}',
                                                     'DEPEND_0': 'Epoch_reduced',
                                                     'DEPEND_1': 'Pitch_Angle',
                                                     'DEPEND_2': 'Energy',
                                                     'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                     'UNITS': 'counts',
                                                     'VALIDMIN': esaData.min(),
                                                     'VALIDMAX': esaData.max(),
                                                     'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}


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
        AlfvenSignatureAnalysis(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer,wDispersions)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            AlfvenSignatureAnalysis(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer,wDispersions)
    else:
        for filesNo in wFiles:
            AlfvenSignatureAnalysis(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer,wDispersions)