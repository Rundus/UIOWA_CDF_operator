# --- LP_ChiSquare_Temp&Density.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Takes in the L2 Langmuir Data and performs the Chi Square analysis to get the Temperature and
# Density from the characteristic curves


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools
import math
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
inputPath_modifier = 'L2'  # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\Langmuir'  # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# --- Fitting Toggles ---
unitConv = 1E9  # converts from A to nA
wSweeps = [100 + i for i in range(10)]  # [] --> all sweeps, [#1,#2,...] specific sweeps

# bad data: High Flyer 36
# auroral case: 260
# calm case:
# Exponential case: 40, 100, 70,90, 150
# Linear Cases: 130



# --- Data Plotting ---
plotEveryFitPerSweep = True

# --- OutputData ---
outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os, scipy
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg, L2_ACES_Quick, q0, m_e, kB, IonMasses
from glob import glob
from os.path import getsize
from decimal import Decimal
from matplotlib import pyplot as plt
setupPYCDF()
from spacepy import pycdf
from tqdm import tqdm
from LP_gridSearch_toggles import Ti,rounding, Vsp_range,voltageStartPoint, transFitParameters,satFitParameters
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
pycdf.lib.set_backward(False)



##############################
# --- FITTED CHI FUNCTIONS ---
##############################
rocketAttrs, b, c = ACES_mission_dicts()
e_coefficient = (((q0 ** 3) * (rocketAttrs.LP_probe_areas[0][0] ** (2))) / (8 * np.pi * m_e)) ** (1 / 2)
i_coefficient = ((Ti * (q0 ** 3) * (rocketAttrs.LP_probe_areas[0][0] ** (2))) / (8 * np.pi * IonMasses[0])) ** (1 / 2)

def transitionFunc(x, a0, a1, a2):
    y = (e_coefficient) * a0 * np.exp((x - a1) / a2)
    return y
def saturationFunc(x, a0, a1, a2):
    y = ((e_coefficient) * a0 - ((a0 / (a2 ** (1 / 2))) * (i_coefficient) * np.exp(-1 * ((x - a1) / Ti))))
    return y

#######################
# --- MAIN FUNCTION ---
#######################
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

    input_names_searchable = [ifile.replace(inputPath_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace(outputPath_modifier.lower() + '_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\',
                                              '')
    fileoutName = dataFile_name.replace('langmuir_', 'langmuir_Temp&Density_')

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made L2: {:3s} '.format(i, input_names_searchable[i],
                                                                       round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict = {}
        with pycdf.CDF(inputFiles[wFile]) as inputDataFile:
            for key, val in inputDataFile.items():
                data_dict = {**data_dict, **{key: [inputDataFile[key][...], {key: val for key, val in inputDataFile[key].attrs.items()}]}}

        data_dict['Epoch_swept_Current'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_swept_Current'][0][i]) for i in (range(len(data_dict['Epoch_swept_Current'][0])))])

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

            if np.abs((data_dict['Epoch_swept_Current'][0][i + 1] / 100000 - data_dict['Epoch_swept_Current'][0][i] / 100000)) >= indvEpochThresh / 100000:

                if counter == 0:
                    sweepIndices.append([0, i])
                    counter += 1
                else:
                    sweepIndices.append([sweepIndices[-1][1] + 1, i])

        # include the last sweep
        sweepIndices.append([sweepIndices[-1][1] + 1, len(data_dict['Epoch_swept_Current'][0])])

        # Add a QA set where we remove anyplace where sweepIndicies[i][0] == sweepIndicies[i][1]
        removeThese = []
        for i in range(len(sweepIndices)):
            if sweepIndices[i][0] == sweepIndices[i][1]:
                removeThese.append(i)

        for thing in removeThese:
            del sweepIndices[thing]

        # Find the index corresponding to the break in upleg/downleg voltage for the sweeps
        breakIndices = []
        qualityCounter = 0
        for i in range(len(sweepIndices)):
            start = sweepIndices[i][0]
            end = sweepIndices[i][1]
            sweeps = np.array(data_dict['swept_Current'][0][start:end])
            sweepMAX, sweepMIN = sweeps.max(), sweeps.min()
            threshold = (sweepMAX - sweepMIN) * 0.75

            for j in range(len(sweeps) - 1):
                if np.abs(sweeps[j + 1] - sweeps[j]) >= threshold:
                    breakIndices.append(j + start)

                    if start < (j + start) < end:
                        qualityCounter += 1

        Done(start_time)

        print(len(sweepIndices), qualityCounter)

        sweepIndices, breakIndices = np.array(sweepIndices), np.array(breakIndices)

        # --- store the sweep data into a single data variables filled with individual sweeps---
        prgMsg('Reorganizing Data')
        sweepsCurrent, sweepsVoltage, sweepsCurrent_epoch = [], [], []

        for sweep in range(len(sweepIndices)):
            start = sweepIndices[sweep][0]
            end = sweepIndices[sweep][1]
            breakPoint = breakIndices[sweep]

            # Get the data and sort it
            xData1 = np.array(data_dict['swept_Voltage'][0][start:breakPoint])
            yData1 = np.array(data_dict['swept_Current'][0][start:breakPoint])
            yData1 = np.array([x for _, x in sorted(zip(xData1, yData1))])
            xData1 = np.array(sorted(xData1))

            xData2 = data_dict['swept_Voltage'][0][breakPoint:end]
            yData2 = np.array(data_dict['swept_Current'][0][breakPoint:end])

            yData2 = np.array([x for _, x in sorted(zip(xData2, yData2))])
            xData2 = np.array(sorted(xData2))

            # append data to lists
            sweepsVoltage.append(xData1)
            sweepsVoltage.append(xData2)
            sweepsCurrent_epoch.append(data_dict['Epoch_swept_Current'][0][start])
            sweepsCurrent.append(yData1)
            sweepsCurrent.append(yData2)
            sweepsCurrent_epoch.append(data_dict['Epoch_swept_Current'][0][breakPoint])

        sweepsVoltage = np.array(sweepsVoltage, dtype='object')
        sweepsCurrent = np.array(sweepsCurrent, dtype='object')
        sweepsCurrent_epoch = np.array(sweepsCurrent_epoch, dtype='object')
        Done(start_time)

        ###########################
        # --- APPLY GRID SEARCH ---
        ###########################

        prgMsg('Performing Grid Search')
        print('\n')

        # data that will be outputed
        sweeps_Te = []
        sweeps_Vsp = []
        sweeps_n0 = []
        sweeps_Epoch = []

        error_ProbeVoltage = data_dict['errorInProbeVoltage'][0]
        error_ProbeCurrent = data_dict['errorInCaldCurrent'][0]/unitConv

        # analyze the specific sweeps in wSweeps or do all of them
        theseSweeps = [i for i in range(len(sweepsCurrent))] if wSweeps == [] else wSweeps

        # loop through all the sweeps
        for sweepNo in tqdm(theseSweeps):

            # for each sweep, store the results of the grid search. format: [transition,saturation]
            gridSearchResults = {
                'n0': [0, 0],
                'Te': [0, 0],
                'Vsp': [0, 0],
                'timeStamp': pycdf.lib.tt2000_to_datetime(sweepsCurrent_epoch[sweepNo]).strftime("%H:%M:%S.%f"),
                'Epoch':sweepsCurrent_epoch[sweepNo],
                'ChiSquare': [10000000000, 10000000000],
                'transData' : [0,0],
                'satData': [0,0],
                'transParams' : [0],
                'satParams' : [0],
                'Ii0' : 0,
                'wSweep' : sweepNo
            }

            # reduce the x and y Data to only consider past Vprobe > 0V
            zeroPoint = np.abs(sweepsVoltage[sweepNo] - voltageStartPoint).argmin()
            xData, yData = sweepsVoltage[sweepNo][zeroPoint:], np.array(sweepsCurrent[sweepNo][zeroPoint:])/unitConv

            #####################
            # --- GRID SEARCH ---
            #####################

            bestChiChecker = 10000000000000000000000000000000

            # --- LOOP THROUGH Vsp VALUES ---
            for i in range(len(Vsp_range)):

                # Break data into Transition and Saturation
                breakHere = np.abs(xData - Vsp_range[i]).argmin()
                xDataTrans, yDataTrans = np.array(xData[:breakHere]), np.array(yData[:breakHere])
                xDataSat, yDataSat = np.array(xData[breakHere:]), np.array(yData[breakHere:])

                # Ion saturation Current
                Ii0 = min(yDataTrans)
                yDataTrans = [y - Ii0 for y in yDataTrans]  # subtract Ii0 from transition data

                # Fit the transition region and calculate ChiSquare
                transParams, tCov = scipy.optimize.curve_fit(transitionFunc, xDataTrans, yDataTrans, **transFitParameters)
                nu = len(xDataTrans) - 3

                ChiSquareTrans = sum(
                    [
                        (yDataTrans[i] - transitionFunc(xDataTrans[i], transParams[0], transParams[1], transParams[2])) ** 2 / (error_ProbeCurrent**2)
                        for i in range(len(xDataTrans))
                    ]) / nu

                # Fit the saturation region and calculate ChiSquare
                satParams, sCov = scipy.optimize.curve_fit(saturationFunc, xDataSat, yDataSat, **satFitParameters)

                nu = len(xDataTrans) - 3
                ChiSquareSat = sum([(yDataSat[i] - saturationFunc(xDataSat[i], satParams[0], satParams[1], satParams[2])) ** 2 / (error_ProbeCurrent ** 2) for i in range(len(xDataSat))])/nu

                # Reset transition data back
                yDataTrans = [y + Ii0 for y in yDataTrans]  # subtract Ii0 from transition data


                if (ChiSquareSat + ChiSquareTrans) <  bestChiChecker:
                    gridSearchResults['n0'] = [transParams[0] / np.sqrt(transParams[2]), satParams[0] / np.sqrt(satParams[2])]
                    gridSearchResults['Te'] = [transParams[2], satParams[2]]
                    gridSearchResults['Vsp'] = [transParams[1], satParams[1]]
                    gridSearchResults['ChiSquare'] = [ChiSquareTrans, ChiSquareSat]
                    gridSearchResults['transData'] = [xDataTrans, yDataTrans]
                    gridSearchResults['satData'] = [xDataSat, yDataSat]
                    gridSearchResults['transParams'] = transParams
                    gridSearchResults['satParams'] = satParams
                    gridSearchResults['Ii0'] = Ii0

                    # update the best Chi value
                    bestChiChecker = ChiSquareSat + ChiSquareTrans



            # --- OUTPUT THE DATA TO THE FINAL STORAGE OBJECT ---

            if plotEveryFitPerSweep:

                # --- plot test data ---
                xDataTrans = gridSearchResults['transData'][0]
                yDataTrans = gridSearchResults['transData'][1]
                xDataSat = gridSearchResults['satData'][0]
                yDataSat = gridSearchResults['satData'][1]
                transParams = gridSearchResults['transParams']
                satParams = gridSearchResults['satParams']
                Ii0 = gridSearchResults['Ii0']
                ChiSquareTrans = gridSearchResults['ChiSquare'][0]
                ChiSquareSat = gridSearchResults['ChiSquare'][1]

                xDataTest_trans = np.linspace(min(xDataTrans), max(xDataTrans))
                xDataTest_sat = np.linspace(min(xDataSat), max(xDataSat))
                yDataTest_trans = [transitionFunc(x, transParams[0], transParams[1], transParams[2]) + Ii0 for x in xDataTest_trans]
                yDataTest_sat = [saturationFunc(x, satParams[0], satParams[1], satParams[2]) for x in xDataTest_sat]

                fig, ax = plt.subplots(2)
                ax[0].scatter(xDataTrans, yDataTrans)
                ax[0].plot(xDataTest_trans, yDataTest_trans, color='red')
                ax[1].scatter(xDataSat, yDataSat)
                ax[1].plot(xDataTest_sat, yDataTest_sat)
                ax[0].set_title(r'$\chi_{\nu} ^2$ ' + f'{ChiSquareTrans}')
                ax[1].set_title(r'$\chi_{\nu} ^2$ ' + f'{ChiSquareSat}')
                ax[0].set_xlabel('Probe Voltage [V]')
                ax[0].set_ylabel('Probe Current [A]')
                ax[1].set_ylabel('Probe Current [A]')
                ax[1].set_xlabel('Probe Voltage [V]')
                plt.suptitle(str(gridSearchResults['wSweep']) + '. ' + gridSearchResults['timeStamp'])
                plt.tight_layout()
                a0 = [transParams[0], satParams[0]]
                a1 = [transParams[1], satParams[1]]
                a2 = [transParams[2], satParams[2]]
                n_0 = [a0[0] / np.sqrt(a2[0]), a0[1] / np.sqrt(a2[1])]

                for i in range(2):
                    ax[i].legend([
                        '$a_{0} = n_{0} T_{e[eV]}^{1/2}$: ' + f'{a0[i]} \n' +
                        '$a_{1} = V_{sp}$: ' + f'{a1[i]}' + ' V' + '\n' +
                        '$a_{2} = T_{e[eV]}$: ' + f'{a2[i]}' + ' eV' + '\n' +
                        '$T_{i[eV]}$= ' + f'{Ti}' + ' eV' + '\n' +
                        '$n_{0} =$' + f'{"{:.2E}".format(Decimal(str(round((10 ** (-6)) * n_0[i], rounding))))} ' + ' $cm^{-3}$' + ''])
                plt.show()

            # determine which data should be output'd
            outputThis = 0 if (gridSearchResults['ChiSquare'][0] > gridSearchResults['ChiSquare'][1]) else 1
            outputThis = 1

            # output the data to the final data storage variables
            sweeps_Epoch.append(gridSearchResults['Epoch'])
            sweeps_Vsp.append(gridSearchResults['Vsp'][outputThis])
            sweeps_Te.append(gridSearchResults['Te'][outputThis])
            sweeps_n0.append(gridSearchResults['n0'][outputThis])


        #####################
        # --- OUTPUT DATA ---
        #####################
        if outputData:
            # --- --- --- --- --- --- ---
            # --- WRITE OUT THE DATA ---
            # --- --- --- --- --- --- ---
            prgMsg('Creating output file')


            # remove unwanted data from output file
            sweeps_Epoch = np.array(sweeps_Epoch)
            sweeps_Vsp = np.array(sweeps_Vsp)
            sweeps_n0 = np.array(sweeps_n0)
            sweeps_Te = np.array(sweeps_Te)

            data_dict = {**data_dict, **{'sweeps_Epoch': [sweeps_Epoch, {'LABLAXIS': 'sweeps_Epoch',
                                                                                      'DEPEND_0': 'sweeps_Epoch',
                                                                                      'DEPEND_1': None,
                                                                                      'DEPEND_2': None,
                                                                                      'FILLVAL': -9223372036854775808,
                                                                                      'FORMAT': 'E12.2', 'UNITS': 'ns',
                                                                                      'VALIDMIN': sweeps_Epoch.min(),
                                                                                      'VALIDMAX': sweeps_Epoch.max(),
                                                                                      'VAR_TYPE': 'support_data',
                                                                                      'MONOTON': 'INCREASE',
                                                                                      'TIME_BASE': 'J2000',
                                                                                      'TIME_SCALE': 'Terrestrial Time',
                                                                                      'REFERENCE_POSITION': 'Rotating Earth Geoid',
                                                                                      'SCALETYP': 'linear'}]}}

            data_dict = {**data_dict, **{'Vsp': [sweeps_Vsp, {'LABLAXIS': 'Vsp',
                                                                               'DEPEND_0': 'sweeps_Epoch',
                                                                               'DEPEND_1': None,
                                                                               'DEPEND_2': None,
                                                                               'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                                               'UNITS': 'Volts',
                                                                               'VALIDMIN': sweeps_Vsp.min(),
                                                                               'VALIDMAX': sweeps_Vsp.max(),
                                                                               'VAR_TYPE': 'data',
                                                                               'SCALETYP': 'linear'}]}}

            data_dict = {**data_dict, **{'n0': [sweeps_n0, {'LABLAXIS': 'n0',
                                                                     'DEPEND_0': 'sweeps_Epoch',
                                                                     'DEPEND_1': None,
                                                                     'DEPEND_2': None,
                                                                     'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                                     'UNITS': 'cm!U-3',
                                                                     'VALIDMIN': sweeps_n0.min(),
                                                                     'VALIDMAX': sweeps_n0.max(),
                                                                     'VAR_TYPE': 'data',
                                                                     'SCALETYP': 'linear'}]}}

            data_dict = {**data_dict, **{'Te': [sweeps_Te, {'LABLAXIS': 'Te',
                                                                    'DEPEND_0': 'sweeps_Epoch',
                                                                    'DEPEND_1': None,
                                                                    'DEPEND_2': None,
                                                                    'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                                    'UNITS': 'eV',
                                                                    'VALIDMIN': sweeps_Te.min(),
                                                                    'VALIDMAX': sweeps_Te.max(),
                                                                    'VAR_TYPE': 'data',
                                                                    'SCALETYP': 'linear'}]}}

            del data_dict['swept_Current'],data_dict['swept_Voltage'],data_dict['errorInProbeVoltage'],data_dict['errorInCaldCurrent'],data_dict['Epoch_step'],data_dict['step_Voltage']

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
                    if 'Epoch' in varKey:  # epoch data
                        L2File.new(varKey, data=varVal[0], type=33)
                    else:  # other data
                        L2File.new(varKey, data=varVal[0], type=pycdf.const.CDF_REAL8)

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
elif wRocket == 1:  # ACES II Integration Low
    rocketFolderPath = Integration_data_folder
    wflyer = 1
elif wRocket == 2:  # TRICE II High
    rocketFolderPath = TRICE_data_folder
    wflyer = 0
elif wRocket == 3:  # TRICE II Low
    rocketFolderPath = TRICE_data_folder
    wflyer = 1
elif wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5:  # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        main(wRocket, 0, rocketFolderPath, justPrintFileNames, wflyer)
    else:
        main(wRocket, 0, rocketFolderPath, justPrintFileNames, wflyer)
