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
inputPath_modifier = 'L2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\Langmuir' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# --- Fitting Toggles ---
unitConv = 10**(9) # converts from A to nA

wSweeps = [] # [] --> all sweeps, [#1,#2,...] specific sweeps


# bad data: High Flyer 36
# auroral case: 260
# calm case:
# Exponential case: 40, 100, 70,90, 150
# Linear Cases: 130

rounding = 9

# --- Data Plotting ---
plotBestFit = True
plotFitParamsOverTime = False

# --- OutputData ---
outputData = False



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os, scipy
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg, L2_ACES_Quick,q0,m_e,kB,IonMasses
from glob import glob
from os.path import getsize
from matplotlib import pyplot as plt
setupPYCDF()
from spacepy import pycdf
from tqdm import tqdm
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
pycdf.lib.set_backward(False)



#############################
# --- GRID SEARCH TOGGLES ---
#############################
Ti = 0.2
Vsp_range = np.linspace(0.8, 1.2, 10)
voltageStartPoint = 0 # Only take data > than this voltage

##############################
# --- FITTED CHI FUNCTIONS ---
##############################
rocketAttrs, b, c = ACES_mission_dicts()

e_coefficient = (((q0 ** 3) * (rocketAttrs.LP_probe_areas[0][0] ** (2))) / (8 * np.pi * m_e)) ** (1 / 2)
i_coefficient = ((Ti * (q0 ** 3) * (rocketAttrs.LP_probe_areas[0][0] ** (2))) / (8 * np.pi * IonMasses[0])) ** (1 / 2)
def transitionFunc(x, a0, a1, a2):
    y = unitConv * (e_coefficient) * a0 * np.exp((x - a1) / a2)
    return y

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

    input_names_searchable = [ifile.replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

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

        data_dict['Epoch_swept_Current'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_swept_Current'][0][i]) for i in (range(len(data_dict['Epoch_swept_Current'][0])  )   )])

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
            threshold = (sweepMAX - sweepMIN)*0.75

            for j in range(len(sweeps)-1):
                if np.abs(sweeps[j+1] - sweeps[j]) >= threshold:
                    breakIndices.append(j+start)

                    if  start < (j+start) < end:
                        qualityCounter += 1

        Done(start_time)

        print(len(sweepIndices), qualityCounter)


        sweepIndices = np.array(sweepIndices)
        breakIndices = np.array(breakIndices)

        # --- store the sweep data into a single data variables filled with individual sweeps---
        prgMsg('Reorganizing Data')
        sweepsCurrent = []
        sweepsVoltage = []
        sweepsCurrent_epoch = []

        for sweep in range(len(sweepIndices)):
            start = sweepIndices[sweep][0]
            end = sweepIndices[sweep][1]
            breakPoint = breakIndices[sweep]

            # Get the data and sort it
            xData1 = np.array(data_dict['swept_Voltage'][0][start:breakPoint])
            yData1 = np.array(data_dict['swept_Current'][0][start:breakPoint]) / unitConv
            yData1 = np.array([x*unitConv for _, x in sorted(zip(xData1, yData1))])
            xData1 = np.array(sorted(xData1))

            xData2 = data_dict['swept_Voltage'][0][breakPoint:end]
            yData2 = np.array(data_dict['swept_Current'][0][breakPoint:end]) / unitConv

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


        # analyze the specific sweeps in wSweeps or do all of them
        numDensity = [[],[]]
        electronTemp = [[],[]]
        plasmaPotential = [[],[]]
        parameterEpoch = [[],[]]
        chiSquares = [[],[]]
        theseSweeps = [i for i in range(len(sweepsCurrent))] if wSweeps == [] else wSweeps

        for sweepNo in tqdm(theseSweeps):
            chiParamsTrans = []
            chiParamsSat = []
            xData = sweepsVoltage[sweepNo]
            zeroPoint = np.abs(xData - voltageStartPoint).argmin() # only consider data past Vprobe > 0
            def saturationFunc(x, a0, a1, a2):
                y = unitConv * ((e_coefficient) * a0 - ((a0 / (a2 ** (1 / 2))) * (i_coefficient) * np.exp(-1 * ((x - a1) / Ti))))
                return y

            data = {
                'timeStamp': pycdf.lib.tt2000_to_datetime(sweepsCurrent_epoch[sweepNo]).strftime("%H:%M:%S.%f"),
                'Functions': [transitionFunc, saturationFunc],
                'xData': [[],[]],
                'yData': [[],[]],
                'fitGuess': [[2.11E7, 1.5, 0.08296951], [2.11E7, 1.5, 0.2]],
                'fitBounds': [([0, 0, 0], [np.Inf, 3, np.Inf]),
                              ([0, 0, 0], [np.Inf, 3, np.Inf])],
                'fitParams': [[], []],
                'ChiSquare': [[], []]
            }

            # --- LOOP THROUGH chi^2 VALUES ---
            for i in range(len((Vsp_range))):
                # Break data into Transition and Saturation
                breakHere = np.abs(xData - Vsp_range[i]).argmin()
                data['xData'][0].append(sweepsVoltage[sweepNo][zeroPoint:breakHere])
                data['xData'][1].append(sweepsVoltage[sweepNo][breakHere:])
                data['yData'][0].append(sweepsCurrent[sweepNo][zeroPoint:breakHere])
                data['yData'][1].append(sweepsCurrent[sweepNo][breakHere:])

                # Ion saturation Current
                Ii0 = min(data['yData'][0][i])
                data['yData'][0][i] = [y - Ii0 for y in data['yData'][0][i]] # subtract Ii0 from transition data

                # Fit the transition and saturation regions
                for k in range(2):
                    # Fit the data
                    if k == 0:
                        Params, cov = scipy.optimize.curve_fit(transitionFunc, data['xData'][0][i], data['yData'][0][i], maxfev=100000, bounds=data['fitBounds'][0], p0=data['fitGuess'][0])
                    else:
                        Params, cov = scipy.optimize.curve_fit(saturationFunc, data['xData'][1][i], data['yData'][1][i], maxfev=100000, bounds=data['fitBounds'][1], p0=data['fitGuess'][1])

                    data['fitParams'][k].append(Params)

                    # Calculate Chi^2
                    nu = (len(data['xData'][k][i]) - 3)
                    data['ChiSquare'][k].append((1/nu)*sum([(data['yData'][k][i][j] - transitionFunc(data['xData'][k][i][j], Params[0], Params[1], Params[2]))**(2) /(np.sqrt(np.abs(data['yData'][k][i][j]) + np.abs(data['xData'][k][i][j])) ) for j in range(len(data['yData'][k][i]))]))

                    # Reset transition data back
                    data['yData'][k][i] = [y + Ii0 for y in data['yData'][k][i]] if k == 0 else data['yData'][k][i]# subtract Ii0 from transition data

            ###########################
            # --- FIND THE BEST FIT ---
            ###########################

            # --- BEST TRANSITION ---
            minChi = 50
            for i, val in enumerate(data['ChiSquare'][0]):
                if np.abs(val - 1) <= np.abs(minChi - 1):
                    transData = [data['ChiSquare'][0][i], data['fitParams'][0][i], data['fitParams'][1][i],i,Vsp_range[i]]

            # --- BEST SATURATION ---
            minChi = 50
            for i, val in enumerate(data['ChiSquare'][1]):
                if np.abs(val - 1) <= np.abs(minChi - 1):
                    satData = [data['ChiSquare'][1][i], data['fitParams'][0][i], data['fitParams'][1][i],i,Vsp_range[i]]

            # Store/collect fit data
            parameters = [transData[1], satData[2]]
            chiSq = [transData[0], satData[0]]
            for i in range(2):
                params = parameters[i]
                a0 = round(params[0], rounding)
                a1 = round(params[1], rounding)
                a2 = round(params[2], rounding)
                n_0 = round(10 ** (-6) * (params[0] / np.sqrt(params[2])), rounding)

                numDensity[i].append(n_0)
                electronTemp[i].append(a2)
                plasmaPotential[i].append(a1)
                chiSquares[i].append(chiSq[i])
                parameterEpoch[i].append(pycdf.lib.tt2000_to_datetime(sweepsCurrent_epoch[sweepNo]))

            # PLOT THE BEST FITS OF BOTH CURVES
            if plotBestFit:
                datum = [transData,satData]
                titles = ['transition','saturation']

                for i, bestFitData in enumerate(datum):

                    # --- PLOT BEST FIT ---
                    fig, ax = plt.subplots(2)
                    #transition plot
                    ax[0].set_title(data['timeStamp'] + ' UTC\n' +f'{titles[i]} '+'$\chi$'+f' : {bestFitData[0]}' + '\n $V_{sp}$: '+ f'{bestFitData[4]} V')
                    ax[0].scatter(data['xData'][0][bestFitData[3]],data['yData'][0][bestFitData[3]])
                    params = bestFitData[1]
                    yData_fitted = [transitionFunc(x, params[0], params[1], params[2])+ Ii0 for x in data['xData'][0][bestFitData[3]]]
                    ax[0].plot(data['xData'][0][bestFitData[3]], yData_fitted)

                    a0 = round(params[0], rounding)
                    a1 = round(params[1], rounding)
                    a2 = round(params[2], rounding)
                    n_0 = round(10 ** (-6) * (params[0] / np.sqrt(params[2])), rounding)
                    ax[0].legend(['$a_{0} = n_{0} T_{e[eV]}^{1/2}$: ' + f'{a0} \n' + '$a_{1} = V_{sp}$: ' + f'{a1}' + ' V' + '\n' + '$a_{2} = T_{e[eV]}$: ' + f'{a2}' + ' eV' + '\n' + '$T_{i[eV]}$= ' + f'{Ti}' + ' eV' + '\n' + '$n_{0} =$' + f'{n_0} ' + ' $cm^{-3}$' + ''])
                    ax[0].set_xlabel('Probe Voltage [V]')
                    ax[0].set_ylabel('Probe Current [nA]')
                    ax[0].xaxis.set_major_locator(MultipleLocator(0.1))

                    # saturation plot
                    ax[1].scatter(data['xData'][1][bestFitData[3]], data['yData'][1][bestFitData[3]])
                    params = bestFitData[2]
                    yData_fitted = [saturationFunc(x, params[0], params[1], params[2]) for x in data['xData'][1][bestFitData[3]]]
                    ax[1].plot(data['xData'][1][bestFitData[3]], yData_fitted)
                    a0 = round(params[0], rounding)
                    a1 = round(params[1], rounding)
                    a2 = round(params[2], rounding)
                    n_0 = round(10 ** (-6) * (params[0] / np.sqrt(params[2])), rounding)
                    ax[1].legend(['$a_{0} = n_{0} T_{e[eV]}^{1/2}$: ' + f'{a0} \n' + '$a_{1} = V_{sp}$: ' + f'{a1}' + ' V' + '\n' + '$a_{2} = T_{e[eV]}$: ' + f'{a2}' + ' eV' + '\n' + '$T_{i[eV]}$= ' + f'{Ti}' + ' eV' + '\n' + '$n_{0} =$' + f'{n_0} ' + ' $cm^{-3}$' + ''])
                    ax[1].set_xlabel('Probe Voltage [V]')
                    ax[1].set_ylabel('Probe Current [nA]')
                    ax[1].xaxis.set_major_locator(MultipleLocator(0.1))
                    plt.tight_layout()
                    plt.show()



        # PLOT the fit data over the whole flight
        if plotFitParamsOverTime:

            titles=['Transition','Saturation']

            # Separate Transition and Saturation Plots
            for i in range(2):
                fig, ax = plt.subplots(4, 1, sharex=True)
                fig.set_size_inches(15, 15)

                fig.suptitle(f'{titles[i]} BestFit Parameters \n ACESII {rocketID}')
                ax[0].plot(parameterEpoch[i],chiSquares[i])
                ax[0].set_ylabel('$\chi ^{2}$ ')
                ax[1].plot(parameterEpoch[i],numDensity[i])
                ax[1].set_ylabel('$n_{0}$ [cm$^{-3}$]')
                ax[2].plot(parameterEpoch[i], electronTemp[i])
                ax[2].set_ylabel('$T_{e}$ [eV]')
                ax[3].plot(parameterEpoch[i], plasmaPotential[i])
                ax[3].set_ylabel('$V_{sp}$ [V]')
                plt.savefig(rf'D:\Data\ACESII\science\Langmuir\plots\{fliers[wRocket - 4]}\LP_Parameters_{fliers[wRocket - 4]}Flyer_{titles[i]}.png')

            # Overlay Transition and Saturation Data
            fig, ax = plt.subplots(4, 1, sharex=True)
            fig.set_size_inches(15, 15)
            for i in range(2):
                fig.suptitle(f'{titles[i]} BestFit Parameters \n ACESII {rocketID}')
                ax[0].plot(parameterEpoch[i], chiSquares[i])
                ax[0].set_ylabel('$\chi ^{2}$ ')
                ax[1].plot(parameterEpoch[i], numDensity[i])
                ax[1].set_ylabel('$n_{0}$ [cm$^{-3}$]')
                ax[2].plot(parameterEpoch[i], electronTemp[i])
                ax[2].set_ylabel('$T_{e}$ [eV]')
                ax[3].plot(parameterEpoch[i], plasmaPotential[i])
                ax[3].set_ylabel('$V_{sp}$ [V]')

            for i in range(4):
                ax[i].legend(titles)

            plt.savefig(rf'D:\Data\ACESII\science\Langmuir\plots\{fliers[wRocket - 4]}\LP_Parameters_{fliers[wRocket - 4]}Flyer_OVERLAY.png')

        #####################
        # --- OUTPUT DATA ---
        #####################
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