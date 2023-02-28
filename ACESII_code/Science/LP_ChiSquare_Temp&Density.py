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
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'L2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\Langmuir' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# --- Fitting Toggles ---
unitConv = 10**(9) # converts from A to nA

wSweeps = [70] # [] --> all sweeps, [#1,#2,...] specific sweeps

plotSaturationRegion = True
plotTransitionRegion = False

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
N_Ti, N_Vsp, N_Te, N_n0 = 20, 5, 100, 100
Ti_range = np.linspace(0.001, 1, N_Ti)
Vsp_range = np.linspace(0, 2.0, N_Vsp)
Te_range = np.linspace(1, 10000, N_Te)
n0_range = np.linspace(0.01, 1000, N_n0)
gridSearchRanges = [n0_range, Te_range]

##############################
# --- FITTED CHI FUNCTIONS ---
##############################
rocketAttrs, b, c = ACES_mission_dicts()

e_coefficient = (((q0 ** 3) * (rocketAttrs.LP_probe_areas[0][0] ** (2))) / (8 * np.pi * m_e)) ** (1 / 2)
i_coefficient = ((1 * (q0 ** 3) * (rocketAttrs.LP_probe_areas[0][0] ** (2))) / (8 * np.pi * IonMasses[0])) ** (1 / 2)

# def transitionFunc(x, a0, a1, a2, a3):
#     y = (e_coefficient) * (a1 ** (1 / 2)) * (a2) * np.exp((x - a0) / a1) - a2 * (
#                 (a3 * (q0 ** 3) * (rocketAttrs.LP_probe_areas[0][0] ** (2))) / (8 * np.pi * IonMasses[0])) ** (
#                     1 / 2)
#     return y

def transitionFunc(x, a0, a1, a2):
    y = unitConv*(e_coefficient) * a0 * np.exp((x - a1) / a2)
    return y
def saturationFunc(x, a0, a1, a2, a3):
    y = unitConv*((e_coefficient) * a0  - (a0/(a2**(1/2))) * ((a3 * (q0 ** 3) * (rocketAttrs.LP_probe_areas[0][0] ** (2))) / (8 * np.pi * IonMasses[0])) ** (1 / 2) * np.exp(-((x - a1) / a3)))
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

        print(len(sweepIndices),qualityCounter)


        sweepIndices = np.array(sweepIndices)
        breakIndices = np.array(breakIndices)

        # --- store the sweep data into a single data variables filled with individual sweeps---
        prgMsg('Reorganizing Data')
        sweepsCurrent = []
        sweepsVoltage = []
        sweepsCurrent_errors = []

        for sweep in range(len(sweepIndices)):
            start = sweepIndices[sweep][0]
            end = sweepIndices[sweep][1]
            breakPoint = breakIndices[sweep]

            # Get the data and sort it
            xData1 = np.array(data_dict['swept_Voltage'][0][start:breakPoint])
            yData1 = np.array(data_dict['swept_Current'][0][start:breakPoint]) / unitConv
            yData1_errors = np.array(data_dict['swept_Current_errors'][0][start:breakPoint] / unitConv)
            yData1 = np.array([x*unitConv for _, x in sorted(zip(xData1, yData1))])
            yData1_errors = np.array([x for _, x in sorted(zip(xData1, yData1_errors))])
            xData1 = np.array(sorted(xData1))

            xData2 = data_dict['swept_Voltage'][0][breakPoint:end]
            yData2 = np.array(data_dict['swept_Current'][0][breakPoint:end]) / unitConv
            yData2_errors = np.array(data_dict['swept_Current_errors'][0][breakPoint:end] / unitConv)

            yData2 = np.array([x for _, x in sorted(zip(xData2, yData2))])
            yData2_errors = np.array([x*unitConv for _, x in sorted(zip(xData2, yData2_errors))])
            xData2 = np.array(sorted(xData2))

            # append data to lists
            sweepsVoltage.append(xData1)
            sweepsVoltage.append(xData2)
            sweepsCurrent.append(yData1)
            sweepsCurrent.append(yData2)
            sweepsCurrent_errors.append(yData1_errors)
            sweepsCurrent_errors.append(yData2_errors)

        sweepsVoltage = np.array(sweepsVoltage,dtype='object')
        sweepsCurrent = np.array(sweepsCurrent,dtype='object')
        sweepsCurrent_errors = np.array(sweepsCurrent_errors,dtype='object')
        Done(start_time)

        ###########################
        # --- APPLY GRID SEARCH ---
        ###########################

        prgMsg('Performing Grid Search')
        print('\n')

        # analyze the specific sweeps in wSweeps or do all of them
        theseSweeps = [i for i in range(len(sweepsCurrent))] if wSweeps == [] else wSweeps

        for sweepNo in theseSweeps:

            # Get all the x and y data of the two probe sweeps
            xData = sweepsVoltage[sweepNo]
            breakHere = np.abs(xData - 1.1).argmin()
            zeroPoint = np.abs(xData - 0).argmin()
            xData = sweepsVoltage[sweepNo][zeroPoint:breakHere]
            yData = sweepsCurrent[sweepNo][zeroPoint:breakHere]
            yData = [y - min(yData) for y in yData]
            yData_errors = sweepsCurrent_errors[sweepNo]
            chiParams = []


            # transition Region
            fig, ax = plt.subplots()
            ax.scatter(xData, yData)

            a = 1.5E7
            b = 0.7
            c=0.2

            bounds =([np.NINF,-3,np.NINF],[np.Inf,3,np.Inf])
            p0 = [1.77758776,-0.52435768,  0.08296951]
            params,cov = scipy.optimize.curve_fit(transitionFunc,xData,yData,maxfev=100000,bounds = bounds,p0=p0)
            yData_fitted = [transitionFunc(x, params[0], params[1], params[2]) for x in xData]
            plt.plot(xData,yData_fitted)
            print(params)
            plt.show()
            # [1.78120679 - 0.52416683  0.08296832]

            # Saturation Region
            xData = sweepsVoltage[sweepNo][breakHere:]
            yData = sweepsCurrent[sweepNo][breakHere:]
            # yData = [y - min(yData) for y in yData]
            fig, ax = plt.subplots()
            ax.scatter(xData,yData)
            print(xData[::int(len(xData)/5)])
            print(yData[::int(len(yData)/5)])
            bounds = ([np.NINF, -3, np.NINF,0.001], [np.Inf, 3, np.Inf,3])
            p0 = [2.11E8,2,0.2,0.3]
            params,cov = scipy.optimize.curve_fit(saturationFunc,xData,yData,maxfev=100000,bounds=bounds,p0=p0)
            yData_fitted = [saturationFunc(x, params[0], params[1], params[2],params[3]) for x in xData]
            print(params)
            ax.plot(xData,yData_fitted)
            plt.show()






            # Determine the transition/saturation region
            for Vsp in tqdm(Vsp_range):
                # Break the data into transition and saturation via current = 0:
                curveBreakpoint = np.abs(xData - Vsp).argmin()

                data = {'xSat': xData[curveBreakpoint:],
                        'ySat': yData[curveBreakpoint:],
                        'ySat_errors': yData_errors[curveBreakpoint:],
                        'xTrans': xData[0:curveBreakpoint],
                        'yTrans': yData[0:curveBreakpoint],
                        'yTrans_errors': yData_errors[0:curveBreakpoint]}

                # Determine and subtract off the ion current contribution
                Ii0 = min(data['yTrans'])
                data['yTrans'] = np.array([y - Ii0 for y in data['yTrans']])

                # just do the transition region
                nu = len(data['xTrans']) - 3
                for Te, n0 in itertools.product(*gridSearchRanges):
                    chisquare = (1/nu) * sum(np.array([((data['yTrans'][i] - transitionFunc(data['yTrans'][i], Vsp, Te, n0))**(2))/((data['yTrans_errors'][i])**(2)) for i in range(len(data['xTrans']))]))
                    chiParams.append([chisquare, Vsp, Te, n0])

        Done(start_time)

        ###########################
        # --- PLOT THE BEST FIT ---
        ###########################
        minChi = [500,0]
        for i,chi in enumerate(chiParams):

            oneDistance = np.abs(1 - chi[0])

            if oneDistance < minChi[0]:
                minChi[0] = oneDistance
                minChi[1] = i

        print(minChi)
        print(chiParams[minChi[1]])



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