# --- L1_to_L1magCalESA.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Get the "true" pitch angle that was calculated by the RingCore magnetometer and
# then use it to re-bin the EEPAA data and collect ChiSquare-Calibrations data
# reduce dataset to only include data after ACS despun rocket since this is where the magPitchAngles are nominal.
# Value is chosen based off of the attitude data variable "SpinAlignmentToB"

# The ChiSquare Calibration values are chosen only AFTER the
# ACS has despun the rocket for nominal data collection, otherwise the
# magnetic field will be very out of line with the ESA detector



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

inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_magPitch = 'science\ESA_magPitch_calibration' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\ESA_magPitch_calibration' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
modifier = ''

# --- Primary Function TOGGLES ---
useTesseractData = False
AnalyzeTruePitchAngle = False # set to False to perform magCalibration

# --- PITCH BIN RESORTING TOGGLES ---
degWidth = 5 # controls the size of the acceptance range for re-bining

# --- CHI SQUARE CALIBRATION TOGGLES ---
outputPath_modifier_chiCal = 'science/ESA_ChiSquare_calibration'

# --- Output Data TOGGLES ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.optimize import curve_fit
from ACESII_code.class_var_func import loadCDFdata,outputCDFdata,L1_TRICE_Quick
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg
from glob import glob
from os.path import getsize
from datetime import datetime

setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


def L1_to_L1magCalESA(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L1'
    outputModelData = L1_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')
    inputFiles_magPitch = glob(f'{rocketFolderPath}{inputPath_modifier_magPitch}\{fliers[wflyer]}{modifier}\*magPitch*')

    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    input_names_magPitch = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]

    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            wInstr = [index, instr]

    fileoutName = f'ACESII_{rocketID}_l1_{wInstr[1]}_magCal.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made L2: {:3s} '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the l1 ESA file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')

        data_dict_esa = loadCDFdata(inputFiles, wFile)

        data_dict_esa['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_esa['Epoch_esa'][0][i]) for i in range(len(data_dict_esa['Epoch_esa'][0]))])

        Done(start_time)

        # --- get the data from the MagPitch file ---
        prgMsg(f'Loading data from MagPitch Files')

        this_magPitchFile = ''
        for file in inputFiles_magPitch:
            if wInstr[1] in file:
                this_magPitchFile = [file]

        data_dict_magPitch = loadCDFdata(this_magPitchFile, 0)

        data_dict_magPitch['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_magPitch['Epoch_esa'][0][i]) for i in (range(len(data_dict_magPitch['Epoch_esa'][0])))])

        Done(start_time)

        if AnalyzeTruePitchAngle:

            prgMsg('Studying true pitch angle')

            Mag_Calculated_Pitch_Angle = data_dict_magPitch['Mag_Calculated_Pitch_Angle'][0]

            # reduce dataset to only data after point ACS down-spun the data

            targetTime = pycdf.lib.datetime_to_tt2000(datetime(2022,11,20,17,21,50,000))
            Mag_Calculated_Pitch_Angle = Mag_Calculated_Pitch_Angle[np.abs(np.array(data_dict_magPitch['Epoch_esa'][0]) - targetTime).argmin():]

            # --- get statistics on each pin bin ---
            binMeans = []
            binMedians = []
            binStdDevs = []

            for ptch in range(len(Mag_Calculated_Pitch_Angle[0])): # loop over all pitch angles

                # get all the data for one pitch bin
                dataSet = [Mag_Calculated_Pitch_Angle[tme][ptch][engy] for engy in range(len(Mag_Calculated_Pitch_Angle[0][0])) for tme in range(len(Mag_Calculated_Pitch_Angle))]
                binMeans.append(round(np.mean(dataSet), 1))
                binMedians.append(round(np.median(dataSet), 1))
                binStdDevs.append(round(np.std(dataSet), 1))

            print('\n\n')
            for i in range(len(binMeans)):
                print(i,binMeans[i],binMedians[i],binStdDevs[i])

            means = []
            stdDevs = []
            medians = []
            for i in range(len(Mag_Calculated_Pitch_Angle[0])):
                means.append(f'{i} {binMeans[i]}')
                stdDevs.append(f'{i} {binStdDevs[i]}')
                medians.append(f'{i} {binMeans[i]}')

            fig, ax = plt.subplots(3)
            ax[0].bar(means, binMeans)
            ax[0].set_title('means')
            ax[1].bar(stdDevs, binStdDevs)
            ax[1].set_title('stdDevs')
            ax[2].bar(medians, binMedians)
            ax[2].set_title('medians')
            for i in range(3):
                ax[i].grid()

            plt.tight_layout()
            plt.show()

            Done(start_time)

        else:
            #############################
            # --- Re-Bin the ESA data ---
            #############################
            prgMsg('Rebinnig the ESA data')

            # reduce dataset to only include data after ACS despun rocket since this is where the magPitchAngles are nominal. Value is chosen based off of the attitude data variable "SpinAlignmentToB"
            targetIndexTimes = [
                np.abs(data_dict_esa['Epoch_esa'][0] - pycdf.lib.datetime_to_tt2000(datetime(2022, 11, 20, 17, 21, 48, 300))  ).argmin()  ,
                np.abs(data_dict_esa['Epoch_esa'][0] - pycdf.lib.datetime_to_tt2000(datetime(2022, 11, 20, 17, 23, 46, 750))  ).argmin()
            ]

            esaRanges = [range(len(data_dict_esa[wInstr[1]][0])), range(len(data_dict_esa[wInstr[1]][0][0])), range(len(data_dict_esa[wInstr[1]][0][0][0]))]
            esaData = data_dict_esa[wInstr[1]][0]
            magPitch = data_dict_magPitch['Mag_Calculated_Pitch_Angle'][0]
            padAngle = data_dict_esa['Pitch_Angle'][0]

            # --- prepare data for re-binning ---
            esaDataSorted = [[[[] for engy in esaRanges[2]] for ptch in esaRanges[1]] for tme in esaRanges[0]]

            # --- --- --- --- --- --- --- --- --
            # --- GET CHI^2 CALIBRATION DATA ---
            # --- --- --- --- --- --- --- --- --

            #appended data will look like: [uncalPoint, prinPoint, uncalPad_index, prinPad_index]
            calData = []

            # for two datapoints at the same time/energy but different pitch pads whose "True" pitch angle are <5deg differenet, store this data for calibrations

            # --- rebin the data and collect calibration ---
            for tme, ptch, engy in tqdm(itertools.product(*esaRanges)):

                pitchVal = magPitch[tme][ptch][engy]
                esaVal = esaData[tme][ptch][engy]

                if (0 < ptch < 20): # 0deg - 180 deg case
                    pad = padAngle[ptch]
                elif ptch == 20: # 190 deg case
                    pad = 170
                elif ptch == 0: # -10 deg case
                    pad = 10

                if (pad - degWidth) <= pitchVal <= (pad + degWidth): # if you're within the right pad, place it where it was
                    esaDataSorted[tme][ptch][engy].append(esaVal)
                elif np.abs(pitchVal - pad) > degWidth: # if pitchVal is further than "degWidth" away from its assigned pad
                    index = np.abs(padAngle - pitchVal).argmin() # find which bin to put the value in
                    esaDataSorted[tme][index][engy].append(esaVal)  # put it in the right bin

                    # CALIBRATION: check if pitchval_uncal needs to be moved but pitchval_prin doesn't
                    if np.abs(magPitch[tme][index][engy] - padAngle[index]) <= degWidth:

                        # check if values occurs during the time after ACS despun rocket for maiden voyage.
                        if tme >= targetIndexTimes[wRocket - 4]:
                            # [uncalPoint, prinPoint, uncalPad_index, prinPad_index]
                            calData.append([esaData[tme][ptch][engy], esaData[tme][index][engy], ptch, index])

            Done(start_time)

            # --- --- --- --- --- --- --- --- --- --- -
            # --- REDUCE DATALISTS TO SINGLE VALUES ---
            # --- --- --- --- --- --- --- --- --- --- -
            prgMsg('Reducing DataSet')
            for tme, ptch, engy in itertools.product(*esaRanges): # take the average of all the lists in the data
                if len(esaDataSorted[tme][ptch][engy]) != 0:
                    esaDataSorted[tme][ptch][engy] = int(sum(esaDataSorted[tme][ptch][engy])/len(esaDataSorted[tme][ptch][engy]))
                else:
                    esaDataSorted[tme][ptch][engy] = -1e30

            esaDataSorted = np.array(esaDataSorted)

            Done(start_time)

            # --- --- --- --- --- --- --- ---
            # --- PREPARE DATA FOR OUTPUT ---
            # --- --- --- --- --- --- --- ---

            if outputData:

                data_dict = {}
                data_dict = {**data_dict, **{f'{wInstr[1]}':
                                                 [esaDataSorted, {'LABLAXIS': f'{wInstr[1]}',
                                                              'DEPEND_0': 'Epoch_esa', 'DEPEND_1': 'Pitch_Angle',
                                                              'DEPEND_2': 'Energy',
                                                              'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                              'UNITS': 'counts',
                                                              'VALIDMIN': esaDataSorted.min(), 'VALIDMAX': esaDataSorted.max(),
                                                              'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
                # get Energy and Pitch angle into data_dict
                for key, val in data_dict_esa.items():
                    if key not in [f'{wInstr[1]}']:
                        data_dict = {**data_dict,**{key:data_dict_esa[key]}}

                # --- --- --- --- --- --- --- --- ---
                # --- WRITE OUT THE L1 MagCal DATA ---
                # --- --- --- --- --- --- --- --- ---
                prgMsg('Creating output file')
                outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'
                globalAttrsMod['Descriptor'] = rocketAttrs.InstrNames_Full[wInstr[0]]
                outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, wInstr[1])

                Done(start_time)

                # --- --- --- --- --- --- --- --- --- -
                # --- WRITE OUT THE CHISQR CAL DATA ---
                # --- --- --- --- --- --- --- --- --- -
                # store the raw ChiSquare calibration data
                del data_dict
                calData = np.array(calData)
                data_dict = {f'ChiSquare_calData': [calData, {'LABLAXIS': 'ChiSquareCalData',
                                                            'DEPEND_0': None, 'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                            'UNITS': 'counts',
                                                            'VALIDMIN': calData.min(),
                                                            'VALIDMAX': calData.max(),
                                                            'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}

                fileoutName = f'ACESII_{rocketID}_{wInstr[1]}_ChiSquareData.cdf'
                outputPath = f'{rocketFolderPath}{outputPath_modifier_chiCal}\{fliers[wflyer]}\\{fileoutName}'
                globalAttrsMod['Descriptor'] = rocketAttrs.InstrNames_Full[wInstr[0]]
                outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, wInstr[1])







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
        L1_to_L1magCalESA(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            L1_to_L1magCalESA(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            L1_to_L1magCalESA(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)