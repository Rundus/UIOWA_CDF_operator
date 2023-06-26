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
wFiles = [1,2]

modifier = ''
inputPath_modifier = 'mag' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science/Magnetometer_filtered' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# if using ENU coordniates
useENUcoordinates = True

##################
# FILTER TOGGLES #
##################
filtOrder = 5 # order of filter
cutoff_Freq = 2 # cut off frequency where gain has reached -3dB
dataSampleFreq = 128 # sample rate of the data

# Plot the filtered data
plotFilteredData = False
plotFreqSpectrogram = False

# output the data
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg, L2_TRICE_Quick, outputCDFdata
from glob import glob
from os.path import getsize
from scipy.signal import butter, filtfilt,spectrogram
import matplotlib.pyplot as plt

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

    # determine which instrument the file corresponds to:
    descriptiorNam = ['eepaa', 'leesa', 'iepaa', 'lp','Tesseract_Geo', 'RingCore_Geo','Tesseract_spun_instrFrame','RingCore_spun_instrFrame','Tesseract', 'RingCore']
    for index, instr in enumerate(descriptiorNam):
        if instr in dataFile_name:
            wInstr = [index, instr, descriptiorNam[index]]
            break

    fileoutName = f'ACESII_{rocketID}_{wInstr[1]}_HighPass.cdf'

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

        data_dict['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][i]) for i in (range(len(data_dict['Epoch'][0])))])

        Done(start_time)

        ##################################
        # --- FITLER AND PLOT MAG DATA ---
        ##################################

        def butter_highpass(cutoff, fs, order):
            nyq = 0.5 * fs
            normal_cutoff = cutoff / nyq
            b, a = butter(order, normal_cutoff, btype='high', analog=False)
            return b, a

        def butter_highpass_filter(data, cutoff, fs, order):
            b, a = butter_highpass(cutoff, fs, order=order)
            y = filtfilt(b, a, data)
            return y


        # format: [Bx,By,Bz]
        if useENUcoordinates:
            magAxes = ['B_east', 'B_north', 'B_up']
        else:
            magAxes = ['Bx','By','Bz']
        conditionedData = np.array([butter_highpass_filter(data_dict[axes][0], cutoff_Freq, dataSampleFreq, filtOrder) for axes in magAxes])
        time = data_dict['Epoch'][0]

        # plot the data
        if plotFilteredData:
            for i in range(len(conditionedData)):

                fig,ax = plt.subplots(2)
                plt.title(f'{magAxes[i]}_filtered vs {magAxes[i]}\nCutoff Freq: {cutoff_Freq} Hz \nOrder: {filtOrder}')
                ax[0].plot(time, conditionedData[i])
                ax[0].set_ylabel(f'{magAxes[i]}_filtered [nT]')
                ax[1].plot(time,data_dict[magAxes[i]][0])
                ax[1].set_ylabel(f'{magAxes[i]} [nT]')
                ax[1].set_xlabel('Epoch [ns]')
                plt.show()

        if plotFreqSpectrogram:
            for i in range(len(conditionedData)):

                fig, ax = plt.subplots()

                f, t, Sxx = spectrogram(conditionedData[i], dataSampleFreq)
                cmap = ax.pcolormesh(t, f, Sxx,cmap='turbo',vmin=1E-3,vmax=1E1,norm='log')

                ax.set_ylabel('Frequency [Hz]')
                ax.set_xlabel('Time')
                ax.set_title(f'{magAxes[i]}_filtered Spectrogram\nCutoff Freq: {cutoff_Freq} Hz \nOrder: {filtOrder}')
                cbar = fig.colorbar(cmap)
                cbar.set_label('Spectral Density')
                plt.show()


        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            # Add filtered data to data dict

            for i,data in enumerate(conditionedData):
                data_dict = {**data_dict, **{f'{magAxes[i]}_filtered':
                                                 [data, {'LABLAXIS': f'{magAxes[i]}_filtered',
                                                                        'DEPEND_0': 'Epoch',
                                                                        'DEPEND_1': None,
                                                                        'DEPEND_2': None,
                                                                        'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                                        'UNITS': 'nT',
                                                                        'VALIDMIN': data.min(),
                                                                        'VALIDMAX': data.max(),
                                                                        'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
            # --- --- --- --- --- --- ---
            # --- WRITE OUT THE DATA ---
            # --- --- --- --- --- --- ---
            prgMsg('Creating output file')

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            globalAttrsMod['Descriptor'] = wInstr[1]
            outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, wInstr[1])

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