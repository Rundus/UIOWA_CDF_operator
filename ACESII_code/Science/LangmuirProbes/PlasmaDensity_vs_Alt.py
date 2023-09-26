# --- PlasmaDensity_vs_Alt.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Using the derived plasma density from the Fixed Langmuir Probe, we can
# see how this density looks when plotted vs altitude



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal

from ACESII_code.myImports import *

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
wFiles = []

modifier = ''
inputPath_modifier = 'science/Langmuir' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_attitude = 'attitude'
outputPath_modifier = '' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.interpolate import CubicSpline


def PlasmaDensity_vs_Alt(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*Langmuir*')
    inputFiles_attitude = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[wflyer]}{modifier}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]

    fileoutName = f'ACESII_{rocketID}_PlasmaDensity_vs_Alt.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Plotting Density vs Altutde information for ACESII_{rocketID}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict_LP = loadDictFromFile(inputFiles[wFile], {})
        Epoch_LP_fixed = np.array([pycdf.lib.datetime_to_tt2000(data_dict_LP['fixed_Epoch'][0][i]) for i in (range(len(data_dict_LP['fixed_Epoch'][0])))])
        ni = np.array(data_dict_LP['fixed_ni_density'][0])
        Done(start_time)

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading data from {inputPath_modifier_attitude} Files')
        data_dict_attitude = loadDictFromFile(inputFiles_attitude[wFile], {})
        data_dict_attitude['Epoch'][0] = np.array(data_dict_attitude['Epoch'][0])
        Done(start_time)

        #########################################################
        # --- Interpolate the Attitude data up to the LP data ---
        #########################################################

        # Attitude
        Epoch_attitude = np.array([pycdf.lib.datetime_to_tt2000(data_dict_attitude['Epoch'][0][i]) for i in range(len(data_dict_attitude['Epoch'][0]))])
        Alt = np.array(data_dict_attitude['Alt'][0])  # in meters!
        Lat = data_dict_attitude['Latgd'][0]
        Long = data_dict_attitude['Long'][0]

        normalData = [Epoch_attitude, Alt / 1000, Lat, Long]

        newData = {
            'newEpoch': Epoch_LP_fixed,
            'newAlt': [],
            'newLat': [],
            'newLong': [],
        }

        # --- DO THE INTERPOLATION ---
        counter = 0
        for key, newDataList in tqdm(newData.items()):
            if key != 'newEpoch':
                # --- cubic interpolation ---
                splCub = CubicSpline(Epoch_attitude, normalData[counter])

                # evaluate the interpolation at all the epoch_mag points
                newData[key] = np.array([splCub(timeVal) for timeVal in Epoch_LP_fixed])

            counter += 1

        Done(start_time)

        # do some quality assurance steps
        badIndicies = []
        for i in range(len(ni)):
            if np.abs(ni[i]) > 1E8:
                badIndicies.append(i)

            if np.abs(newData['newAlt'][i]) > 1E8:
                badIndicies.append(i)

        newAlt = np.delete(newData['newAlt'],badIndicies)
        newNi = np.delete(ni,badIndicies)


        # --- Smooth the density data to remove spin effects ---
        newNi = scipy.signal.savgol_filter(newNi, window_length=10000, polyorder=2, mode='nearest')


        # --- Split the data into upleg and downleg ---

        # find the max in the altitude
        splitIndex = np.abs(newAlt -  newAlt.max()).argmin()

        newAlt_upleg = newAlt[:splitIndex]
        newNi_upleg = newNi[:splitIndex]

        newAlt_downleg = newAlt[splitIndex:]
        newNi_downleg = newNi[splitIndex:]

        #################################
        # --- Plot the Density vs Alt ---
        #################################


        xData = [newNi_upleg,newNi_downleg]
        yData = [newAlt_upleg,newAlt_downleg]
        labels= ['Upleg','Downleg']
        for i in range(2):
            fig, ax = plt.subplots()
            fig.suptitle(f'ACESII {rocketID} {labels[i]}\n'
                         r'$T_{i}$ = 0.1eV')
            ax.scatter(xData[i],yData[i])
            ax.set_xlabel('Density [cm$^{-3}$]')
            ax.set_ylabel('Altitude [km]')
            ax.set_ylim(70, 410) if wRocket ==4 else ax.set_ylim(70, 200)
            ax.set_xlim(1E1, 1E7)
            ax.set_xscale('log')
            ax.grid(True)
        plt.show()






        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('Creating output file')

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict, L2_TRICE_Quick, globalAttrsMod, "Langmuir")

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
        PlasmaDensity_vs_Alt(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            PlasmaDensity_vs_Alt(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            PlasmaDensity_vs_Alt(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)