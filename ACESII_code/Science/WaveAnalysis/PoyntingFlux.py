# --- PoyntingFlux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Determine the PoyntingFLux of the data using E-Field and B-Field Measurements



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from ACESII_code.myImports import *

start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = True
printMagFiles = True
printElecFiles = True

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles_elec = [0]
wFiles_mag = [0]

useSpunData = True

modifier = ''
inputPath_modifier = 'l1' if useSpunData else 'l2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science/PoyntingFlux' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = False



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def PoyntingFlux(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles_elec = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*E_Field*')
    inputFiles_mag = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*RingCore_rktFrm*')

    fileoutName = f'ACESII_{rocketID}_PoyntingFlux_spun.cdf' if useSpunData else f'ACESII_{rocketID}_PoyntingFlux_ENU.cdf'

    if justPrintFileNames:
        if printMagFiles:
            for i, file in enumerate(inputFiles_mag):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_mag[i], round(getsize(file) / (10 ** 6), 1)))
        elif printElecFiles:
            for i, file in enumerate(inputFiles_elec):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_elec[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Calculating Poynting flux for ACESII {rocketID}' + color.END)

        # --- get the data from the mag file ---
        prgMsg(f'Loading data from mag Files')
        data_dict_mag = loadDictFromFile(inputFiles_mag[wFiles_mag[0]],{})
        data_dict_mag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag['Epoch'][0]])
        Done(start_time)

        # --- get the data from the mag file ---
        prgMsg(f'Loading data from Electric Field Files')
        data_dict_elec = loadDictFromFile(inputFiles_elec[wFiles_elec[0]], {})
        data_dict_elec['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_elec['Epoch'][0]])
        Done(start_time)

        ##################################
        # --- INTERPOLATE MAG ONTO EFI ---
        ##################################

        #############################
        # --- FILTER B-Field Data ---
        #############################

        #################################
        # --- CALCULATE POYNTING FLUX ---
        #################################



        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('Creating output file')

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod,wInstr[1])

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
        PoyntingFlux(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            PoyntingFlux(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            PoyntingFlux(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)