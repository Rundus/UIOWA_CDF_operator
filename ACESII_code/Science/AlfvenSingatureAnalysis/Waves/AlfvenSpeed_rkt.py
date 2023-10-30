# --- AlfvenSpeed_rkt.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Determine the Alfven speed over the flight of the payloads using
# the density data and magnetometer data



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt

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
wRocket = 4

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = []

modifier = ''
inputPath_modifier_density = 'science/Langmuir' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
wFile_ni = 0
inputPath_modifier_mag = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
wFile_mag = 0
outputPath_modifier = 'science/AlfvenSpeed_rkt' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
outputData = True


# --- reduce data ---
reduceData = True
targetTimes = [dt.datetime(2022,11,20,17,24,20,00),dt.datetime(2022,11,20,17,28,00,00)]

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import u0,IonMasses

def AlfvenSpeed_rkt(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)


    inputFiles_mag = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}{modifier}\*RingCore_rktFrm*')
    inputFiles_density = glob(f'{rocketFolderPath}{inputPath_modifier_density}\{fliers[wflyer]}{modifier}\*fixed*.cdf')

    fileoutName = f'ACESII_{rocketID}_AlfvenSpeed_flight.cdf'

    if justPrintFileNames:
        print('--- MAG ---')
        for i, file in enumerate(inputFiles_mag):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_mag[i], round(getsize(file) / (10 ** 6), 1)))
        print('\n--- Density ---')
        for i, file in enumerate(inputFiles_density):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_density[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Calculating Alfven Speed' + color.END)

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading data from mag File')
        data_dict_mag = loadDictFromFile(inputFiles_mag[wFile_mag],{},reduceData=reduceData,targetTimes=targetTimes)
        Done(start_time)

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading data from density File')
        data_dict_density = loadDictFromFile(inputFiles_density[wFile_ni], {},reduceData=reduceData,targetTimes=targetTimes)
        Done(start_time)

        # reduce data
        lowCut, highCut = np.abs(data_dict_mag['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_mag['Epoch'][0] - targetTimes[1]).argmin()

        for key, val in data_dict_mag.items():
            data_dict_mag[key][0] = data_dict_mag[key][0][lowCut:highCut]

        lowCut, highCut = np.abs(data_dict_density['Epoch'][0] - targetTimes[0]).argmin(), np.abs( data_dict_density['Epoch'][0] - targetTimes[1]).argmin()
        for key, val in data_dict_density.items():
            data_dict_density[key][0] = data_dict_density[key][0][lowCut:highCut]


        ###################################
        # --- DOWNSAMPLE DENSITY TO MAG ---
        ###################################
        prgMsg('Downsampling Density')
        from ACESII_code.class_var_func import InterpolateDataDict,dateTimetoTT2000

        data_dict_densityInterp = InterpolateDataDict(InputDataDict=data_dict_density,
                                                      InputEpochArray=dateTimetoTT2000(data_dict_density['Epoch'][0], False),
                                                      wKeys=['ni', 'Epoch'],
                                                      targetEpochArray=dateTimetoTT2000(data_dict_mag['Epoch'][0],False))
        Done(start_time)



        ####################################
        # --- CALCULATE THE ALFVEN SPEED ---
        ####################################

        AlfvenSpeed = np.array([
            ((1E-9)*data_dict_mag['Bmag'][0][i])/np.sqrt(u0*(100**3)*data_dict_densityInterp['ni'][0][i]*IonMasses[0]) for i in range(len(data_dict_mag['Epoch'][0]))
        ])


        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('creating output file')

            data_dict_output = {}

            data = AlfvenSpeed
            varAttrs = {'LABLAXIS': "AlfvenSpeed", 'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None,
                        'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'm/s',
                        'VALIDMIN': data.min(), 'VALIDMAX': data.max(),
                        'VAR_TYPE': 'data', 'SCALETYP': 'linear'}

            data_dict_output = {**data_dict_output, **{'AlfvenSpeed': [data, varAttrs]}}

            Epoch_output = deepcopy(data_dict_mag['Epoch'])

            Epoch_output[1]['VAR_TYPE'] = 'support_data'

            data_dict_output = {**data_dict_output, **{'Epoch': Epoch_output}}

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict_output, outputModelData, globalAttrsMod, 'Attitude')

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

if len(glob(f'{rocketFolderPath}{inputPath_modifier_density}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        AlfvenSpeed_rkt(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        AlfvenSpeed_rkt(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)