# --- RingCore_Despin.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Use the calibrated Spun RingCore data and use IGRF model to
# remove the background field and coning in the magnetometer data.

# To do this we need to interpolate the attitude solution to match the time series
# of the magnetometer data, then spin the IGRF into the rocket's frame, THEN account
# for the position of the magnetometer itself to adjust for coning



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

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


modifier = ''
inputPath_modifier = 'mag'
inputPath_modifier_attitude = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'mag' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import pyIGRF

def RingCore_Despin(wRocket, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    ModelData = L0_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*_Spun*')
    inputFiles_attitude = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[wflyer]}{modifier}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    input_names_attitude = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles_attitude]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    input_names_searchable_attitude = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier_attitude.lower() + '_', '').replace('_v00', '') for ifile in input_names_attitude]

    fileoutName = f'ACESII_{rocketID}_RingCore_DeSpun_c.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {fileoutName}' + color.END)
        print('[0]   ' + str(round(getsize(inputFiles[0]) / (10 ** 6), 1)) + 'MiB')

        ###########################
        # --- Get the Variables ---
        ###########################

        # --- get the Magnetometer Data ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict_mag = loadDictFromFile(inputFiles[0], {})
        Done(start_time)

        # --- get the Attitude Data ---
        prgMsg(f'Loading data from {inputPath_modifier_attitude} Files')
        data_dict_attitude = loadDictFromFile(inputFiles_attitude[0], {})
        Done(start_time)


        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ########################################
        #### DEP-SPIN THE MAGNETOMETER DATA ####
        ########################################
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        # --- [1] Match Attitude and Magnetometer Epoch ---
        prgMsg('Collecting Variables')
        # assign all the variables
        alt = np.array(data_dict_attitude['Alt'][0])/1000
        long = data_dict_attitude['Long'][0]
        lat = data_dict_attitude['Latgd'][0]
        Epoch_attitude = data_dict_attitude['Epoch'][0]
        Epoch_mag = data_dict_mag['Epoch'][0]
        B = np.array([np.array([data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]]) for i in range(len(Epoch_mag))])

        # create new altitude that accounts for the position of the magnetometer
        r = 0/1000 # corresponds to the length (in km) that the magnetometer was away from the central gyroscope/thing that measured altitude

        X_Az = data_dict_attitude['X_Az'][0]
        X_El = data_dict_attitude['X_El'][0]

        from ACESII_code.class_var_func import lat_to_meter,meter_to_long

        adjustedAlt = np.array([
            (r*np.cos(np.radians(90 - X_El[i]))) + alt[i] for i in range(len(Epoch_attitude))
        ])
        adjustedLat = np.array([
            lat[i] + (1/lat_to_meter)*r*np.sin(np.radians(90 - X_El[i]))*np.cos(np.radians(90 - X_Az[i])) for i in range(len(Epoch_attitude))
        ])

        adjustedLong = np.array([
            long[i] + meter_to_long(r * np.sin(np.radians(90 - X_El[i])) * np.sin(np.radians(90 - X_Az[i]))) for i in range(len(Epoch_attitude))
        ])


        # -- Output order forpyIGRF.igrf_value ---
        # [0] Declination (+ E | - W)
        # [1] Inclination (+ D | - U)
        # [2] Horizontal Intensity
        # [3] North Comp (+ N | - S)
        # [4] East Comp (+ E | - W)
        # [5] Vertical Comp (+ D | - U)
        # [6] Total Field

        date = 2022 + 323 / 365  # Corresponds to 11/20/2022
        IGRF = np.array([pyIGRF.igrf_value(adjustedLat[i], adjustedLong[i], adjustedAlt[i], date) for i in range(len(Epoch_attitude))])
        IGRF_B = np.array([[IGRF[i][4], IGRF[i][3], -1*IGRF[i][5]] for i in range(len(Epoch_attitude))])

        DCMMatrix = np.array([
            [
                [data_dict_attitude['a11'][0][i], data_dict_attitude['a12'][0][i], data_dict_attitude['a13'][0][i]],
                [data_dict_attitude['a21'][0][i], data_dict_attitude['a22'][0][i], data_dict_attitude['a23'][0][i]],
                [data_dict_attitude['a31'][0][i], data_dict_attitude['a32'][0][i], data_dict_attitude['a33'][0][i]]
            ]
            for i in range(len(Epoch_attitude))
        ])

        DCMMatrix_inverse = np.array([np.linalg.inv(DCMMatrix[i]) for i in range(len(DCMMatrix))])

        # --- Apply inverse to get IGRF in rocket frme ---
        IGRF_B_rocket = np.array([np.matmul(DCMMatrix_inverse[i], IGRF_B[i]) for i in range(len(Epoch_attitude))])
        Done(start_time)

        # adjust the time series of the IGRF field
        # Down-sample the mag epoch to the ACS, do a linear chi-fit between the ACS and mag
        # then apply the fit to the ACS data.

        #make sure both datasets are tt2000
        Epoch_mag = np.array([pycdf.lib.datetime_to_tt2000(val) for val in Epoch_mag])
        Epoch_attitude = np.array([pycdf.lib.datetime_to_tt2000(val) for val in Epoch_attitude])

        prgMsg('Downsampling Mag')
        magIndicies = [np.abs(Epoch_mag - Epoch_attitude[i]).argmin() for i in range(len(Epoch_attitude))]
        Epoch_mag_ds = np.array([ pycdf.lib.tt2000_to_datetime(Epoch_mag[index]) for index in magIndicies])
        Epoch_attitude = np.array([pycdf.lib.tt2000_to_datetime(tme) for tme in Epoch_attitude])
        B_ds = np.array([B[index] for index in magIndicies])
        Done(start_time)



        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('Creating output file')
            data_dict = {}
            data_dict = {**data_dict, **{'Epoch_mag': [Epoch_mag_ds, {'DEPEND_0': 'Epoch_mag', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                           'FORMAT': 'I5', 'UNITS': 'ns', 'VALIDMIN': Epoch_mag_ds.min(), 'VALIDMAX': Epoch_mag_ds.max(), 'VAR_TYPE': 'support_data',
                           'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time',
                           'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear'}]}}


            varName = ['Bx', 'By', 'Bz']
            varData = np.array([B_ds[:, 0], B_ds[:, 1], B_ds[:, 2]])
            for i in range(3):
                data_dict = {**data_dict, **{varName[i]: [varData[i], {'LABLAXIS': varName[i],
                                                                                       'DEPEND_0': 'Epoch_mag',
                                                                                       'DEPEND_1': None,
                                                                                       'DEPEND_2': None,
                                                                                       'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                                       'FORMAT': 'E12.2',
                                                                                       'UNITS': 'nT',
                                                                                       'VALIDMIN': varData[i].min(),
                                                                                       'VALIDMAX': varData[i].max(),
                                                                                       'VAR_TYPE': 'data',
                                                                                       'SCALETYP': 'linear'}]}}

            varName = ['Bx_IGRF', 'By_IGRF', 'Bz_IGRF']
            varData = np.array([IGRF_B_rocket[:, 0], IGRF_B_rocket[:, 1], IGRF_B_rocket[:, 2]])

            for i in range(3):
                data_dict = {**data_dict, **{varName[i]: [varData[i], {'LABLAXIS': varName[i],
                                                              'DEPEND_0': 'Epoch_Attitude',
                                                              'DEPEND_1': None,
                                                              'DEPEND_2': None,
                                                              'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                              'UNITS': 'nT',
                                                              'VALIDMIN': varData[i].min(),
                                                              'VALIDMAX': varData[i].max(),
                                                              'VAR_TYPE': 'data',
                                                              'SCALETYP': 'linear'}]}}

            data_dict = {**data_dict, **{'Epoch_Attitude': [Epoch_attitude, {'DEPEND_0': 'Epoch_Attitude',
                                                                             'DEPEND_1': None, 'DEPEND_2': None,
                                                                             'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                             'FORMAT': 'I5', 'UNITS': 'ns',
                                                                             'VALIDMIN':None, 'VALIDMAX': None,
                                                                             'VAR_TYPE':'support_data',
                                                                             'MONOTON':'INCREASE', 'TIME_BASE':'J2000',
                                                                             'TIME_SCALE':'Terrestrial Time',
                                                                             'REFERENCE_POSITION':'Rotating Earth Geoid', 'SCALETYP':'linear'}]}}

            fileoutName = f'ACESII_{rocketID}_RingCore_IGRF_model.cdf'

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict, ModelData, globalAttrsMod, "RingCore")

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

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*Spun*')) == 0:
    print(color.RED + 'There are no RingCore_Spun.cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        RingCore_Despin(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        RingCore_Despin(wRocket, rocketFolderPath, justPrintFileNames,wflyer)