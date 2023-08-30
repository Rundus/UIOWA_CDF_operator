# --- RingCore_L1_to_L2_Despin&IGRF_interpolation.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Interpolate attitude data (alt,lat,long,DCM) onto ringCore Epoch and use
# it to Despin the ringCore data. The interpolation is done by drawing a line between
# successive points in the attitude data, seeing how many magnetometer epoch values fall
# between the two attitude epoch values, then evaluating the line at the magnetomter epoch values
# and storing the new attitude data.


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
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4
wFiles = [0]

inputPath_modifier = 'l1'
inputPath_modifier_attitude = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier_IGRF = 'science\IGRF_interpolated' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
outputPath_modifier_L2 = 'L2' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# --- constant offset over whole flight ---
determineScalarOffset = False
N = 20

offsetGridSearchRange = [[0.1049*1E9,0.1051*1E9], [0.105*1E9,0.106*1E9]] # + offsets applied to attitude solution timeseries. Only use positive offsets since ACS was always faster than MAG. Format: [linspace LOW, linspace HIGH]
# targetTimes = [[dt.datetime(2022,11,20,17,22,00,000000), dt.datetime(2022,11,20,17,30,00,000000)],
#                [dt.datetime(2022,11,20,17,23,45,000000), dt.datetime(2022,11,20,17,27,55,000000)]]# reduce data l2 to only regions of science interest

targetTimes = [[dt.datetime(2022,11,20,17,24,30,000000), dt.datetime(2022,11,20,17,25,20,000000)],
               [dt.datetime(2022,11,20,17,23,45,000000), dt.datetime(2022,11,20,17,27,55,000000)]]# reduce data l2 to only regions of science interest

# --- results of determineScalarOffset ---
offsetResults = [104944444.44444445, 117000000]

# HF results: [104944444.44444445, 390227194712994.5] format: [time in ns, BestChiSquare]
# LF results: [105789473.68421052, 1715649142403403.8]

# --- Dynamic De-spin Correction ---
dynamicDespin = True
plotRkTData_vs_IGRF = False # plots de-spun RKT data using the interpolated DCM with ONLY the scalar offset correction
lookAtFilteredB = True # need to ensure the filtered data ONLY removed the science wiggles. Here we can compare
plotIGRFvsNewlyDespun = False
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from pyIGRF import igrf_value
from numpy import array,matmul,dot,abs
from numpy.linalg import inv,norm
from scipy.signal import butter, filtfilt
from scipy.interpolate import CubicSpline

def butterworth(lowcutoff, highcutoff, fs, order, filtertype):
    if filtertype.lower() == 'bandpass':
        return butter(N = order, Wn= [lowcutoff, highcutoff], fs=fs, btype='bandpass')
    elif filtertype.lower() == 'lowpass':
        return butter(N=order, Wn = highcutoff, fs=fs, btype='lowpass')
    elif filtertype.lower() == 'highpass':
        return butter(N=order, Wn= lowcutoff, fs=fs, btype='highpass')
def butter_filter(data, lowcutoff, highcutoff, fs, order,filtertype):
    b, a = butterworth(lowcutoff, highcutoff, fs, order, filtertype)
    y = filtfilt(b, a, data)
    return y

def RingCore_Despin(wRocket,wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    ModelData = L0_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*RingCore*')
    inputFiles_attitude = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[wflyer]}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\\', '') for ifile in inputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]


    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Interpolating L1 data to DeSpun L2 data' + color.END)
        print('[0] ' + str(round(getsize(inputFiles[0]) / (10 ** 6), 1)) + 'MiB')

        ###########################
        # --- Get the Variables ---
        ###########################

        # --- get the Magnetometer Data ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict_mag = loadDictFromFile(inputFiles[wFile], {})
        data_dict_mag['Epoch'][0] = array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag['Epoch'][0]])
        Done(start_time)

        # --- get the Attitude Data ---
        prgMsg(f'Loading data from {inputPath_modifier_attitude} Files')
        data_dict_attitude = loadDictFromFile(inputFiles_attitude[0], {})
        data_dict_attitude['Epoch'][0] = array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude['Epoch'][0]])
        Done(start_time)

        # --- reduce dataset to science region ---
        prgMsg('Reducing data to science region')
        tTimes = [pycdf.lib.datetime_to_tt2000(targetTimes[wRocket - 4][0]),
                  pycdf.lib.datetime_to_tt2000(targetTimes[wRocket - 4][1])]

        lowCutoff, highCutoff = np.abs(data_dict_attitude['Epoch'][0] - tTimes[0]).argmin(), np.abs(data_dict_attitude['Epoch'][0] - tTimes[1]).argmin()
        for key, val in data_dict_attitude.items(): # apply reduction to attitude data
            if key == 'Alt':
                data_dict_attitude[key][0] = array(data_dict_attitude[key][0][lowCutoff:highCutoff])/1000 # convert to km
            else:
                data_dict_attitude[key][0] = array(data_dict_attitude[key][0][lowCutoff:highCutoff])

        lowCutoff, highCutoff = np.abs(data_dict_mag['Epoch'][0] - tTimes[0]).argmin(), np.abs(data_dict_mag['Epoch'][0] - tTimes[1]).argmin()
        for key, val in data_dict_mag.items(): # apply reduction to mag data
            data_dict_mag[key][0] = array(data_dict_mag[key][0][lowCutoff:highCutoff])

        Done(start_time)

        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ### INTERPOLATE ATTITUDE DATA ###
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        # --- [0] Collect Variables ---
        prgMsg('Collecting Variables')

        # define some storage variables
        dataKeys = ['Epoch', 'Alt', 'Latgd', 'Long', 'a11', 'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33']
        dataKeysVal = [deepcopy(data_dict_mag['Epoch'][0]), [], [], [], [], [], [], [], [], [], [], [], []]
        attitudeData = [deepcopy(data_dict_attitude[key][0]) for key in dataKeys] # a list to contain the attitude only the data that I care about
        dataInterp_dict_attitude = {key : value for key, value in zip(dataKeys, dataKeysVal)}

        # Define the Rocket Magnetic Field Vector
        B_rkt = array([[data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]] for i in range(len(data_dict_mag['Epoch'][0]))])

        if determineScalarOffset:

            ############################
            # --- FIND SCALAR OFFSET ---
            ############################
            prgMsg('Determining best scalar offset')
            offsetsVals = np.linspace(offsetGridSearchRange[wRocket-4][0],offsetGridSearchRange[wRocket-4][1], N)
            bestOffset = [0, 1E30] # format [offsetVal (in secs), ChiSquareVal]. Dummy values to be replaced

            # --- Loop over the temporal offsets ---
            for loopIndex in tqdm(range(len(offsetsVals))):

                # --- [1] apply offset to attitude data ---
                Epoch_attitude_loop = array([attitudeData[0][i] + offsetsVals[loopIndex] for i in range(len(attitudeData[0]))])

                # --- [2] Spline interpolate the time-adjusted attitude data onto the mag data ---
                counter = 0
                for key, newDataList in dataInterp_dict_attitude.items():
                    if key != 'Epoch':
                        # --- cubic interpolation ---
                        splCub = CubicSpline(Epoch_attitude_loop, attitudeData[counter])

                        # evaluate the interpolation at all the epoch_mag points
                        dataInterp_dict_attitude[key] = array([splCub(timeVal) for timeVal in dataInterp_dict_attitude['Epoch']])
                    counter += 1

                # EVALUATE SPUN-UP IGRF MODEL

                date = 2022 + 323 / 365  # Corresponds to 11/20/2022

                # --- [3] Evaluate IGRF at new attitude coordinates ---
                ### IGRF info ###
                # [3] North Comp (+ N | - S)
                # [4] East Comp (+ E | - W)
                # [5] Vertical Comp (+ D | - U)
                # [6] Total Field
                IGRF = array([igrf_value(dataInterp_dict_attitude['Latgd'][i], dataInterp_dict_attitude['Long'][i], dataInterp_dict_attitude['Alt'][i], date) for i in range(len(dataInterp_dict_attitude['Epoch']))])
                IGRF_ENU = array([ [vec[4], vec[3], -1 * vec[5]] for vec in IGRF])

                # --- [4] Despin the mag data ---
                DCM = array([
                    [[dataInterp_dict_attitude['a11'][i], dataInterp_dict_attitude['a12'][i], dataInterp_dict_attitude['a13'][i]],
                     [dataInterp_dict_attitude['a21'][i], dataInterp_dict_attitude['a22'][i], dataInterp_dict_attitude['a23'][i]],
                     [dataInterp_dict_attitude['a31'][i], dataInterp_dict_attitude['a32'][i], dataInterp_dict_attitude['a33'][i]]]
                    for i in range(len(dataInterp_dict_attitude['Epoch']))
                ]) # construct the new, interpolated DCM matrix

                # --- [5] get inverse of DCM, spin up IGRF ---
                IGRF_spun = array([matmul(inv(DCM[i]),IGRF_ENU[i]) for i in range(len(dataInterp_dict_attitude['Epoch']))])

                # --- [6] Analyze results of loop to see if chosen offset is good ---

                sub = IGRF_spun - B_rkt
                ChiSquare = (1 / len(dataInterp_dict_attitude['Epoch'])) *(sum([dot(sub[i], sub[i]) for i in range(len(sub))]) ** 2)
                print([offsetsVals[loopIndex], ChiSquare])
                if ChiSquare <= bestOffset[1]:
                    bestOffset = [offsetsVals[loopIndex], ChiSquare]

            Done(start_time)

            print(f'The Best Offset was: {bestOffset}')

        elif dynamicDespin:

            # --- [1] apply determined scalar offset to attitude data ---
            Epoch_attitude_offset = array([attitudeData[0][i] + offsetResults[wRocket-4] for i in range(len(attitudeData[0]))])

            # --- [2] Collect the spline interpolation information into a dictionary ---
            spline_dict = {}
            counter = 0

            for key, newDataList in dataInterp_dict_attitude.items():
                if key != 'Epoch':
                    # --- cubic interpolation ---
                    splCub = CubicSpline(Epoch_attitude_offset, attitudeData[counter])

                    # store the spline information for later
                    spline_dict = {**spline_dict, **{key: splCub}}

                    # evaluate the interpolation at all the epoch_mag points
                    dataInterp_dict_attitude[key] = array([splCub(timeVal) for timeVal in dataInterp_dict_attitude['Epoch']])

                counter += 1


            # apply the interpolation DCM
            interpDCM = array([[
                [spline_dict['a11'](tme), spline_dict['a12'](tme), spline_dict['a13'](tme)],
                [spline_dict['a21'](tme), spline_dict['a22'](tme), spline_dict['a23'](tme)],
                [spline_dict['a31'](tme), spline_dict['a32'](tme), spline_dict['a33'](tme)]
            ] for tme in dataInterp_dict_attitude['Epoch']])
            B_rkt_despun = array([matmul(interpDCM[i],B_rkt[i]) for i in range(len(B_rkt))])


            # --- [3] Low-pass Filter the "despun" data ---
            fs = 128
            order = 8
            cutoff = [3, 0]
            B_rkt_X_filtered = butter_filter(B_rkt_despun[:, 0], lowcutoff=cutoff[0], highcutoff=cutoff[1], fs=fs, order=order, filtertype = 'highpass')
            B_rkt_Y_filtered = butter_filter(B_rkt_despun[:, 1], lowcutoff=cutoff[0], highcutoff=cutoff[1], fs=fs, order=order, filtertype = 'highpass')
            B_rkt_Z_filtered = butter_filter(B_rkt_despun[:, 2], lowcutoff=cutoff[0], highcutoff=cutoff[1], fs=fs, order=order, filtertype = 'highpass')
            B_filtered = array([[B_rkt_X_filtered[i], B_rkt_Y_filtered[i], B_rkt_Z_filtered[i]] for i in range(len(B_rkt_Z_filtered))])

            if lookAtFilteredB:

                labels = ['X-axis', 'Y-axis', 'Z-axis']

                xDataPlot = array([pycdf.lib.tt2000_to_datetime(tme) for tme in dataInterp_dict_attitude['Epoch']])

                fig, ax = plt.subplots(3)
                ax[0].scatter(xDataPlot, B_filtered[:, 0])
                ax[0].set_ylabel('E-axis')
                ax[0].set_ylim(-20,20)
                ax[1].scatter(xDataPlot, B_filtered[:, 1])
                ax[1].set_ylabel('N-axis')
                ax[1].set_ylim(-20, 20)
                ax[2].scatter(xDataPlot, B_filtered[:, 2])
                ax[2].set_ylabel('U-axis')
                ax[2].set_ylim(-20, 20)
                plt.show()

                # fig, ax = plt.subplots(3)
                # ax[0].scatter(xDataPlot, B_rkt[:, 0])
                # ax[0].scatter(xDataPlot, B_filtered[:, 0])
                # ax[0].set_ylabel('X-axis')
                # ax[1].scatter(xDataPlot, B_rkt[:, 1])
                # ax[1].scatter(xDataPlot, B_filtered[:, 1])
                # ax[1].set_ylabel('Y-axis')
                # ax[2].scatter(xDataPlot, B_rkt[:, 2])
                # ax[2].scatter(xDataPlot, B_filtered[:, 2])
                # ax[2].set_ylabel('Z-axis')
                # plt.show()

            # --- [4] Loop through the interpolated attitude data to determine best DCM for each mag point ---
            # NOTE: The GOAL is to find a new DCM matrix for each mag point that minimizes the de-spin

            newDCM = []
            lookPercent = 1  # distance to look left/right for best MAG value
            M = 300 # number of points to look through between +/- deltaT on a single point
            date = 2022 + 323 / 365  # Corresponds to 11/20/2022

            # loop through all the interpolated points
            for i in tqdm(range(len(dataInterp_dict_attitude['Epoch']))):

                # determine the DeltaT around the mag point to check the interpolation

                if i == 0: # first index
                    deltaT = lookPercent*(dataInterp_dict_attitude['Epoch'][i+1] - dataInterp_dict_attitude['Epoch'][i+1])
                    deltaTspace = np.linspace(-1*deltaT, deltaT, M)
                elif i == len(dataInterp_dict_attitude['Epoch']) - 1: # last index
                    deltaT = lookPercent*(dataInterp_dict_attitude['Epoch'][i] - dataInterp_dict_attitude['Epoch'][i-1])
                    deltaTspace = np.linspace(-1 * deltaT, deltaT, M)
                else:
                    deltaT_p = lookPercent*(dataInterp_dict_attitude['Epoch'][i+1] - dataInterp_dict_attitude['Epoch'][i])
                    deltaT_n = lookPercent*(dataInterp_dict_attitude['Epoch'][i] - dataInterp_dict_attitude['Epoch'][i-1])
                    deltaTspace = np.linspace(-1 * deltaT_n, deltaT_p, M)

                # evaluate all the interpolated DCM around all the deltaT points and apply the DCM to the FILTERED rocket B-data
                bestDCM = [array([[1, 2, 3], [1, 2, 3], [1, 2, 3]]), 1E20]
                for timeVal in deltaTspace:
                    evalTime = dataInterp_dict_attitude['Epoch'][i] + timeVal

                    searchingDCM = array([
                    [spline_dict['a11'](evalTime), spline_dict['a12'](evalTime), spline_dict['a13'](evalTime)],
                    [spline_dict['a21'](evalTime), spline_dict['a22'](evalTime), spline_dict['a23'](evalTime)],
                    [spline_dict['a31'](evalTime), spline_dict['a32'](evalTime), spline_dict['a33'](evalTime)]
                    ])

                    IGRF = array(igrf_value(dataInterp_dict_attitude['Latgd'][i], dataInterp_dict_attitude['Long'][i], dataInterp_dict_attitude['Alt'][i], date))

                    ringCoreIGRFresidual = norm(matmul(searchingDCM, B_rkt[i]) - array([IGRF[4], IGRF[3], -1 * IGRF[5]]))

                    # for each points, find the DCM for which there's a minimum in IGRF - DCM*rkt
                    if ringCoreIGRFresidual < bestDCM[1]:
                        bestDCM = [searchingDCM, ringCoreIGRFresidual]


                # store this new DCM as the DCM for the i-th Epoch_mag point
                newDCM.append(bestDCM[0])


            # --- --- --- --- --- --- --- ---
            # --- PREPARE DATA FOR OUTPUT ---
            # --- --- --- --- --- --- --- ---
            # apply the new DCM to the ringcore data
            B_rkt_despun = array([matmul(newDCM[i], B_rkt[i]) for i in range(len(B_rkt))])

            # get the IGRF field at the magnetometer cadence
            date = 2022 + 323 / 365  # Corresponds to 11/20/2022
            IGRF = array(
                [
                    igrf_value(data_dict_attitude['Latgd'][0][i], data_dict_attitude['Long'][0][i],
                               data_dict_attitude['Alt'][0][i], date) for i in
                    range(len(Epoch_attitude_offset))
                ]
            )
            IGRF_ENU = array([[vec[4], vec[3], -1 * vec[5]] for vec in IGRF])

            # spin-up the IGRF using the interpolated DCM but NOT the "optimized" DCM
            DCM_interpNotOptimized = array([[
                [spline_dict['a11'](tme), spline_dict['a12'](tme), spline_dict['a13'](tme)],
                [spline_dict['a21'](tme), spline_dict['a22'](tme), spline_dict['a23'](tme)],
                [spline_dict['a31'](tme), spline_dict['a32'](tme), spline_dict['a33'](tme)]
            ] for tme in Epoch_attitude_offset])

            IGRF_spun = array([matmul(inv(DCM_interpNotOptimized[i]), IGRF_ENU[i]) for i in range(len(Epoch_attitude_offset))])



        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:

            prgMsg('Creating output file')

            # items I need for diangostics
            # [1] rocket data that's shifted by the scalar amount
            # [2] spun IGRF field
            # [3] Attempted RKT Despun Data

            # # construct magnetometer data B-vector and despin the magnetometer data
            # B_ENU = matmul(array([data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]]), IGRF_ENU)

            IGRF_east, IGRF_north, IGRF_up = array(IGRF_ENU[:, 0]), array(IGRF_ENU[:, 1]), array(IGRF_ENU[:, 2])
            IGRF_rkt_X, IGRF_rkt_Y, IGRF_rkt_Z = array(IGRF_spun[:, 0]), array(IGRF_spun[:, 1]), array(IGRF_spun[:, 2])
            IGRF_mag = array([norm(IGRF_ENU[i]) for i in range(len(IGRF_ENU))])

            def buildDataDict(labels, dataToOutput):
                data_dict = {}
                for i, data in enumerate(dataToOutput):
                    if labels[i] != 'Epoch': # non-epoch case
                        data_dict = {**data_dict, **{labels[i]:[data, {'LABLAXIS': labels[i], 'DEPEND_0': 'Epoch',
                                                                  'DEPEND_1': None,
                                                                  'DEPEND_2': None,
                                                                  'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                                  'UNITS': 'nT',
                                                                  'VALIDMIN': data.min(),
                                                                  'VALIDMAX': data.max(),
                                                                  'VAR_TYPE': 'data',
                                                                  'SCALETYP': 'linear'}]}}
                    else: # Epoch case

                        # data = array([datum + offset[wRocket-4] for datum in data])
                        data_dict = {**data_dict, **{labels[i]: [data, {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None,
                                                                                 'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                                 'FORMAT': 'I5', 'UNITS': 'ns',
                                                                                 'VALIDMIN':data.min(), 'VALIDMAX': data.max(),
                                                                                 'VAR_TYPE':'support_data',
                                                                                 'MONOTON':'INCREASE', 'TIME_BASE':'J2000',
                                                                                 'TIME_SCALE':'Terrestrial Time',
                                                                                 'REFERENCE_POSITION':'Rotating Earth Geoid', 'SCALETYP':'linear'}]}}
                return data_dict

            ########################
            ### IGRF OUTPUT DATA ###
            ########################

            # --- Construct the data dict for IGRF ---
            labels = ['Epoch', 'IGRF_east', 'IGRF_north', 'IGRF_up', 'IGRF_X', 'IGRF_Y', 'IGRF_Z', 'IGRF_mag']
            dataToOutput = [array(Epoch_attitude_offset), IGRF_east, IGRF_north, IGRF_up, IGRF_rkt_X, IGRF_rkt_Y, IGRF_rkt_Z, IGRF_mag]
            data_dict_IGRF = buildDataDict(labels, dataToOutput)

            fileoutName = f'ACESII_{rocketID}_IGRF_interpolated_model.cdf'
            outputPath = f'{rocketFolderPath}{outputPath_modifier_IGRF}\{fliers[wflyer]}\\{fileoutName}'
            outputCDFdata(outputPath, data_dict_IGRF, ModelData, globalAttrsMod, "RingCore")

            # -- Construct the data dict for despun mag data ---
            labels = ['Epoch', 'B_east', 'B_north', 'B_up', 'B_x', 'B_y', 'B_z']
            dataToOutput = [dataInterp_dict_attitude['Epoch'], B_rkt_despun[:, 0], B_rkt_despun[:, 1], B_rkt_despun[:, 2], B_rkt[:, 0], B_rkt[:, 1], B_rkt[:, 2]]
            data_dict_Despun = buildDataDict(labels, dataToOutput)

            fileoutName = f'ACESII_{rocketID}_l2_RingCore_Despun.cdf'
            outputPath = f'{rocketFolderPath}{outputPath_modifier_L2}\{fliers[wflyer]}\\{fileoutName}'
            outputCDFdata(outputPath, data_dict_Despun, ModelData, globalAttrsMod, "RingCore")

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

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*_RingCore*')) == 0:
    print(color.RED + 'There are no RingCore.cdf files in the specified directory' + color.END)
else:
    for filesNo in wFiles:
        RingCore_Despin(wRocket, filesNo, rocketFolderPath, justPrintFileNames, wflyer)
