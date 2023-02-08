#--- L1_to_Langmuir.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert the engineering Langmuir data to scientifically useful units.
# Also renames "Boom_Monitor_1/2" --> "Fixed_Boom_Monitor_1/2" etc
# it was discovered that ni_swept for the low AND high Flyer start with an extra value that should not be there. Fixed in L0_to_L1.py


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

# --- FIXED VOLTAGE TOGGLES ---
unit_conversion = 10**(9) # 1 for Amps 10**9 for nano amps, etc
unitLabel = 'nA'
# IDeal parameters for ACES High and Low [no_of_splits=2,splice_value = 2]
fixedDigitalToVoltage = False
no_of_splits = 2 # 0, 1 or 2
split_value = 3 # If no_of_splits == 1: calibrated curve data is split in TWO, starting at the <--- split_value, If no_of_splits == 2: calibrated curve data is split in THREE, starting at the <--- splice value, then +1 to +2
splice_value = 2
rounding = 20 # displays the error parameters up to this many decimal points
plotCalCurve = False

# --- SWEPT TO VOLTAGE TOGGLES ---
# Does a little analysis on the "step" variable to determine the average values of the steps.
# Needed to calculate from "Step" variable to voltage
# MUST BE TRUE WHEN SWEPT PROBE == TRUE
stepToVoltage = False

# --- SWEPT PROBE CAL TOGGLES ---
sweptProbeCal = True

# --- SWEPT PROBE ROCKET TOGGLES ---
sweptProbeRocketData = False
breakIntoIndCurves = True
targetVoltage_min = -1 # only care about voltage sweeps above this voltage value. Nominally -1
digitalVariance = 5 # how much the digitized step point can vary when looking for the top and bottom of curves. nominally = 5
indvEpochThresh = 15000000 # Value, in tt2000, that determines the time diff needed between epoch points to identify particular sweeps

# --- DATA OUTPUT ---
outputData = False


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy

from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, L2_ACES_Quick,L2_TRICE_Quick, prgMsg
from glob import glob
from os.path import getsize
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)
from collections import Counter



def calFunction(x, A, B, C):
    y = -A * ((x) ** (B)) - C
    return y

def L1_to_Langmuir(wRocket, rocketFolderPath, justPrintFileNames, wflyer):

    outputFolderPath = rocketFolderPath + r'L2\\'

    # --- Get ACESII rocket Attributes ---
    rocketAttrs,b,c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'Langmuir'
    L2ModelData = L2_ACES_Quick(wflyer)

    # Set the paths for the file names
    L1Files = glob(f'{rocketFolderPath}L1\{fliers[wflyer]}\*_lp_*')
    LangmuirFiles = glob(f'{outputFolderPath}\{fliers[wflyer]}\*_langmuir_*')
    LangmuirSweptCalFiles = glob(f'{rocketFolderPath}\science\LP_calibration\{fliers[wflyer]}\*_lp_*')

    L1_names = [ifile.replace(f'{rocketFolderPath}L1\{fliers[wflyer]}\\', '') for ifile in L1Files]
    Langmuir_names = [ofile.replace(f'{outputFolderPath}\{fliers[wflyer]}\\', '') for ofile in LangmuirFiles]

    L1_names_searchable = [ifile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace('l1_', '').replace('_v00', '') for ifile in L1_names]

    dataFile_name = L1_names_searchable[0].replace(f'{rocketFolderPath}L1\{fliers[wflyer]}\\', '').replace('lp_','').replace('ni_','').replace('ne_swept_','').replace('step_','').replace('ni_swept_','').replace('deltaNdivN_','')
    fileoutName = rf'ACESII_{rocketAttrs.rocketID[wRocket-4]}_langmuir_{dataFile_name}'

    if justPrintFileNames:
            for i, file in enumerate(L1Files):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, L1_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Creating langmuir data for {fliers[wRocket-4]} flyer' + color.END)

        # --- get the data from the tmCDF file ---
        prgMsg('Loading data from L1Files')

        # Collect the LP data except deltaNdivN into a data dict
        data_dict = {}
        for file in L1Files:
            if 'deltaNdivN' not in file:
                with pycdf.CDF(file) as L1DataFile:
                    for key, val in L1DataFile.items():
                        if key not in data_dict:
                            data_dict = {**data_dict, **{key : [L1DataFile[key][...] , {key:val for key,val in L1DataFile[key].attrs.items()  }  ]  }  }

        # Collect the LP data except deltaNdivN into a data dict
        data_dict_cal = {}
        for file in LangmuirSweptCalFiles:
            if 'deltaNdivN' not in file:
                with pycdf.CDF(file) as LangmuirSweptCalFiles:
                    for key, val in LangmuirSweptCalFiles.items():
                        if key not in data_dict:
                            data_dict_cal = {**data_dict, **{key: [LangmuirSweptCalFiles[key][...], {key: val for key, val in LangmuirSweptCalFiles[key].attrs.items()}]}}

        Done(start_time)
        print(LangmuirSweptCalFiles)
        for key,val in data_dict_cal.items():
            print(key)

        #####################
        # --- Fixed Probe ---
        #####################

        if fixedDigitalToVoltage:
            prgMsg('Converting fixed probe to Voltage')
            fixedCalResistances = rocketAttrs.LPFixed_calResistances[wRocket - 4]
            probeBias = rocketAttrs.LPFixedProbeBias[wRocket - 4]
            calibrationCurrents = []
            digital_vals = []
            for key, val in fixedCalResistances.items():
                digital_vals.append(val)
                if key == 'Open':
                    calibrationCurrents.append(probeBias/(10*unit_conversion)) #10 gigaOhms
                else:
                    calibrationCurrents.append((unit_conversion)*probeBias/key)

            digital_vals = np.array(digital_vals)
            calibrationCurrents = np.array(calibrationCurrents)

            if no_of_splits == 2:
                no_of_loops = 3
                digital_vals = [np.array(digital_vals[0:splice_value+1]),
                                np.array(digital_vals[splice_value:splice_value+3]),
                                np.array(digital_vals[splice_value+2:])]
                calibrationCurrents = [np.array(calibrationCurrents[0:splice_value+1]),
                                       np.array(calibrationCurrents[splice_value:splice_value+3]),
                                       np.array(calibrationCurrents[splice_value+2:])]
                testData = [np.linspace(0, max(digital_vals[0]), 1000),
                            np.linspace(max(digital_vals[0]), max(digital_vals[1]), 1000),
                            np.linspace(max(digital_vals[1]), 4095 , 1000)]
                fit_params = []
            elif no_of_splits == 1:
                no_of_loops = 2
                digital_vals = [np.array(digital_vals[0:split_value]),np.array(digital_vals[split_value:])]
                calibrationCurrents = [np.array(calibrationCurrents[0:split_value]),np.array(calibrationCurrents[split_value:])]
                testData = [np.linspace(0, max(digital_vals[0]), 1000),np.linspace(max(digital_vals[0]), 4095 , 1000)]
                fit_params = []
            elif no_of_splits == 0:
                no_of_loops = 1
                digital_vals = [np.array(digital_vals)]
                calibrationCurrents = [np.array(calibrationCurrents)]
                testData = [np.linspace(0, 4095, 4096)]
                fit_params = []

            for k in range(no_of_loops):
                parameters, covariance = scipy.optimize.curve_fit(calFunction, digital_vals[k],calibrationCurrents[k],maxfev=100000)
                fit_params.append(parameters)

                if plotCalCurve:
                    linData = [calFunction(testData[k][i],parameters[0],parameters[1],parameters[2]) for i in range(len(testData[k]))]
                    plt.scatter(digital_vals[k], calibrationCurrents[k])
                    plt.plot(testData[k],linData)

                # Quantify Fit errors
                fixedCalErrors = {
                    'digital_value       ':[digital_vals[k][i] for i in range(len(digital_vals[k]))],
                    'calCurrent          ':[calibrationCurrents[k][i] for i in range(len(digital_vals[k]))],
                    'calFunc             ':[round(calFunction(digital_vals[k][i],parameters[0],parameters[1],parameters[2]),rounding) for i in range(len(digital_vals[k]))],
                    'Difference          ':[round(calibrationCurrents[k][i] - calFunction(digital_vals[k][i],parameters[0],parameters[1],parameters[2]),rounding) for i in range(len(digital_vals[k]))],
                    'calCurrent % diff.  ':[str(round(100*(calibrationCurrents[k][i] - calFunction(digital_vals[k][i],parameters[0],parameters[1],parameters[2]))/calibrationCurrents[k][i],rounding))+'%' for i in range(len(digital_vals[k]))]
                }

                if plotCalCurve:
                    print(f'--- CALIBRATION ERRORS FOR FIXED, LOOP NO. {k + 1}:---')

                    for key, val in fixedCalErrors.items():
                        print(key, val)

                    print('FIT PARAMETERS for A*(x**(B)) - C', parameters, '\n')

            if plotCalCurve:
                plt.xlabel('Digital')
                plt.ylabel(f'Current ({unitLabel})')
                plt.show()

            # Apply the calibration function curve
            ni = np.zeros(shape = (len(data_dict['ni'][0])))
            sign_of_data = -1 # better to keep the current as positive. People can understand what to do as long as it's labelled ni
            if no_of_loops == 1:
                for j in range(len(data_dict['ni'][0])):
                    ni[j] = sign_of_data*calFunction(data_dict['ni'][0][j],fit_params[0][0],fit_params[0][1],fit_params[0][2])
            elif no_of_loops == 2:
                for j in range(len(data_dict['ni'][0])):
                    if data_dict['ni'][0][j] <= max(digital_vals[0]):
                        ni[j] = sign_of_data*calFunction(data_dict['ni'][0][j],fit_params[0][0],fit_params[0][1],fit_params[0][2])
                    elif data_dict['ni'][0][j] > max(digital_vals[0]):
                        ni[j] = sign_of_data*calFunction(data_dict['ni'][0][j],fit_params[1][0],fit_params[1][1],fit_params[1][2])
            elif no_of_loops == 3:
                for j in range(len(data_dict['ni'][0])):
                    if data_dict['ni'][0][j] <= max(digital_vals[0]):
                        ni[j] = sign_of_data*calFunction(data_dict['ni'][0][j],fit_params[0][0],fit_params[0][1],fit_params[0][2])
                    elif data_dict['ni'][0][j] > max(digital_vals[0]) and data_dict['ni'][0][j] <= max(digital_vals[1]):
                        ni[j] = sign_of_data*calFunction(data_dict['ni'][0][j],fit_params[1][0],fit_params[1][1],fit_params[1][2])
                    elif data_dict['ni'][0][j] > max(digital_vals[1]):
                        ni[j] = sign_of_data*calFunction(data_dict['ni'][0][j], fit_params[2][0], fit_params[2][1], fit_params[2][2])


            data_dict = {**data_dict, **{'n_i': [ni, {'LABLAXIS': 'ni',
                                                            'DEPEND_0': 'Epoch_ni',
                                                            'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': -1e30,
                                                            'FORMAT': 'E12.2',
                                                            'UNITS': unitLabel,
                                                            'VALIDMIN': ni.min(),
                                                            'VALIDMAX': ni.max(),
                                                            'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
            Done(start_time)

        #########################
        # --- STEP to Voltage ---
        #########################

        if stepToVoltage:

            prgMsg('Calculating swept step voltage')

            # calculate the epoch index of beginning/end of sample range
            Epoch_step = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_step'][0][i]) for i in range(len(data_dict['Epoch_step'][0]))])
            sampleStart = np.abs(np.array(Epoch_step - pycdf.lib.datetime_to_tt2000(rocketAttrs.Epoch_range_to_determine_stepDAC[wRocket-4][0]))).argmin()
            sampleEnd = np.abs(np.array(Epoch_step - pycdf.lib.datetime_to_tt2000(rocketAttrs.Epoch_range_to_determine_stepDAC[wRocket-4][1]))).argmin()


            adjustments = [[2,2],[10,11]]

            # determine the values of step for each step
            sampleData = data_dict['step'][0][sampleStart+adjustments[wRocket-4][0]:sampleEnd+adjustments[wRocket-4][1]]
            stepsDigitalVals = [round(sum(sampleData[i*10:(i+1)*10])/10) for i in range(round(len(sampleData)/10))]

            counted = Counter(stepsDigitalVals)
            counted_reduced = dict(sorted(counted.items()))

            # filter out values of the step so that only 101 remain
            loopDict = [counted_reduced for i in range(4)]
            for j in range(len(loopDict)):
                keys = []
                removekeys = []

                for key, val in loopDict[j].items():
                    keys.append(key)

                for i in range(len(keys)-1):
                    if np.abs(keys[i+1] - keys[i]) < 9:
                        if loopDict[j][keys[i+1]] > loopDict[j][keys[i]]:
                            removekeys.append(keys[i])
                        elif loopDict[j][keys[i+1]] < loopDict[j][keys[i]]:
                            removekeys.append(keys[i+1])
                        elif loopDict[j][keys[i]] == loopDict[j][keys[i]]:
                            removekeys.append(keys[i+1])
                    if loopDict[j][keys[i]] < 5:
                        removekeys.append(keys[i])

                for thing in removekeys:
                    counted_reduced.pop(thing,None)

            voltageRange = rocketAttrs.LPswept_voltage_range[wRocket - 4]
            finalStep_digital = max(counted_reduced)
            initialStep_digital = min(counted_reduced)

            # --- linear conversion ---
            slope = (voltageRange[1] - voltageRange[0]) / (finalStep_digital - initialStep_digital)
            intercept = voltageRange[0]

            def step_to_Voltage(analog_voltage):
                return slope*analog_voltage + intercept


            stepVoltage = np.array([slope*data_dict['step'][0][i] + intercept if data_dict['step'][0][i] not in [-1,65535] else data_dict['step'][0][i] for i in range(len(data_dict['step'][0]))])
            data_dict = {**data_dict, **{'step_Voltage': [stepVoltage, {'LABLAXIS': 'Step Voltage',
                                                         'DEPEND_0': 'Epoch_step',
                                                         'DEPEND_1': None,
                                                         'DEPEND_2': None,
                                                         'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                         'UNITS': 'Volts',
                                                         'VALIDMIN': stepVoltage.min(), 'VALIDMAX': stepVoltage.max(),
                                                         'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}


            Done(start_time)

        #################################
        # --- Swept Probe Calibration ---
        #################################

        if sweptProbeCal:
            print('potato')

        #####################
        # --- Swept Probe ---
        #####################

        if sweptProbeRocketData:
            prgMsg('Calculating swept LP current')

            # Calculate the swept current
            n_i_swept = np.array(data_dict['ni_swept'][0])
            n_e_swept = np.array(data_dict['ne_swept'][0])
            Epoch_sweptCurrent_temp = np.array(data_dict['Epoch_step'][0])
            sweptCurrent_temp = np.zeros(shape=(len(data_dict['ne_swept'][0])))

            for i in range(len(n_i_swept)):
                sweptCurrent_temp[i] = (n_e_swept[i] - n_i_swept[i])/unit_conversion
                if np.abs(sweptCurrent_temp[i]) > 10:
                    sweptCurrent_temp[i] = 0


            if breakIntoIndCurves:
                Epoch_sweptCurrent = []
                sweptCurrent = []

                # Reduce the sweptCurrent data to only include data with stepVoltage >= targetVoltage_min
                for i in range(len(data_dict['step_Voltage'][0])):
                    if data_dict['step_Voltage'][0][i] >= targetVoltage_min:
                        sweptCurrent.append(sweptCurrent_temp[i])
                        Epoch_sweptCurrent.append(Epoch_sweptCurrent_temp[i])

                # Reduce data to only specified epoch range that contains good data
                epochIndexLow = np.abs(np.array(Epoch_sweptCurrent) - rocketAttrs.startEndLangmuirBreakIntoCurves[wRocket - 4][0]).argmin()
                epochIndexHigh = np.abs(np.array(Epoch_sweptCurrent) - rocketAttrs.startEndLangmuirBreakIntoCurves[wRocket - 4][1]).argmin()
                sweptCurrent = sweptCurrent[epochIndexLow:epochIndexHigh+1]
                Epoch_sweptCurrent = np.array([pycdf.lib.datetime_to_tt2000(Epoch_sweptCurrent[i]) for i in range(epochIndexLow,epochIndexHigh+1)])

                # Identify Individual sweeps first by the large gap in Epoch value
                indvSweepIndicies = []
                sweepNo = 0

                for i in range(len(Epoch_sweptCurrent)-1):

                    if (Epoch_sweptCurrent[i+1] - Epoch_sweptCurrent[i]) > indvEpochThresh:

                        if sweepNo == 0:
                            indvSweepIndicies.append([0,i])
                        else:
                            indvSweepIndicies.append([indvSweepIndicies[-1][1]+1, i])

                        sweepNo += 1

                # the Above algorithim doesn't include the final sweep, so I add it in here
                indvSweepIndicies.append([indvSweepIndicies[-1][1]+1, len(Epoch_sweptCurrent)])

                # Separate each individual sweep into upleg and downleg
                sweptCurrent_New = []
                Epoch_sweptCurrent_New = []
                for index in range(len(indvSweepIndicies)):
                    indvCurrent = sweptCurrent[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]
                    indvEpoch = Epoch_sweptCurrent[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]

                    # take the top 9 values and split the curve based on that
                    countedIndicies = sorted(range(len(indvCurrent)), key=lambda i: indvCurrent[i])[-9:]
                    countedIndicies.sort()
                    middleIndex = countedIndicies[int((len(countedIndicies) - 1) / 2)]

                    # Break up the curve
                    currentUpLeg = indvCurrent[0:middleIndex]
                    epochUpLeg = indvEpoch[0:middleIndex]

                    currentDownLeg = indvCurrent[middleIndex:]
                    epochDownLeg = indvEpoch[middleIndex:]

                    # Store the broken up curve data
                    sweptCurrent_New.append(currentUpLeg)
                    sweptCurrent_New.append(currentDownLeg[::-1])
                    Epoch_sweptCurrent_New.append(epochUpLeg)
                    Epoch_sweptCurrent_New.append(epochDownLeg)

                if stepToVoltage:
                    # Separate the up-leg and down-leg parts of the curves
                    step_values_Voltage = np.array([slope * data + intercept for data in counted_reduced])
                    curveBottom = step_values_Voltage[np.abs(step_values_Voltage- targetVoltage_min).argmin()]
                    curveTop = step_values_Voltage.max()
                    stepVoltageVariance = slope*digitalVariance


                sweptCurrent = []
                Epoch_sweptCurrent = []

                # Flatten the data. np.flatten does not work for some reason
                for i in range(len(sweptCurrent_New)):
                    for thing in sweptCurrent_New[i]:
                        sweptCurrent.append(thing)
                    for thing in Epoch_sweptCurrent_New[i]:
                        Epoch_sweptCurrent.append(thing)

                sweptCurrent = np.array(sweptCurrent)
                Epoch_sweptCurrent = np.array(Epoch_sweptCurrent)



            else:
                sweptCurrent = np.array(sweptCurrent_temp)
                Epoch_sweptCurrent = np.array(Epoch_sweptCurrent_temp)

            data_dict = {**data_dict, **{'swept_Current': [sweptCurrent, {'LABLAXIS': 'swept_Current',
                                                            'DEPEND_0': 'Epoch_swept_Current',
                                                            'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                            'UNITS': unitLabel,
                                                            'VALIDMIN': sweptCurrent.min(),
                                                            'VALIDMAX': sweptCurrent.max(),
                                                            'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
            data_dict = {**data_dict, **{'Epoch_swept_Current': [Epoch_sweptCurrent, {'LABLAXIS':'Epoch_swept_Current',
                                                                    'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None,
                                                                    'FILLVAL': -9223372036854775808, 'FORMAT': 'E12.2', 'UNITS': 'ns',
                                                                    'VALIDMIN': Epoch_sweptCurrent.min(), 'VALIDMAX': Epoch_sweptCurrent.max(), 'VAR_TYPE': 'support_data',
                                                                    'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000',
                                                                    'TIME_SCALE': 'Terrestrial Time','REFERENCE_POSITION': 'Rotating Earth Geoid',
                                                                    'SCALETYP': 'linear'}]}}
            Done(start_time)



        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        if outputData:

            # --- small processing ---

            # remove uneeded data
            removeThese = ['ni','ni_swept','ne_swept',"Epoch_ni_swept",'Epoch_ne_swept','step','EXP_Current','28V_Monitor','Epoch_monitors','Epoch_ni']
            for thing in removeThese:
                data_dict.pop(thing)

            # rename Boom_Monitors
            data_dict['Swept_Boom_Monitor'] = data_dict.pop('Boom_Monitor_1')
            data_dict['Swept_Boom_Monitor'][1]['DEPEND_0'] = "Epoch_Swept_Boom_Monitor"
            data_dict['Epoch_Swept_Boom_Monitor'] = data_dict.pop('Epoch_monitor_1')

            data_dict['Fixed_Boom_Monitor'] = data_dict.pop('Boom_Monitor_2')
            data_dict['Fixed_Boom_Monitor'][1]['DEPEND_0'] = "Epoch_Fixed_Boom_Monitor"
            data_dict['Epoch_Fixed_Boom_Monitor'] = data_dict.pop('Epoch_monitor_2')

            prgMsg('Creating output file')
            outputPath = f'{outputFolderPath}\\' + f'{fliers[wflyer]}\{fileoutName}'

            # --- delete output file if it already exists ---
            if os.path.exists(outputPath):
                os.remove(outputPath)

            # --- open the output file ---
            with pycdf.CDF(outputPath, '') as LangmuirFile:
                LangmuirFile.readonly(False)

                # --- write out global attributes ---
                inputGlobDic = L2ModelData.cdfFile.globalattsget()
                for key, val in inputGlobDic.items():
                    if key in globalAttrsMod:
                        LangmuirFile.attrs[key] = globalAttrsMod[key]
                    else:
                        LangmuirFile.attrs[key] = val

                # --- WRITE OUT DATA ---
                for varKey, varVal in data_dict.items():
                    if 'Epoch' in varKey: # epoch data
                        LangmuirFile.new(varKey, data=varVal[0], type=33)
                    else: # other data
                        LangmuirFile.new(varKey, data=varVal[0],type=pycdf.const.CDF_REAL8)

                    # --- Write out the attributes and variable info ---
                    for attrKey, attrVal in data_dict[varKey][1].items():
                        if attrKey == 'VALIDMIN':
                            LangmuirFile[varKey].attrs[attrKey] = varVal[0].min()
                        elif attrKey == 'VALIDMAX':
                            LangmuirFile[varKey].attrs[attrKey] = varVal[0].max()
                        elif attrVal != None:
                            LangmuirFile[varKey].attrs[attrKey] = attrVal

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

if len(glob(f'{rocketFolderPath}L1\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        L1_to_Langmuir(wRocket, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        L1_to_Langmuir(wRocket, rocketFolderPath, justPrintFileNames,wflyer)