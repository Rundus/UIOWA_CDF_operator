#--- L1_to_Langmuir.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert the engineering Langmuir data to scientifically useful units.
# Also renames "Boom_Monitor_1/2" --> "Fixed_Boom_Monitor_1/2" etc
# it was discovered that ni_swept for the low AND high Flyer start with an extra value that should not be there. This is fixed in L0_to_L1.py


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

#################
# --- TOGGLES ---
#################

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

# D:\Data\\ACESII\\science\\\low\ACESII_36364_langmuir_playback_BSS_Card0_processedtm_20221120T170824.cdf

###############################
# --- FIXED VOLTAGE TOGGLES ---
###############################
unit_conversion = 10**(9) # 1 for Amps 10**9 for nano amps, etc

# Ideal parameters for ACES High and Low [no_of_splits=2,splice_value = 2]
fixedProbeCal = True
no_of_splits = 2 # 0, 1 or 2
split_value = 3 # If no_of_splits == 1: calibrated curve data is split in TWO, starting at the <--- split_value, If no_of_splits == 2: calibrated curve data is split in THREE, starting at the <--- splice value, then +1 to +2
splice_value = 2
rounding = 20 # displays the error parameters up to this many decimal points

##################################
# --- SWEPT TO VOLTAGE TOGGLES ---
##################################
# Does a little analysis on the "step" variable to determine the average values of the steps.
# Needed to calculate from "Step" variable to voltage
# MUST BE ==TRUE WHEN SWEPT PROBE ==TRUE
stepToVoltage = True

# the capacitive effect between the probe and plasma cause an RC decay on the data.
# This toggle only uses the 10th (the last) datapoint of each voltage setpoint to eliminate this effect
downSample_RCeffect = True

#################################
# --- SWEPT PROBE CAL TOGGLES ---
#################################
sweptProbeCal_Andoya = True
applySweptCalCurve = False # False - data is in analog values, True - data is in current
useSingleSweeptData = True # use data from a single good sweep or multiple sweeps and average. THIS GIVES BETTER DATA IF TRUE
plotAllCalCurves = False
plotmasterData = False
plotFittedCurves = False

# the thresholds of the different regions for BOTH flyers
threshAnalog = [[-3960, 1850], [-3900, 2000]]  # Break curve up first by analog value into 2 regions
threshCurrent = [[-2, 1], [-5, 0]]  # separate data by current value, in NANOAMPS


#################################
# --- SWEPT PROBE CAL TOGGLES ---
#################################
sweptProbeCal_Iowa = True

####################################
# --- SWEPT PROBE ROCKET TOGGLES ---
####################################
sweptProbeData = True
breakIntoCurves = True
targetVoltage_min = -1 # only care about voltage sweeps above this voltage value. Nominally -1
digitalVariance = 5 # how much the digitized step point can vary when looking for the top and bottom of curves. nominally = 5
indvEpochThresh = 15000000 # Value, in tt2000, that determines the time diff needed between epoch points to identify particular sweeps



#####################
# --- DATA OUTPUT ---
#####################
outputData = True

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

def calFunction_fixed(x, A, B, C):
    y = -A * ((x) ** (B)) - C
    return y

def linear(x, A, B):
    y = A * x + B
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
    LangmuirSweptCalFiles = glob(f'{rocketFolderPath}\science\LP_calibration\{fliers[wflyer]}\*_345deg_*')

    L1_names = [ifile.replace(f'{rocketFolderPath}L1\{fliers[wflyer]}\\', '') for ifile in L1Files]

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
                        if key not in data_dict_cal:
                            data_dict_cal = {**data_dict_cal, **{key: [LangmuirSweptCalFiles[key][...], {key: val for key, val in LangmuirSweptCalFiles[key].attrs.items()}]}}

        Done(start_time)

        #####################
        # --- Fixed Probe ---
        #####################

        if fixedProbeCal:
            plotFixedCalCurve = False

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
                parameters, covariance = scipy.optimize.curve_fit(calFunction_fixed, digital_vals[k],calibrationCurrents[k],maxfev=100000)
                fit_params.append(parameters)

                if plotFixedCalCurve:
                    linData = [calFunction_fixed(testData[k][i],parameters[0],parameters[1],parameters[2]) for i in range(len(testData[k]))]
                    plt.scatter(digital_vals[k], calibrationCurrents[k])
                    plt.plot(testData[k],linData)

                # Quantify Fit errors
                fixedCalErrors = {
                    'digital_value       ':[digital_vals[k][i] for i in range(len(digital_vals[k]))],
                    'calCurrent          ':[calibrationCurrents[k][i] for i in range(len(digital_vals[k]))],
                    'calFunc             ':[round(calFunction_fixed(digital_vals[k][i],parameters[0],parameters[1],parameters[2]),rounding) for i in range(len(digital_vals[k]))],
                    'Difference          ':[round(calibrationCurrents[k][i] - calFunction_fixed(digital_vals[k][i],parameters[0],parameters[1],parameters[2]),rounding) for i in range(len(digital_vals[k]))],
                    'calCurrent % diff.  ':[str(round(100*(calibrationCurrents[k][i] - calFunction_fixed(digital_vals[k][i],parameters[0],parameters[1],parameters[2]))/calibrationCurrents[k][i],rounding))+'%' for i in range(len(digital_vals[k]))]
                }

                if plotFixedCalCurve:
                    print(f'--- CALIBRATION ERRORS FOR FIXED, LOOP NO. {k + 1}:---')

                    for key, val in fixedCalErrors.items():
                        print(key, val)

                    print('FIT PARAMETERS for A*(x**(B)) - C', parameters, '\n')

            if plotFixedCalCurve:
                plt.xlabel('Digital')
                plt.ylabel(f'Current [nA])')
                plt.show()

            # Apply the calibration function curve
            ni = np.zeros(shape = (len(data_dict['ni'][0])))
            sign_of_data = -1 # better to keep the current as positive. People can understand what to do as long as it's labelled ni
            if no_of_loops == 1:
                for j in range(len(data_dict['ni'][0])):
                    ni[j] = sign_of_data*calFunction_fixed(data_dict['ni'][0][j],fit_params[0][0],fit_params[0][1],fit_params[0][2])
            elif no_of_loops == 2:
                for j in range(len(data_dict['ni'][0])):
                    if data_dict['ni'][0][j] <= max(digital_vals[0]):
                        ni[j] = sign_of_data*calFunction_fixed(data_dict['ni'][0][j],fit_params[0][0],fit_params[0][1],fit_params[0][2])
                    elif data_dict['ni'][0][j] > max(digital_vals[0]):
                        ni[j] = sign_of_data*calFunction_fixed(data_dict['ni'][0][j],fit_params[1][0],fit_params[1][1],fit_params[1][2])
            elif no_of_loops == 3:
                for j in range(len(data_dict['ni'][0])):
                    if data_dict['ni'][0][j] <= max(digital_vals[0]):
                        ni[j] = sign_of_data*calFunction_fixed(data_dict['ni'][0][j],fit_params[0][0],fit_params[0][1],fit_params[0][2])
                    elif data_dict['ni'][0][j] > max(digital_vals[0]) and data_dict['ni'][0][j] <= max(digital_vals[1]):
                        ni[j] = sign_of_data*calFunction_fixed(data_dict['ni'][0][j],fit_params[1][0],fit_params[1][1],fit_params[1][2])
                    elif data_dict['ni'][0][j] > max(digital_vals[1]):
                        ni[j] = sign_of_data*calFunction_fixed(data_dict['ni'][0][j], fit_params[2][0], fit_params[2][1], fit_params[2][2])


            for i,data in enumerate(ni):
                if np.abs(data)>=5000:
                    ni[i] = -1e30


            data_dict = {**data_dict, **{'fixed_ni': [ni, {'LABLAXIS': 'ni',
                                                            'DEPEND_0': 'fixed_Epoch',
                                                            'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': -1e30,
                                                            'FORMAT': 'E12.2',
                                                            'UNITS': 'nA',
                                                            'VALIDMIN': ni.min(),
                                                            'VALIDMAX': ni.max(),
                                                            'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
            data_dict['fixed_Epoch'] = data_dict.pop('Epoch_ni') # rename to fixed
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


            ###############################################
            # --- REMOVE RC EFFECT BY /10  DOWNSAMPLING ---
            ###############################################
            if downSample_RCeffect:

                # Find where the bottom of the first sweep occurs:
                stepStart = np.abs(data_dict['step'][0] - list(counted_reduced.keys())[0]).argmin()
                n_i_swept = np.array(data_dict['ni_swept'][0][stepStart::])
                n_e_swept = np.array(data_dict['ne_swept'][0][stepStart::])
                sweptStep = np.array(data_dict['step'][0][stepStart::])
                Epoch_sweptCurrent = np.array(data_dict['Epoch_step'][0][stepStart::])
                n_i_swept_Div10, n_e_Swept_Div10, sweptStep_Div10, Epoch_sweptCurrent_Div10 = [], [], [], []

                # Downsample the data
                for i in range(0, len(sweptStep)-10, 10):
                    n_i_swept_Div10.append(n_i_swept[i+9])
                    n_e_Swept_Div10.append(n_e_swept[i+9])
                    sweptStep_Div10.append(sweptStep[i+9])
                    Epoch_sweptCurrent_Div10.append(Epoch_sweptCurrent[i+9])

                data_dict['ni_swept'][0] = np.array(n_i_swept_Div10)
                data_dict['ne_swept'][0] = np.array(n_e_Swept_Div10)
                data_dict['step'][0] = np.array(sweptStep_Div10)
                data_dict['Epoch_step'][0] = np.array(Epoch_sweptCurrent_Div10)

            # --- linear conversion ---
            voltageRange = rocketAttrs.LPswept_voltage_range[wRocket - 4]
            finalStep_digital = max(counted_reduced)
            initialStep_digital = min(counted_reduced)
            slope = (voltageRange[1] - voltageRange[0]) / (finalStep_digital - initialStep_digital)
            intercept = voltageRange[0]

            def step_to_Voltage(analog_voltage):
                return slope*analog_voltage + intercept

            stepVoltage = np.array([slope*data_dict['step'][0][i] + intercept if data_dict['step'][0][i] not in [-1,65535] else data_dict['step'][0][i] for i in range(len(data_dict['step'][0]))])
            data_dict = {**data_dict, **{'step_Voltage': [stepVoltage, {'LABLAXIS': 'step Voltage', 'DEPEND_0': 'Epoch_step', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1e30, 'FORMAT': 'E12.2', 'UNITS': 'Volts', 'VALIDMIN': stepVoltage.min(), 'VALIDMAX': stepVoltage.max(), 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}

            Done(start_time)

        ###########################################
        # --- Swept Probe Calibration at Andoya ---
        ###########################################

        if sweptProbeCal_Andoya:

            prgMsg('Collecting Swept Calibration Data')

            if useSingleSweeptData:
                sweptCalRanges = rocketAttrs.LPswept_cal_epoch_ranges_single_sweep[wflyer]
            else:
                sweptCalRanges = rocketAttrs.LPswept_cal_epoch_ranges[wflyer]

            # --- --- --- --- --- --- --- --- ---
            # --- COLLECT/CALC CAL RANGE DATA ---
            # --- --- --- --- --- --- --- --- ---
            calEpoch = np.array([ pycdf.lib.datetime_to_tt2000(data_dict_cal['Epoch_ne_swept'][0][i]) for i in range(len(data_dict_cal['Epoch_ne_swept'][0]))])

            sweptCal_voltage = [] # voltage of the step
            sweptCal_calCurrent = [] # corresponding current based on voltage
            sweptCal_analog_current = [] # current determined from the analag calibration data
            sweptCal_Epoch = [] # epoch values of the start/end of this cal range
            for i in range(len(sweptCalRanges)):

                # find the indicies of the calibration range
                start = np.abs(calEpoch - pycdf.lib.datetime_to_tt2000(sweptCalRanges[i][0])).argmin()
                end = np.abs(calEpoch - pycdf.lib.datetime_to_tt2000(sweptCalRanges[i][1])).argmin()

                # calculate the analog current, voltage and calibration current
                sweptCal_analog_current.append(np.array(data_dict_cal['ne_swept'][0][start:end]) - np.array(data_dict_cal['ni_swept'][0][start:end]))
                sweptCal_Epoch.append(calEpoch[start:end])
                resistance = rocketAttrs.LPswept_cal_resistances[i]
                sweptCal_step = data_dict_cal['step'][0][start:end]
                sweptCal_voltage.append(np.array([step_to_Voltage(sweptCal_step[i]) for i in range(len(sweptCal_step))]))
                sweptCal_calCurrent.append(sweptCal_voltage[-1]/resistance)

            Done(start_time)

            # --- --- --- --- --- --- --- --- ---
            # --- Remove outliers from data ---
            # --- --- --- --- --- --- --- --- ---

            #surprisingly, only the high flyer cal needed needs to be filtered for outliers
            if wRocket == 4:
                removeOutliers = True
            elif wRocket == 5:
                removeOutliers = False

            if removeOutliers:
                prgMsg('Removing Outliers')
                threshold_percent = 0.05 # if the next point is >= 50% the the distance between max/min of the dataset, it is an outlier
                repeat_process = 20 # how many times to repeat outlier search

                for process in range(repeat_process):
                    for i in range(len(sweptCalRanges)):
                        maxV = sweptCal_analog_current[i].max()
                        minV = sweptCal_analog_current[i].min()
                        threshDistance = np.abs(maxV - minV) * threshold_percent

                        # check every value and record the indices
                        remove_indices = []

                        for j in range(len(sweptCal_analog_current[i]) - 1):
                            # take an average of data to see where I am:
                            no_of_points = 10
                            AVG = sum([sweptCal_analog_current[i][j - k] for k in range(no_of_points) ])/no_of_points

                            if AVG < -3400: #apply a stricter threshold for things close to when the current is low
                                if np.abs(sweptCal_analog_current[i][j+1] - sweptCal_analog_current[i][j]) >= threshDistance:
                                    remove_indices.append(j+1)
                            elif np.abs(sweptCal_analog_current[i][j+1]) >= 4095: # remove any points at absolute maximum
                                remove_indices.append(j + 1)
                            elif sweptCal_Epoch[i][j+1] < 3178576.8709884263:
                                remove_indices.append(j+1)

                        sweptCal_calCurrent[i] = np.delete(sweptCal_calCurrent[i], remove_indices, axis=0)
                        sweptCal_voltage[i] = np.delete(sweptCal_voltage[i], remove_indices, axis=0)
                        sweptCal_analog_current[i] = np.delete(sweptCal_analog_current[i], remove_indices, axis=0)
                        sweptCal_Epoch[i] = np.delete(sweptCal_Epoch[i], remove_indices, axis=0)

                Done(start_time)

            # PLOT ALL THE RAW CAL CURVES
            if plotAllCalCurves:
                for i in range(1, 7):
                    plt.scatter(sweptCal_analog_current[i], sweptCal_calCurrent[i], label=str(rocketAttrs.LPswept_cal_resistances[i]))

                plt.legend()
                plt.show()

            #################################################
            # --- Create MASTER calibration curve dataset ---
            #################################################

            prgMsg('Creating Master Dataset')

            # Find where the min/max in the calibration current occurs
            minCurrent_calCur = [sweptCal_calCurrent[i].min() for i in range(len(sweptCal_analog_current))]
            maxCurrent_calCur = [sweptCal_calCurrent[i].max() for i in range(len(sweptCal_analog_current))]

            # NOTE: We will not be including the "open" case during calibration.
            sweptCal_master_analog_current_temp = [  [[], []] for i in range(7)]
            sweptCal_master_calCurrent_temp = [  [[], []] for i in range(7)]

            for i in range(1, len(sweptCal_analog_current)): # loop through all curves EXCEPT the "open" case

                # loop through all values in a sweep
                for j in range(len(sweptCal_analog_current[i])): # loop through all values within a curve
                    value = sweptCal_calCurrent[i][j]

                    if i ==1: # special case for the final set of data. Collect all data here but break it in half
                        if j <= int(len(sweptCal_analog_current[i])/2):
                            sweptCal_master_calCurrent_temp[i][0].append(value)
                            sweptCal_master_analog_current_temp[i][0].append(sweptCal_analog_current[i][j])
                        elif j > int(len(sweptCal_analog_current[i])/2):
                            sweptCal_master_calCurrent_temp[i][1].append(value)
                            sweptCal_master_analog_current_temp[i][1].append(sweptCal_analog_current[i][j])
                    else:
                        # only collect data from a curve if you in the special range
                        if (minCurrent_calCur[i] <= value < minCurrent_calCur[i-1]): # Minimum range
                            sweptCal_master_calCurrent_temp[i][0].append(value)
                            sweptCal_master_analog_current_temp[i][0].append(sweptCal_analog_current[i][j])

                        elif (maxCurrent_calCur[i-1] < value <= maxCurrent_calCur[i]): #Maximum range
                            sweptCal_master_calCurrent_temp[i][1].append(value)
                            sweptCal_master_analog_current_temp[i][1].append(sweptCal_analog_current[i][j])

            # PLOT THE MASTER DATA
            if plotmasterData:
                colors = ['brown', 'red', 'gold', 'orange', 'green', 'blue', 'purple']
                for i in range(1,7):
                    for j in range(2):
                        plt.scatter(sweptCal_master_analog_current_temp[i][j], sweptCal_master_calCurrent_temp[i][j], color = colors[i], label=f"{rocketAttrs.LPswept_cal_resistances[i]/(10**6)} M")

                # plt.plot(xData_fit, yData_fit )
                plt.xlabel('$I_{analog}$')
                plt.ylabel('$I_{cal}$ [nA]')
                plt.xscale('symlog')
                plt.show()

            # Flatten the master data into one dataset
            sortVals = [[6,0],[5,0],[4,0],[3,0],[2,0],[1,0],[1,1],[2,1],[3,1],[4,1],[5,1],[6,1]]

            sweptCal_master_analog_current = []
            sweptCal_master_calCurrent = []

            # Put the data into the master dataset in the order from lowest resistance to highest, then flip at middlepoint
            for val in sortVals:
                sweptCal_master_analog_current.append(sweptCal_master_analog_current_temp[val[0]][val[1]])
                sweptCal_master_calCurrent.append(sweptCal_master_calCurrent_temp[val[0]][val[1]])

            # Flatten all the data AND convert it to nanoAmperes
            sweptCal_master_analog_current = np.array([val for sublist in sweptCal_master_analog_current for val in sublist])
            sweptCal_master_calCurrent = np.array([val*unit_conversion for sublist in sweptCal_master_calCurrent for val in sublist])

            Done(start_time)


            ##############################
            # --- FIT THE MASTER CURVE ---
            ##############################

            prgMsg('Calculating Calibration curve')

            bottomData = [[],[]]
            middleData= [[],[]]
            topData = [[],[]]

            #Cut up the master data into pieces to fit curves to
            for i in range(len(sweptCal_master_calCurrent)):
                xVal = sweptCal_master_analog_current[i]
                yVal = sweptCal_master_calCurrent[i]

                # TOP DATA, analog >= 2000
                if xVal >= threshAnalog[wRocket-4][1]:
                    if yVal >= threshCurrent[wRocket-4][1]: # TOP, analog 2000+
                        topData[0].append(xVal)
                        topData[1].append(yVal)
                    elif yVal < threshCurrent[wRocket-4][1]: # Middle TOP, analog 2000+
                        middleData[0].append(xVal)
                        middleData[1].append(yVal)

                # MIDDLE TOP, -4000 < analog < 2000
                elif threshAnalog[wRocket-4][0] < xVal < threshAnalog[wRocket-4][1]:
                    if yVal >= threshCurrent[wRocket - 4][0]:  # Middle BOTTOM, -4095 < analog < 0, close to 0 current
                        middleData[0].append(xVal)
                        middleData[1].append(yVal)
                    elif yVal < threshCurrent[wRocket - 4][0]:  # BOTTOM, -4095 < analog < 0 and Negative current
                        bottomData[0].append(xVal)
                        bottomData[1].append(yVal)

            N = 1000

            # Middle data DOWN
            paramsMid, covarianceMid = scipy.optimize.curve_fit(linear, middleData[0], middleData[1], maxfev=100000, p0=[1.5 / 6000, 0.5])
            xData_mid = np.linspace(min(middleData[0]), max(middleData[0]), N)
            fittedMidData = [linear(x, paramsMid[0], paramsMid[1]) for x in xData_mid]

            if plotFittedCurves:

                # --- PLOTTING ---
                rounded = 5
                fig, ax = plt.subplots(2,1)
                ax[0].set_title('$I_{analog}$ vs $I_{cal}$ fitted plots, y = mx+b ')

                ax[0].scatter(middleData[0],middleData[1])
                ax[0].plot(xData_mid,fittedMidData,color='red',label=f"m: {round(paramsMid[0],rounded)}   b: {round(paramsMid[1],rounded)}")
                ax[0].legend()

                ax[1].scatter(middleData[0], middleData[1])
                ax[1].plot(xData_mid, fittedMidData, color='red',label=f"m: {round(paramsMid[0], rounded)}   b: {round(paramsMid[1], rounded)}")
                ax[1].legend()
                ax[1].set_xscale('symlog')


                plt.tight_layout()
                plt.show()


            # Apply Calibration Curve
            def sweptCal_Analog_to_Current(xVal):
                # TOP DATA, -4000 <= analog <= 2000
                if threshAnalog[wRocket - 4][1] <= xVal <= threshAnalog[wRocket - 4][1]:
                    return linear(xVal, paramsMid[0], paramsMid[1])
                else:
                    return linear(xVal, paramsMid[0], paramsMid[1])


            Done(start_time)

        ###########################################
        # --- Swept Probe Calibration at Iowa ---
        ###########################################

        if sweptProbeCal_Iowa:

            prgMsg('Collecting Swept Calibration Data')

            if useSingleSweeptData:
                sweptCalRanges = rocketAttrs.LPswept_cal_epoch_ranges_single_sweep[wflyer]
            else:
                sweptCalRanges = rocketAttrs.LPswept_cal_epoch_ranges[wflyer]

            # --- --- --- --- --- --- --- --- --- ---
            # --- COLLECT/PROCESS CALIBRATION DATA ---
            # --- --- --- --- --- --- --- --- --- ---
            calEpoch = np.array([pycdf.lib.datetime_to_tt2000(data_dict_cal['Epoch_ne_swept'][0][i]) for i in
                                 range(len(data_dict_cal['Epoch_ne_swept'][0]))])

            sweptCal_voltage = []  # voltage of the step
            sweptCal_calCurrent = []  # corresponding current based on voltage
            sweptCal_analog_current = []  # current determined from the analag calibration data
            sweptCal_Epoch = []  # epoch values of the start/end of this cal range
            for i in range(len(sweptCalRanges)):
                # find the indicies of the calibration range
                start = np.abs(calEpoch - pycdf.lib.datetime_to_tt2000(sweptCalRanges[i][0])).argmin()
                end = np.abs(calEpoch - pycdf.lib.datetime_to_tt2000(sweptCalRanges[i][1])).argmin()

                # calculate the analog current, voltage and calibration current
                sweptCal_analog_current.append(np.array(data_dict_cal['ne_swept'][0][start:end]) - np.array(
                    data_dict_cal['ni_swept'][0][start:end]))
                sweptCal_Epoch.append(calEpoch[start:end])
                resistance = rocketAttrs.LPswept_cal_resistances[i]
                sweptCal_step = data_dict_cal['step'][0][start:end]
                sweptCal_voltage.append(
                    np.array([step_to_Voltage(sweptCal_step[i]) for i in range(len(sweptCal_step))]))
                sweptCal_calCurrent.append(sweptCal_voltage[-1] / resistance)

            Done(start_time)

            # --- --- --- --- --- --- --- --- ---
            # --- Remove outliers from data ---
            # --- --- --- --- --- --- --- --- ---

            # surprisingly, only the high flyer cal needed needs to be filtered for outliers
            if wRocket == 4:
                removeOutliers = True
            elif wRocket == 5:
                removeOutliers = False

            if removeOutliers:
                prgMsg('Removing Outliers')
                threshold_percent = 0.05  # if the next point is >= 50% the the distance between max/min of the dataset, it is an outlier
                repeat_process = 20  # how many times to repeat outlier search

                for process in range(repeat_process):
                    for i in range(len(sweptCalRanges)):
                        maxV = sweptCal_analog_current[i].max()
                        minV = sweptCal_analog_current[i].min()
                        threshDistance = np.abs(maxV - minV) * threshold_percent

                        # check every value and record the indices
                        remove_indices = []

                        for j in range(len(sweptCal_analog_current[i]) - 1):
                            # take an average of data to see where I am:
                            no_of_points = 10
                            AVG = sum(
                                [sweptCal_analog_current[i][j - k] for k in range(no_of_points)]) / no_of_points

                            if AVG < -3400:  # apply a stricter threshold for things close to when the current is low
                                if np.abs(sweptCal_analog_current[i][j + 1] - sweptCal_analog_current[i][
                                    j]) >= threshDistance:
                                    remove_indices.append(j + 1)
                            elif np.abs(sweptCal_analog_current[i][
                                            j + 1]) >= 4095:  # remove any points at absolute maximum
                                remove_indices.append(j + 1)
                            elif sweptCal_Epoch[i][j + 1] < 3178576.8709884263:
                                remove_indices.append(j + 1)

                        sweptCal_calCurrent[i] = np.delete(sweptCal_calCurrent[i], remove_indices, axis=0)
                        sweptCal_voltage[i] = np.delete(sweptCal_voltage[i], remove_indices, axis=0)
                        sweptCal_analog_current[i] = np.delete(sweptCal_analog_current[i], remove_indices, axis=0)
                        sweptCal_Epoch[i] = np.delete(sweptCal_Epoch[i], remove_indices, axis=0)

                Done(start_time)

            # PLOT ALL THE RAW CAL CURVES
            if plotAllCalCurves:
                for i in range(1, 7):
                    plt.scatter(sweptCal_analog_current[i], sweptCal_calCurrent[i],
                                label=str(rocketAttrs.LPswept_cal_resistances[i]))

                plt.legend()
                plt.show()

            #################################################
            # --- Create MASTER calibration curve dataset ---
            #################################################

            prgMsg('Creating Master Dataset')

            # Find where the min/max in the calibration current occurs
            minCurrent_calCur = [sweptCal_calCurrent[i].min() for i in range(len(sweptCal_analog_current))]
            maxCurrent_calCur = [sweptCal_calCurrent[i].max() for i in range(len(sweptCal_analog_current))]

            # NOTE: We will not be including the "open" case during calibration.
            sweptCal_master_analog_current_temp = [[[], []] for i in range(7)]
            sweptCal_master_calCurrent_temp = [[[], []] for i in range(7)]

            for i in range(1, len(sweptCal_analog_current)):  # loop through all curves EXCEPT the "open" case

                # loop through all values in a sweep
                for j in range(len(sweptCal_analog_current[i])):  # loop through all values within a curve
                    value = sweptCal_calCurrent[i][j]

                    if i == 1:  # special case for the final set of data. Collect all data here but break it in half
                        if j <= int(len(sweptCal_analog_current[i]) / 2):
                            sweptCal_master_calCurrent_temp[i][0].append(value)
                            sweptCal_master_analog_current_temp[i][0].append(sweptCal_analog_current[i][j])
                        elif j > int(len(sweptCal_analog_current[i]) / 2):
                            sweptCal_master_calCurrent_temp[i][1].append(value)
                            sweptCal_master_analog_current_temp[i][1].append(sweptCal_analog_current[i][j])
                    else:
                        # only collect data from a curve if you in the special range
                        if (minCurrent_calCur[i] <= value < minCurrent_calCur[i - 1]):  # Minimum range
                            sweptCal_master_calCurrent_temp[i][0].append(value)
                            sweptCal_master_analog_current_temp[i][0].append(sweptCal_analog_current[i][j])

                        elif (maxCurrent_calCur[i - 1] < value <= maxCurrent_calCur[i]):  # Maximum range
                            sweptCal_master_calCurrent_temp[i][1].append(value)
                            sweptCal_master_analog_current_temp[i][1].append(sweptCal_analog_current[i][j])

            # PLOT THE MASTER DATA
            if plotmasterData:
                colors = ['brown', 'red', 'gold', 'orange', 'green', 'blue', 'purple']
                for i in range(1, 7):
                    for j in range(2):
                        plt.scatter(sweptCal_master_analog_current_temp[i][j],
                                    sweptCal_master_calCurrent_temp[i][j], color=colors[i],
                                    label=f"{rocketAttrs.LPswept_cal_resistances[i] / (10 ** 6)} M")

                # plt.plot(xData_fit, yData_fit )
                plt.xlabel('$I_{analog}$')
                plt.ylabel('$I_{cal}$ [nA]')
                plt.xscale('symlog')
                plt.show()

            # Flatten the master data into one dataset
            sortVals = [[6, 0], [5, 0], [4, 0], [3, 0], [2, 0], [1, 0], [1, 1], [2, 1], [3, 1], [4, 1], [5, 1],
                        [6, 1]]

            sweptCal_master_analog_current = []
            sweptCal_master_calCurrent = []

            # Put the data into the master dataset in the order from lowest resistance to highest, then flip at middlepoint
            for val in sortVals:
                sweptCal_master_analog_current.append(sweptCal_master_analog_current_temp[val[0]][val[1]])
                sweptCal_master_calCurrent.append(sweptCal_master_calCurrent_temp[val[0]][val[1]])

            # Flatten all the data AND convert it to nanoAmperes
            sweptCal_master_analog_current = np.array(
                [val for sublist in sweptCal_master_analog_current for val in sublist])
            sweptCal_master_calCurrent = np.array(
                [val * unit_conversion for sublist in sweptCal_master_calCurrent for val in sublist])

            Done(start_time)

            ##############################
            # --- FIT THE MASTER CURVE ---
            ##############################

            prgMsg('Calculating Calibration curve')

            bottomData = [[], []]
            middleData = [[], []]
            topData = [[], []]

            # Cut up the master data into pieces to fit curves to
            for i in range(len(sweptCal_master_calCurrent)):
                xVal = sweptCal_master_analog_current[i]
                yVal = sweptCal_master_calCurrent[i]

                # TOP DATA, analog >= 2000
                if xVal >= threshAnalog[wRocket - 4][1]:
                    if yVal >= threshCurrent[wRocket - 4][1]:  # TOP, analog 2000+
                        topData[0].append(xVal)
                        topData[1].append(yVal)
                    elif yVal < threshCurrent[wRocket - 4][1]:  # Middle TOP, analog 2000+
                        middleData[0].append(xVal)
                        middleData[1].append(yVal)

                # MIDDLE TOP, -4000 < analog < 2000
                elif threshAnalog[wRocket - 4][0] < xVal < threshAnalog[wRocket - 4][1]:
                    if yVal >= threshCurrent[wRocket - 4][
                        0]:  # Middle BOTTOM, -4095 < analog < 0, close to 0 current
                        middleData[0].append(xVal)
                        middleData[1].append(yVal)
                    elif yVal < threshCurrent[wRocket - 4][0]:  # BOTTOM, -4095 < analog < 0 and Negative current
                        bottomData[0].append(xVal)
                        bottomData[1].append(yVal)

            N = 1000

            # Middle data DOWN
            paramsMid, covarianceMid = scipy.optimize.curve_fit(linear, middleData[0], middleData[1], maxfev=100000,
                                                                p0=[1.5 / 6000, 0.5])
            xData_mid = np.linspace(min(middleData[0]), max(middleData[0]), N)
            fittedMidData = [linear(x, paramsMid[0], paramsMid[1]) for x in xData_mid]

            if plotFittedCurves:
                # --- PLOTTING ---
                rounded = 5
                fig, ax = plt.subplots(2, 1)
                ax[0].set_title('$I_{analog}$ vs $I_{cal}$ fitted plots, y = mx+b ')

                ax[0].scatter(middleData[0], middleData[1])
                ax[0].plot(xData_mid, fittedMidData, color='red',
                           label=f"m: {round(paramsMid[0], rounded)}   b: {round(paramsMid[1], rounded)}")
                ax[0].legend()

                ax[1].scatter(middleData[0], middleData[1])
                ax[1].plot(xData_mid, fittedMidData, color='red',
                           label=f"m: {round(paramsMid[0], rounded)}   b: {round(paramsMid[1], rounded)}")
                ax[1].legend()
                ax[1].set_xscale('symlog')

                plt.tight_layout()
                plt.show()

            # Apply Calibration Curve
            def sweptCal_Analog_to_Current(xVal):
                # TOP DATA, -4000 <= analog <= 2000
                if threshAnalog[wRocket - 4][1] <= xVal <= threshAnalog[wRocket - 4][1]:
                    return linear(xVal, paramsMid[0], paramsMid[1])
                else:
                    return linear(xVal, paramsMid[0], paramsMid[1])

            Done(start_time)

        #####################
        # --- Swept Probe ---
        #####################

        if sweptProbeData:
            prgMsg('Calculating swept LP current')

            # Calculate the swept current
            n_i_swept = np.array(data_dict['ni_swept'][0])
            n_e_swept = np.array(data_dict['ne_swept'][0])
            sweptStep_temp = np.array(data_dict['step_Voltage'][0])
            Epoch_sweptCurrent_temp = np.array(data_dict['Epoch_step'][0])
            sweptCurrent_temp = np.zeros(shape=(len(data_dict['ne_swept'][0])))

            for i in range(len(n_i_swept)):
                sweptCurrent_temp[i] = (n_e_swept[i] - n_i_swept[i])

            if breakIntoCurves:
                Epoch_sweptCurrent = []
                sweptCurrent = []
                step_sweptVoltage = []

                # Reduce the sweptCurrent data to only include data with stepVoltage >= targetVoltage_min
                for i in range(len(data_dict['step_Voltage'][0])):
                    if data_dict['step_Voltage'][0][i] >= targetVoltage_min:
                        sweptCurrent.append(sweptCurrent_temp[i])
                        Epoch_sweptCurrent.append(Epoch_sweptCurrent_temp[i])
                        step_sweptVoltage.append(sweptStep_temp[i])

                # Reduce data to only specified epoch range that contains good data
                epochIndexLow = np.abs(np.array(Epoch_sweptCurrent) - rocketAttrs.startEndLangmuirBreakIntoCurves[wRocket - 4][0]).argmin()
                epochIndexHigh = np.abs(np.array(Epoch_sweptCurrent) - rocketAttrs.startEndLangmuirBreakIntoCurves[wRocket - 4][1]).argmin()
                sweptCurrent = sweptCurrent[epochIndexLow:epochIndexHigh+1]
                step_sweptVoltage = step_sweptVoltage[epochIndexLow:epochIndexHigh+1]
                Epoch_sweptCurrent = np.array([pycdf.lib.datetime_to_tt2000(Epoch_sweptCurrent[i]) for i in range(epochIndexLow,epochIndexHigh+1)])

                # Identify Individual sweeps first by the large gap in Epoch value
                indvSweepIndicies = []
                sweepNo = 0

                for i in range(len(Epoch_sweptCurrent)-1):

                    if (Epoch_sweptCurrent[i+1] - Epoch_sweptCurrent[i]) > indvEpochThresh:

                        if sweepNo == 0:
                            indvSweepIndicies.append([0, i])
                        else:
                            indvSweepIndicies.append([indvSweepIndicies[-1][1]+1, i])

                        sweepNo += 1

                # the Above algorithim doesn't include the final sweep, so I add it in here
                indvSweepIndicies.append([indvSweepIndicies[-1][1]+1, len(Epoch_sweptCurrent)])

                # Separate each individual sweep into upleg and downleg
                sweptCurrent_New = []
                Epoch_sweptCurrent_New = []
                sweptStep_New = []
                for index in range(len(indvSweepIndicies)):
                    indvCurrent = sweptCurrent[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]
                    indvEpoch = Epoch_sweptCurrent[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]
                    indvStep = step_sweptVoltage[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]

                    if downSample_RCeffect:
                        # Take the top value of "step" and split the curve based on that
                        middleIndex = np.array(indvStep).argmax()

                    else:
                        # take the top 9 values of 'indvCurernt' and split the curve at the middle value
                        countedIndicies = sorted(range(len(indvCurrent)), key=lambda i: indvCurrent[i])[-9:]
                        countedIndicies.sort()
                        middleIndex = countedIndicies[int((len(countedIndicies) - 1) / 2)]

                    # Break up the curve
                    currentUpLeg = indvCurrent[0:middleIndex]
                    epochUpLeg = indvEpoch[0:middleIndex]
                    stepUpLeg = indvStep[0:middleIndex]

                    currentDownLeg = indvCurrent[middleIndex:]
                    epochDownLeg = indvEpoch[middleIndex:]
                    stepDownLeg = indvStep[middleIndex:]

                    # Store the broken up curve data
                    sweptCurrent_New.append(currentUpLeg)
                    sweptCurrent_New.append(currentDownLeg[::-1])
                    Epoch_sweptCurrent_New.append(epochUpLeg)
                    Epoch_sweptCurrent_New.append(epochDownLeg)
                    sweptStep_New.append(stepUpLeg)
                    sweptStep_New.append(stepDownLeg[::-1])

                sweptCurrent = []
                Epoch_sweptCurrent = []
                step_sweptVoltage = []

                # Flatten the data. np.flatten does not work for some reason
                for i in range(len(sweptCurrent_New)):
                    for thing in sweptCurrent_New[i]:
                        sweptCurrent.append(thing)
                    for thing in Epoch_sweptCurrent_New[i]:
                        Epoch_sweptCurrent.append(thing)
                    for thing in sweptStep_New[i]:
                        step_sweptVoltage.append(thing)


                sweptCurrent = np.array(sweptCurrent)
                Epoch_sweptCurrent = np.array(Epoch_sweptCurrent)
                step_sweptVoltage = np.array(step_sweptVoltage)

                ###########################################################
                # --- APPLY THE CALIBRATION CURVES FROM SWEPT_PROBE_CAL ---
                ###########################################################
                if applySweptCalCurve:
                    sweptCurrent = np.array([sweptCal_Analog_to_Current(analogVal) for analogVal in sweptCurrent])

            else:
                sweptCurrent = np.array(sweptCurrent_temp)
                Epoch_sweptCurrent = np.array(Epoch_sweptCurrent_temp)
                step_sweptVoltage = np.array(sweptStep_temp)

                if applySweptCalCurve:
                    sweptCurrent = np.array([sweptCal_Analog_to_Current(analogVal) for analogVal in sweptCurrent])
                else:
                    sweptCurrent = sweptCurrent_temp


            units = 'nA' if applySweptCalCurve else 'Analog'

            data_dict = {**data_dict, **{'swept_Current': [sweptCurrent, {'LABLAXIS': 'swept_Current',
                                                            'DEPEND_0': 'Epoch_swept_Current',
                                                            'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                            'UNITS': units,
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

            data_dict = {**data_dict, **{'swept_Voltage': [step_sweptVoltage, {'LABLAXIS': 'swept_Voltage',
                                                                          'DEPEND_0': 'Epoch_swept_Current',
                                                                          'DEPEND_1': None,
                                                                          'DEPEND_2': None,
                                                                          'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                                          'UNITS': 'Volts',
                                                                          'VALIDMIN': step_sweptVoltage.min(),
                                                                          'VALIDMAX': step_sweptVoltage.max(),
                                                                          'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}

            Done(start_time)

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        if outputData:

            # remove uneeded data.
            removeThese = ['ni', 'ni_swept', 'ne_swept', "Epoch_ni_swept", 'Epoch_ne_swept', 'step', 'EXP_Current', '28V_Monitor', 'Epoch_monitors']
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
            outputPath = f'{outputFolderPath}\\' + f'{fliers[wRocket-4]}\{fileoutName}'

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