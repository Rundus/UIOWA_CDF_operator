# --- L0_to_L1.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert electrostatic analyzer data from sector counts to raw counts
# Truncates all data to only that past 17:20:00.
# Special case given to LP's s.t. its dataset starts on sfid = 0 when reducing dataset


# NOTE: I have aligned the sector counts to start on an sfid 0 and end with an sfid 39 BUT
#  that doesn't mean I've started my energy sweeps at the beginning. It takes 4 minorframes
#  to complete 1 energy record and there are 49 + 1 retrace energy records per instrument sweep.
#  It takes 360,280 and 360 words for the eepaa, iepaa, and leesa to complete a sweep respectively,
#  which means for each major frame eepaa,iepaa,leesa complete 7 sweeps with 10 words left over



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import gc
# --- --- --- --- ---
import time
from ACESII_code.class_var_func import Done, setupPYCDF
start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files in
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
wFiles = [2]

# how many energy values not to keep, starting from the lowest values
energy_adjustment = 1

#########################
# --- Special Toggles ---
#########################
# Truncates all data to everything past 17:20:00
truncateData = True



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os
import datetime as dt

from copy import deepcopy
from gc import collect
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, L1_TRICE_Quick, prgMsg
from glob import glob
from os.path import getsize
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)

def L0_to_L1(wRocket, wFile, rocketFolderPath, justPrintFileNames,wflyer):


    if wRocket in [0,1,4,5]:

        # --- ACES II Flight/Integration Data ---
        rocketAttrs,b,c = ACES_mission_dicts()
        globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
        globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L1'
        L1ModelData = L1_TRICE_Quick(wflyer)


    elif wRocket in [2, 3]:

        # --- TRICE II ---
        globalAttrsMod = {}
        rocketAttrs,b,c = TRICE_mission_dicts()
        L1ModelData = L1_TRICE_Quick(wflyer)


    # Set the paths for the file names
    L0Files = glob(f'{rocketFolderPath}L0\{fliers[wflyer]}\*.cdf')
    L1Files = glob(f'{rocketFolderPath}L1\{fliers[wflyer]}\*.cdf')

    L0_names = [ifile.replace(f'{rocketFolderPath}L0\{fliers[wflyer]}\\', '') for ifile in L0Files]
    L1_names = [ofile.replace(f'{rocketFolderPath}L1\{fliers[wflyer]}\\', '') for ofile in L1Files]

    L0_names_searchable = [ifile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace('l0_', '').replace('_v00', '') for ifile in L0_names]
    L1_names_searchable = [ofile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace('l1_', '').replace('_v00','').replace('__', '_') for ofile in L1_names]

    dataFile_name = L0Files[wFile].replace(f'{rocketFolderPath}L0\{fliers[wflyer]}\\', '')
    fileoutName = dataFile_name.replace('l0', 'l1')

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            wInstr = [index, instr]


    if justPrintFileNames:
        for i, file in enumerate(L0Files):
            anws = ["yes" if L0_names_searchable[i] in L1_names_searchable else "no"]
            print('[{:.0f}] {:70s}{:5.1f} MB   Made L1: {:3s} '.format(i, L0_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to L1 data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(L0Files[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg('Loading data from L0Files')
        data_dict = {}
        with pycdf.CDF(L0Files[wFile]) as L0DataFile:
            for key, val in L0DataFile.items():
                data_dict = {**data_dict, **{key: [L0DataFile[key][...], {key:val for key, val in L0DataFile[key].attrs.items()  }  ]  }  }
        Done(start_time)


        # --- --- --- --- --- --- --- --- ---
        # --- Calculate Instrument Data ---
        # --- --- --- --- --- --- --- --- ---

        prgMsg('\nCreating L1 instrument data')

        # 0 -> EEPAA
        # 1 -> LEESA
        # 2 -> IEPAA
        # 3 -> LP
        if wInstr[0] in [0, 1, 2]:
            data_dict = {**data_dict, **{'Energy':           [[0], {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808, 'FORMAT': 'E12.2', 'UNITS': 'eV', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'log', 'LABLAXIS': 'Energy'}]}}
            data_dict = {**data_dict, **{'geometric_factor': [[0], {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808, 'FORMAT': 'E12.2', 'UNITS': 'cm!U2!N str ev ev!U-1!N', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]}}
            data_dict = {**data_dict, **{'Epoch_esa':        [[],  {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808, 'FORMAT': 'E12.2', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time', 'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear'}]}}
            data_dict['Pitch_Angle'] = [[0], {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808, 'FORMAT': 'E12.2', 'UNITS': 'deg', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear', 'LABLAXIS': 'Pitch_Angle'}]
            data_dict['28V_Monitor'][1]['DEPEND_0'] = 'Epoch_monitors'
            data_dict['Boom_Monitor'][1]['DEPEND_0'] = 'Epoch_monitors'

            sectorCounts = data_dict['Sector_Counts'][0]

            if energy_adjustment == 0:
                energy_adjust = -1*len(rocketAttrs.Instr_Energy[0])
            else:
                energy_adjust = energy_adjustment



            if wInstr[0] in [0]: # EEPAA
                data_dict['Pitch_Angle'][0] =np.array(rocketAttrs.Instr_sector_to_pitch[0])
                pitches = data_dict['Pitch_Angle'][0]
                data_dict['Energy'][0] = np.array(rocketAttrs.Instr_Energy[0][:(-1 * energy_adjust)])
                data_dict['geometric_factor'][0] = np.array(rocketAttrs.geometric_factor[0])

            elif wInstr[0] in [1]: # LEESA
                data_dict['Pitch_Angle'][0] = np.array(rocketAttrs.Instr_sector_to_pitch[1])
                pitches = data_dict['Pitch_Angle'][0]
                data_dict['Energy'][0] = np.array(rocketAttrs.Instr_Energy[1][:(-1 * energy_adjust)])
                data_dict['geometric_factor'][0] = np.array(rocketAttrs.geometric_factor[1])

            elif wInstr[0] in [2]: # IEPAA
                data_dict['Pitch_Angle'][0] = np.array(rocketAttrs.Instr_sector_to_pitch[2])
                pitches = data_dict['Pitch_Angle'][0]
                data_dict['Energy'][0] = np.array(rocketAttrs.Instr_Energy[2][:(-1 * energy_adjust)])
                data_dict['geometric_factor'][0] = np.array(rocketAttrs.geometric_factor[2])

            # --- PROCESS ESA DATA ---
            if energy_adjustment == 0:
                energy_adjust = 0
            else:
                energy_adjust = energy_adjustment

            sweepLength = len(data_dict['Energy'][0]) + (energy_adjust) + 1 # must be a total of 50 values in a sweep

            # --- find index locations of the start and end point of the sweeps ---
            # try:
            start_of_sweeps = np.where(data_dict['sweep_step'][0] == 49)[0] # before its turned on, sweep_steps are just zeros
            sweeps_start = start_of_sweeps[0] + 1 # index for the start of the first sweep, beginning on a retrace
            sweeps_end = start_of_sweeps[-1] + 1
            sectorCounts_reduced = sectorCounts[sweeps_start:sweeps_end]
            epoch_reduced = data_dict['Epoch'][0][sweeps_start:sweeps_end]
            no_of_sweeps = int(len(sectorCounts_reduced) / (sweepLength))

            # --- PLACE THE OUTPUT DATA ---
            counts = np.zeros(shape=(no_of_sweeps, len(pitches), sweepLength - (energy_adjust + 1)  ))

            for i in tqdm(range(no_of_sweeps)):

                # downsample the ESA's epoch so it starts at the beginning of each sweep
                data_dict['Epoch_esa'][0].append(epoch_reduced[i*sweepLength])

                # take subset of sector counts and epoch data with all pitch angles
                subSet = sectorCounts_reduced[i * sweepLength: (1 + i) * sweepLength]

                if energy_adjustment == 0:
                    remthese = [1]
                else:
                    remthese = [1,2]

                for j in (range(len(pitches))): # place the sector count data
                    sweep = [int(subSet[k][j]) for k in (range(sweepLength)) if k not in remthese]
                    sweep.insert(len(sweep)-1, sweep.pop(0)) # move the index 0 value to where the highest
                    counts[i][j] = sweep[::-1] # Invert the order to match the energy to the sweep no.

            data_dict['Epoch_esa'][0] = np.array(data_dict['Epoch_esa'][0])

            data_dict = {**data_dict, **{rocketAttrs.InstrNames_LC[wInstr[0]]:
                                             [counts, {'LABLAXIS': rocketAttrs.InstrNames_LC[wInstr[0]],
                                                       'DEPEND_0': 'Epoch_esa', 'DEPEND_1': 'Pitch_Angle',
                                                       'DEPEND_2': 'Energy',
                                                       'FILLVAL': 65535, 'FORMAT': 'E12.2',
                                                       'UNITS': 'counts', 'VALIDMIN': 0,
                                                       'VALIDMAX': rocketAttrs.esaMaxCounts, 'VAR_TYPE': 'data',
                                                       'SCALETYP': 'linear'}]}}

            del data_dict['Sector_Counts'], data_dict['sfid'], data_dict['sweep_step'],data_dict['minor_frame_counter'],data_dict['major_frame_counter']


            if truncateData:

                # --- --- --- --- --- ---
                # --- REDUCE DATASET ---
                # --- --- --- --- --- ---
                # Nothing interesting happens before 17:20 in the data, find this point and only keep data after it

                targetDate, targetTime  = dt.datetime(2022, 11, 20), dt.time(17, 20, 0, 0)
                targetDateTime_TT2000 = pycdf.lib.datetime_to_tt2000(targetDate.combine(targetDate, targetTime))
                Epoch_targetIndex = np.array([np.abs(pycdf.lib.datetime_to_tt2000(time) - targetDateTime_TT2000) for time in data_dict['Epoch'][0]]).argmin()
                Epoch_esa_targetIndex = np.array([np.abs(pycdf.lib.datetime_to_tt2000(time) - targetDateTime_TT2000) for time in data_dict['Epoch_esa'][0]]).argmin()
                Epoch_monitors_targetIndex = np.array([np.abs(pycdf.lib.datetime_to_tt2000(time) - targetDateTime_TT2000) for time in data_dict['Epoch_monitors'][0]]).argmin()

                noReduction = ['geometric_factor','Energy','Pitch_Angle','Sector_Number']
                needsReduction_esa = [rocketAttrs.InstrNames_LC[wInstr[0]], 'Epoch_esa']
                needsReduction_monitors = ['28V_Monitor','Boom_Monitor','Epoch_monitors']

                for key, val in data_dict.items():
                    if key in needsReduction_esa: # counts data
                        data_dict[key][0] = data_dict[key][0][Epoch_esa_targetIndex:]
                    elif key in needsReduction_monitors: # Monitors data
                        data_dict[key][0] = data_dict[key][0][Epoch_monitors_targetIndex:]
                    elif key not in noReduction: # Everything else
                        data_dict[key][0] = data_dict[key][0][Epoch_targetIndex:]

            dataFailed = False

            # except Exception as e:
            #     print(color.RED + f'Data was not processed. Reason: {e}' + color.END)
            #     dataFailed = True

            # --- --- --- --- --- --- --- ---
            # --- WRITE OUT THE ESA DATA ---
            # --- --- --- --- --- --- --- ---
            if not dataFailed:
                prgMsg('Creating output ESA file')

                outputPath = f'{rocketFolderPath}L1\{fliers[wflyer]}\\{fileoutName}'

                # --- delete output file if it already exists ---
                if os.path.exists(outputPath):
                    os.remove(outputPath)

                # --- open the output file ---
                with pycdf.CDF(outputPath, '') as L1File:
                    L1File.readonly(False)

                    # --- write out global attributes ---
                    inputGlobDic = L1ModelData.cdfFile.globalattsget()
                    for key, val in inputGlobDic.items():
                        if key == 'Descriptor':
                            globalAttrsMod[key] = rocketAttrs.InstrNames_Full[wInstr[0]]
                        if key in globalAttrsMod:
                            L1File.attrs[key] = globalAttrsMod[key]
                        else:
                            L1File.attrs[key] = val

                    # --- WRITE OUT DATA ---
                    for varKey, varVal in data_dict.items():
                        if varKey in ['Epoch', 'Epoch_monitors', 'Epoch_esa']:  # epoch data
                            L1File.new(varKey, data=varVal[0], type=33)
                        elif varKey == wInstr[1]:  # instrument data
                            L1File.new(varKey, data=varVal[0])
                        else:  # support data
                            L1File.new(varKey, data=varVal[0])

                        # --- Write out the attributes and variable info ---
                        for attrKey, attrVal in data_dict[varKey][1].items():
                            if attrKey == 'VALIDMIN':
                                L1File[varKey].attrs[attrKey] = varVal[0].min()
                            elif attrKey == 'VALIDMAX':
                                L1File[varKey].attrs[attrKey] = varVal[0].max()
                            elif attrVal != None:
                                L1File[varKey].attrs[attrKey] = attrVal

            Done(start_time)

        elif wInstr[0] in [3]:

            # --- --- --- --- --- --- --- --- --- ---
            # --- SEPARATE AND COLLECT LP DATA ---
            # --- --- --- --- --- --- --- --- --- ---
            try:
                # --- Get L0 data ---
                L0epoch = [pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][i]) for i in range(len(data_dict['Epoch'][0]))]
                L0ChannelCounts = data_dict['Channel_Counts'][0]
                L0sfid = data_dict['sfid'][0]
                samplePeriod = rocketAttrs.LPSamplePeriod
                minorFrameTime = rocketAttrs.MinorFrameTime # each LP value is really sampled 1 minorframe before its reported value in the PCM matrix

                del data_dict['Epoch'], data_dict['Channel_Counts'],data_dict['Channel_Number'], data_dict['minor_frame_counter'],data_dict['major_frame_counter'], data_dict['sfid']
                collect()

                # --- --- --- --- --- ---
                # --- REDUCE DATASET ---
                # --- --- --- --- --- ---

                if truncateData:
                    # Nothing interesting happens before 17:20 in the data, find this point and only keep data after it
                    targetDate, targetTime = dt.datetime(2022, 11, 20), dt.time(17, 20, 0, 0)
                    targetDateTime_TT2000 = pycdf.lib.datetime_to_tt2000(targetDate.combine(targetDate, targetTime))
                    Epoch_targetIndex = np.array([np.abs(time - targetDateTime_TT2000) for time in L0epoch]).argmin()

                    # --- Backtrack to where sfid = 0 ---
                    dataset_targetIndex = 0
                    for i in range(0,len(L0sfid)):
                        if L0sfid[Epoch_targetIndex-1*i] == 0:
                            dataset_targetIndex = Epoch_targetIndex-1*i
                            break
                else:
                    # if we don't truncate, just find where sfid = 0 and cut everything else
                    dataset_targetIndex = np.where(L0sfid == 0)[0][0]


                # --- reduce the dataset ---
                L0epoch = L0epoch[dataset_targetIndex:]
                L0ChannelCounts = L0ChannelCounts[dataset_targetIndex:]
                L0sfid = L0sfid[dataset_targetIndex:]


                # create the 5 variables in the LP data_dict: deltaNdivN,ni,ne_swept,ni_swept,step
                for i in range(len(rocketAttrs.LP_Variables)):
                    data_dict = {**data_dict, **{f'Epoch_{rocketAttrs.LP_Variables[i]}': [[], {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808, 'FORMAT': 'E12.2', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time', 'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear'}]}}
                    data_dict = {**data_dict, **{f'{rocketAttrs.LP_Variables[i]}': [[], {'LABLAXIS': rocketAttrs.LP_Variables[i], 'DEPEND_0': f"Epoch_{rocketAttrs.LP_Variables[i]}", 'DEPEND_1':None, 'DEPEND_2':None, 'FILLVAL': -1, 'FORMAT': 'E12.2', 'UNITS': 'Digital', 'VALIDMIN': 0, 'VALIDMAX': rocketAttrs.esaMaxCounts, 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}

                no_of_fillvals = []

                for j in tqdm(range(len(L0epoch))):

                    LP_var_counter = L0sfid[j]%4 + 1 # sets which LP_var I'm on. (step, ne_swept, ni_swept, ni)

                    # deltaNdivN takes the first 8 values
                    data_dict[rocketAttrs.LP_Variables[0]][0].append(L0ChannelCounts[j][0:8])

                    if L0epoch[j] != -9223372036854775808:

                        # assign epoch values for each point with a sample rate of 35.25kHz
                        data_dict[f'Epoch_{rocketAttrs.LP_Variables[0]}'][0].append([(L0epoch[j] - minorFrameTime) + l * samplePeriod for l in range(8)])

                        if L0ChannelCounts[j][8] > 4095: # set all values above 4095 to the fillval
                            data_dict[rocketAttrs.LP_Variables[LP_var_counter]][0].append(-1)
                        else:
                            data_dict[rocketAttrs.LP_Variables[LP_var_counter]][0].append(L0ChannelCounts[j][8]) # LP_vars

                        data_dict[f'Epoch_{rocketAttrs.LP_Variables[LP_var_counter]}'][0].append((L0epoch[j] - minorFrameTime) + 8 * samplePeriod)

                    else: # if epoch is a fillval
                        no_of_fillvals.append([rocketAttrs.LP_Variables[LP_var_counter],j])
                        data_dict[f'Epoch_{rocketAttrs.LP_Variables[0]}'][0].append([-9223372036854775808 for l in range(8)]) # deltaNdivN
                        data_dict[rocketAttrs.LP_Variables[LP_var_counter]][0].append(-1)  # LP_vars
                        data_dict[f'Epoch_{rocketAttrs.LP_Variables[LP_var_counter]}'][0].append(-9223372036854775808) # LP_vars


                # Converts  all data to numpy arrays
                for key, var in data_dict.items():
                    data_dict[key][0] = np.array(data_dict[key][0])

                data_dict['deltaNdivN'][0] = np.array(data_dict['deltaNdivN'][0]).flatten()
                data_dict['Epoch_deltaNdivN'][0] = np.array(data_dict['Epoch_deltaNdivN'][0]).flatten()

                dataFailed = False


            except Exception as e:
                dataFailed = True
                print(color.RED + f"{e}" + color.END)

            # --- --- --- --- --- ---
            # --- WRITE OUT DATA ---
            # --- --- --- --- --- ---

            if not dataFailed:

                prgMsg('\nCreating output LP file')

                # LP data is too large to be dumped all at once, must break it up
                for i in range(len(rocketAttrs.LP_Variables)):

                    # determine the outputPath and the variables which need to be output'd
                    outputPath = f'{rocketFolderPath}L1\{fliers[wflyer]}\\{fileoutName.replace("lp", f"lp_{rocketAttrs.LP_Variables[i]}")}'
                    noWriteVarst = ["deltaNdivN","step", "ne_swept", "ni_swept","ni"]
                    noWriteVarst.pop(i)
                    noWriteVars = [thing for thing in noWriteVarst]

                    for thing in noWriteVarst:
                        noWriteVars.append(f'Epoch_{thing}')

                    writeOutVars = [key for key, val in data_dict.items() if key not in noWriteVars]

                    # create another dictionary that contains only the needed variables
                    data_dict_writeOut = {}
                    data_dict_copy = deepcopy(data_dict)
                    for key, val in data_dict_copy.items():
                        if key in writeOutVars:
                            data_dict_writeOut = {**data_dict_writeOut, **{key:val}}

                    del data_dict_copy
                    gc.collect()

                    # Convert everything to numpy array
                    for key, val in data_dict_writeOut.items():
                        data_dict_writeOut[key][0] = np.array(data_dict_writeOut[key][0])

                    # --- --- --- --- --- ---
                    # --- WRITE OUT DATA ---
                    # --- --- --- --- --- ---

                    # --- delete output file if it already exists ---
                    if os.path.exists(outputPath):
                        os.remove(outputPath)

                    with pycdf.CDF(outputPath,'') as outputFile:
                        outputFile.readonly(False)

                        # --- write out global attributes ---
                        inputGlobDic = L1ModelData.cdfFile.globalattsget()
                        for key, val in inputGlobDic.items():
                            if key == 'Descriptor':
                                globalAttrsMod[key] = rocketAttrs.InstrNames_Full[wInstr[0]]

                            if key in globalAttrsMod:
                                outputFile.attrs[key] = globalAttrsMod[key]
                            else:
                                outputFile.attrs[key] = val


                        for varKey, varVal in data_dict_writeOut.items():
                            if 'Epoch' in varKey:
                                outputFile.new(varKey, data=varVal[0], type=33)
                            elif 'ni_swept' in varKey or 'ne_swept' in varKey:
                                outputFile.new(varKey, data=varVal[0], type=pycdf.const.CDF_INT4)
                            else:
                                outputFile.new(varKey, data=varVal[0])

                            # --- Write out the attributes and variable info ---
                            for attrKey, attrVal in data_dict_writeOut[varKey][1].items():
                                if attrKey == 'VALIDMIN':
                                    outputFile[varKey].attrs[attrKey] = varVal[0].min()
                                elif attrKey == 'VALIDMAX':
                                    outputFile[varKey].attrs[attrKey] = varVal[0].max()
                                elif attrVal != None:
                                    outputFile[varKey].attrs[attrKey] = attrVal

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

if len(glob(f'{rocketFolderPath}L0\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        L0_to_L1(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}L0\{fliers[wflyer]}\*.cdf')))):
            L0_to_L1(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            L0_to_L1(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)