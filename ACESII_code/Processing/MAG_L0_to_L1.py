# --- MAG_L0_to_L1.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Magnetometer data to Science Units AND adjust Epoch-range with time-offset

# Truncates all data to only that past 17:20:00. Look at tmCDF_to_L0 to see when Mag data Epoch
# is sampled (should be sfid == 1 & 21)

# Because we forgot to take into account that the coil system at Iowa and at Wallops
# are left-handed to minimize gradients.  We simply flipped the magnetometer axis Y-Axis the
# wrong way rather than flip the reference data axis. So we add a minus sign to the mag Y-axis in
# our code.

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
wFiles = [0]

modifier = ''
inputPath_modifier = 'L0' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'L1' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def MAG_L0_to_L1(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L1'
    outputModelData = L1_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*RingCore*')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]

    input_names_searchable = [ifile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('.cdf','') for ifile in input_names]

    fileoutName = f'ACESII_{rocketID}_{outputPath_modifier.lower()}_{input_names_searchable[wFile]}.cdf'
    fileoutName = fileoutName.replace('_magFrm', '_rktFrm')

    # determine which instrument the file corresponds to:
    mags = ['RingCore', 'Tesseract']
    for index, instr in enumerate(mags):
        if instr.lower() in inputFiles[wFile].lower():
            wInstr = [index, instr, mags[index]]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:35s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {input_names_searchable[wFile]}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict_mag = loadDictFromFile(inputFiles[wFile], {})
        data_dict_mag['Epoch'][0] = [pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag['Epoch'][0]]
        Done(start_time)

        #######################
        # --- TRUNCATE DATA ---
        #######################
        prgMsg('Truncating Data')

        # --- do a quality check on the data to ensure epoch is > 0 or at least fillval ---
        # NOTE: There's some weirdness where the data doesn't recognize a negative epoch as a fillval. Below code
        # removes the weirdness
        badIndicies = []
        for i, tme in enumerate(data_dict_mag['Epoch'][0]):
            if tme < 0:
                badIndicies.append(i)

        for key, val in data_dict_mag.items():
            data_dict_mag[key][0] = np.delete(data_dict_mag[key][0], badIndicies)

        # --- remove data previous to 17:20:00 ---
        targetTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 000000))
        truncIndex = np.abs(np.array(data_dict_mag['Epoch'][0]) - targetTime).argmin()
        for key, val in data_dict_mag.items():
            data_dict_mag[key][0] = data_dict_mag[key][0][truncIndex:]
        Done(start_time)

        # #################################
        # --- CONVERT TO SCIENCE UNITS ---
        # #################################
        prgMsg('Converting to Science Units')

        B_engineering_units = np.array([[data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]] for i in range(len(data_dict_mag['Epoch'][0]))])
        def ringCoreADCtoScience(B):

            # Linear correction coefficents
            linearCorrection = rocketAttrs.ringCoreLinearCorrections[wflyer]
            B = [linearCorrection[0] * B[0], linearCorrection[1] * B[1], linearCorrection[2] * B[2]]

            # non-linear correction 5th, order coeffients. The terms descend in exponent order i.e. fithOrderCorrection[0] --> x^5 term, etc
            fithOrderCorrection = rocketAttrs.ringCore5thOrderCorrections[wflyer]

            for i in range(len(B)):
                B[i] = B[i] * fithOrderCorrection[i][4] + (B[i] ** 2) * fithOrderCorrection[i][3] + (B[i] ** 3) * \
                       fithOrderCorrection[i][2] + (B[i] ** 4) * fithOrderCorrection[i][1] + (B[i] ** 5) * \
                       fithOrderCorrection[i][0]
            return B

        # convert the data and assign it back to the data dict
        B_science_units = np.array([ringCoreADCtoScience(B) for B in B_engineering_units])
        data_dict_mag['Bx'][0], data_dict_mag['By'][0], data_dict_mag['Bz'][0] = B_science_units[:, 0], B_science_units[:, 1], B_science_units[:, 2]
        Done(start_time)

        # ############################
        # --- CORRECT FLIPPED AXIS ---
        # ############################
        # The magnetometer Y-axis was flipped during calibration. Here we add back in the minus sign
        data_dict_mag['By'][0] = -1*data_dict_mag['By'][0]

        # ###############################
        # --- ROTATE INTO ROCKETFRAME ---
        # ###############################
        # The magnetometer frame can be converted to rocket-body frame by a 90deg rotation about Z-axis
        from ACESII_code.class_var_func import Rz
        B_rktFrm = np.array([np.matmul(Rz(90), np.array([data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]])) for i in range(len(data_dict_mag['Epoch'][0]))])
        data_dict_mag['Bx'][0], data_dict_mag['By'][0], data_dict_mag['Bz'][0] = np.array(B_rktFrm[:, 0]), np.array(B_rktFrm[:, 1]), np.array(B_rktFrm[:, 2])

        # ##########################
        # --- APPLY CALIBRATIONS ---
        # ##########################
        # description: Apply the calibration matrix which is a vector cal with the model: orthogonality, offset, and Euler rotation.
        # It is a coupled matrix containing all cals found via robust weighted regression.
        # Columns 1:3 are multiplied by the spun data then column 4 is added (basically a line equation in 3-D)
        prgMsg('Applying Calibrations')
        calMatrix = rocketAttrs.ringCoreCalMatrix[wRocket-4]
        orthoEuler = np.array([calMatrix[i][0:3] for i in range(len(calMatrix))])
        offset = np.array([calMatrix[i][3] for i in range(len(calMatrix))])
        B_cald = np.array([np.matmul(orthoEuler, np.array([data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]])) + offset for i in range(len(data_dict_mag['Epoch'][0]))])
        data_dict_mag['Bx'][0], data_dict_mag['By'][0], data_dict_mag['Bz'][0] = np.array(B_cald[:, 0]), np.array(B_cald[:, 1]), np.array(B_cald[:, 2])
        Done(start_time)

        ###############################
        # --- FIX SENSITIVITY ISSUE ---
        ###############################

        # description: Bx needs to be multiplied by -1. This was discovered after sanity checking
        # the calibration parameters and seeing that the X sensitivity was negative (which is impossible)
        # Bx in rocket body was the original Y term in the magnetometer frame, so that is explained
        # by the flipped axis.

        prgMsg('Fixing Sensitivity')
        data_dict_mag['By'][0] = np.array([-1*data_dict_mag['By'][0][i] for i in range(len(data_dict_mag['By'][0]))])
        Done(start_time)

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('Creating output file')

            for key,val in data_dict_mag.items():
                data_dict_mag[key][0] = np.array(data_dict_mag[key][0])

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'
            outputCDFdata(outputPath, data_dict_mag, outputModelData, globalAttrsMod, wInstr[1])
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
        MAG_L0_to_L1(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            MAG_L0_to_L1(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            MAG_L0_to_L1(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)