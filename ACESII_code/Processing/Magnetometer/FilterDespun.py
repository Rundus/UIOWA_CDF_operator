# --- FilterDespun.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: script to convert .cdfs to .mat so they can be filtered in MATLAB
# OR convert the same files back to .cdfs



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

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

modifier = ''
inputPath_modifier = 'l2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder inside data\ACESII
outputPath_modifier_E = 'science\DespunE_filtered' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder inside data\ACESII\ACESII_matlab
outputPath_modifier_B = 'science\DespunB_filtered' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder inside data\ACESII\ACESII_matlab



wBfile = 0
wEfile =0


SECTION_convertToMAT = True
SECTION_convertFilteredMATToCDF = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.io import savemat


def FilterDespun(wRocket, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L0_ACES_Quick(wflyer)


    inputFilesB = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*RingCore_DeSpun_ENU*')
    inputFilesE = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*RingCore_DeSpun_ENU*')


    if justPrintFileNames:
        print('--- B ---')
        for i, file in enumerate(inputFilesB):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFilesB[i], round(getsize(file) / (10 ** 6), 1)))

        print('\n--- E ---')
        for i, file in enumerate(inputFilesB):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFilesE[i], round(getsize(file) / (10 ** 6), 1)))
    else:


        if SECTION_convertToMAT:

            if wRocket == 4:


                data_dict_mag = loadDictFromFile(inputFilesB[wBfile],{},reduceData=False,targetTimes=[])
                data_dict_mag['Epoch'][0] = dateTimetoTT2000(data_dict_mag['Epoch'][0],False)

                #####################################
                # --- prepare dictionary for .mat ---
                #####################################
                # remove the attributes attached to the dictonary
                mdic = {}
                for key, val in data_dict_mag.items():
                    mdic = {**mdic, **{key: val[0]}}

                # --- --- --- --- --- --- ---
                # --- WRITE OUT THE DATA ---
                # --- --- --- --- --- --- ---
                prgMsg('Creating output file')
                fileoutNameB = 'ACESII_36359_l2_RingCore_DeSpun_ENU.mat'
                outputPath = f'{rocketFolderPath}\ACESII_matlab\{outputPath_modifier_B}\{fliers[wflyer]}\\{fileoutNameB}'
                # --- delete output file if it already exists ---
                if os.path.exists(outputPath):
                    os.remove(outputPath)

                # output the data
                savemat(outputPath, mdict=mdic)

                Done(start_time)

            elif wRocket == 5:

                data_dict_mag = loadDictFromFile(inputFilesB[wBfile], {}, reduceData=False, targetTimes=[])
                data_dict_mag['Epoch'][0] = dateTimetoTT2000(data_dict_mag['Epoch'][0], False)
                data_dict_elec = loadDictFromFile(inputFilesE[wEfile], {}, reduceData=False, targetTimes=[])
                data_dict_elec['Epoch'][0] = dateTimetoTT2000(data_dict_elec['Epoch'][0], False)

                #####################################
                # --- prepare dictionary for .mat ---
                #####################################
                # remove the attributes attached to the dictonary
                mdic_mag = {}
                for key, val in data_dict_mag.items():
                    mdic_mag = {**mdic_mag, **{key: val[0]}}

                mdic_elec = {}
                for key, val in data_dict_elec.items():
                    mdic_elec = {**mdic_elec, **{key: val[0]}}

                # --- --- --- --- --- --- ---
                # --- WRITE OUT THE DATA ---
                # --- --- --- --- --- --- ---
                prgMsg('Creating output file for MAG')

                fileoutNameB = 'ACESII_36364_l2_RingCore_DeSpun_ENU.mat'
                outputPath = f'{rocketFolderPath}\ACESII_matlab\{outputPath_modifier_B}\{fliers[wflyer]}\\{fileoutNameB}'

                # --- delete output file if it already exists ---
                if os.path.exists(outputPath):
                    os.remove(outputPath)

                # output the data
                savemat(outputPath, mdict=mdic_mag)

                Done(start_time)

                prgMsg('Creating output file for E-Field')
                fileoutNameE = 'ACESII_36364_l2_flight_E_Field_ENU.mat'
                outputPath = f'{rocketFolderPath}\ACESII_matlab\{outputPath_modifier_E}\{fliers[wflyer]}\\{fileoutNameE}'
                # --- delete output file if it already exists ---
                if os.path.exists(outputPath):
                    os.remove(outputPath)

                # output the data
                savemat(outputPath, mdict=mdic_elec)

                Done(start_time)



        elif SECTION_convertFilteredMATToCDF:
            print('potoat')


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
    FilterDespun(wRocket, rocketFolderPath, justPrintFileNames,wflyer)

