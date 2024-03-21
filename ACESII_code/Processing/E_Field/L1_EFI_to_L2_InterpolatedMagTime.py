# --- L1_EFI_to_L2_InterpolatedMagTime.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: little script to interpolated the despun EFI data onto the mag frame.
# Justified since the signals I'm interested in are low frequency and this step speeds up
# my analysis a lot



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

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier_elec = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_mag = 'l2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'l2' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def L1_EFI_to_L2_InterpolatedMagTime(wRocket, wFile, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    globalAttrsMod = rocketAttrs.globalAttributes[wRocket-4]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wRocket-4)

    inputFiles_elec = glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wRocket-4]}{modifier}\*E_Field*')
    inputFiles_mag = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wRocket-4]}{modifier}\*RingCore_Despun_ENU*')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wRocket-4]}{modifier}\\', '') for ifile in inputFiles_elec]
    fileoutName = input_names[wFile].replace('l1','l2')


    if justPrintFileNames:
        for i, file in enumerate(inputFiles_elec):
            print('[{:.0f}] {:80s}'.format(i, input_names[i], round(getsize(file) / (10 ** 6))))
    else:

        # --- get the data from the B-Field file ---
        prgMsg(f'Loading data from {inputPath_modifier_mag} Files')
        data_dict_mag = loadDictFromFile(inputFiles_mag[wFile], {},reduceData=False,targetTimes=[])
        Done(start_time)

        # --- get the data from the E-Field file ---
        prgMsg(f'Loading data from {inputPath_modifier_elec} Files')
        data_dict_elec = loadDictFromFile(inputFiles_elec[0], {},reduceData=False,targetTimes=[])
        Done(start_time)

        ##########################################################
        # --- interpolate E-Field data onto magnetometer epoch ---
        ##########################################################
        prgMsg('Interpolating E-Field Data')
        data_dict_elecInterp = InterpolateDataDict(InputDataDict=data_dict_elec,
                                                       InputEpochArray=data_dict_elec['Epoch'][0],
                                                       wKeys=[],
                                                       targetEpochArray=data_dict_mag['Epoch'][0])

        Done(start_time)


        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('Creating output file')

            # create the output data_dict
            data_dict_output = deepcopy(data_dict_elecInterp)

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wRocket-4]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict_output, outputModelData, globalAttrsMod, 'EFI')

            Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder

if len(glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        L1_EFI_to_L2_downsampled(wRocket, 0, rocketFolderPath, justPrintFileNames,wRocket-4)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wRocket-4]}\*E_Field*')))):
            L1_EFI_to_L2_downsampled(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wRocket-4)
    else:
        for filesNo in wFiles:
            L1_EFI_to_L2_downsampled(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wRocket-4)