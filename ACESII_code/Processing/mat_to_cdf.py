# --- template.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:



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
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0,1]

modifier = ''
inputPath_modifier = 'mag_formatted' # e.g. 'L1' or 'L1'. It's the name of the broader input folder inside data\ACESII
outputPath_modifier = 'mag' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder inside data\ACESII\ACESII_matlab


rotateIntoRingCoreFrame = False
flipYAxis = True


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os
import scipy
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg, L0_ACES_Quick,L0_TRICE_Quick
from glob import glob
from os.path import getsize

setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)
from mat73 import loadmat
from copy import deepcopy


###########################
# --- SPECIAL VARIABLES ---
###########################

# Program default: recognize Epoch variable, turn it into "support data". Turn other data into "data"
# Below you can add additional specifications to known variables
#
# special_mods = {'Epoch': {'UNITS':'ns','LABLAXIS':'Epoch'},
#                 'Bx': {'UNITS': 'nT', 'LABLAXIS': 'Bx'},
#                 'By': {'UNITS': 'nT', 'LABLAXIS': 'By'},
#                 'Bz': {'UNITS': 'nT', 'LABLAXIS': 'Bz'},
#                 'Bx_Model': {'UNITS': 'nT', 'LABLAXIS': 'Bx'},
#                 'By_Model': {'UNITS': 'nT', 'LABLAXIS': 'By'},
#                 'Bz_Model': {'UNITS': 'nT', 'LABLAXIS': 'Bz'}
#                 }
special_mods = {'Epoch': {'UNITS':'ns','LABLAXIS':'Epoch'},
                'Bx': {'UNITS': 'nT', 'LABLAXIS': 'Bx'},
                'By': {'UNITS': 'nT', 'LABLAXIS': 'By'},
                'Bz': {'UNITS': 'nT', 'LABLAXIS': 'Bz'},
                }


def main(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    inputrocketFolderPath = rocketFolderPath + 'ACESII_matlab'
    outputrocketFolderPath = rocketFolderPath

    # --- ACESII ---
    if wRocket in [0, 1, 4, 5]:
        # --- ACES II Flight/Integration Data ---
        rocketAttrs, b, c = ACES_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
        outputModelData = L0_ACES_Quick(wflyer)

    # --- TRICE II ---
    elif wRocket in [2, 3]:
        globalAttrsMod = {}
        rocketAttrs, b, c = TRICE_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        outputModelData = L0_TRICE_Quick(wflyer)

    inputFiles = glob(f'{inputrocketFolderPath}\{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.mat')
    outputFiles = glob(f'{outputrocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    input_names = [ifile.replace(f'{inputrocketFolderPath}\{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    output_names = [ofile.replace(f'{outputrocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace('ACES_', '').replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{inputrocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')

    fileoutName = dataFile_name.replace(f'{inputrocketFolderPath}\\{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', "").replace('.mat','.cdf')

    for index, instr in enumerate(['RingCore','Tesseract']):
        if instr.lower() in dataFile_name.lower():
            wInstr = [index, instr]


    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made cdf: {:3s} '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        ##########################
        # --- DEFINE VARIABLES ---
        ##########################
        prgMsg('Converting variables to .cdf format')

        exampleVar = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                        'FORMAT': 'I5', 'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                        'SCALETYP': 'linear', 'LABLAXIS': None}
        exampleEpoch = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                  'FORMAT': 'I5', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data',
                  'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time',
                  'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear', 'LABLAXIS': 'Epoch'}


        # DETERMINE THE FILETYPE
        loadThisFile = inputFiles[wFile]

        try:
            # Load in data to data_dict
            mat = loadmat(loadThisFile)

        except Exception as e:
            print('\n')
            print(e,'\n')
            prgMsg('Loading with scipy.io instead')
            mat = scipy.io.loadmat(loadThisFile)

        # --- convert data into dictionary ---
        data_dict = {}
        iterAble = {}
        counter = 0
        for key, val in special_mods.items():
            iterAble = {**iterAble,**{key: np.array(mat['output'][0][0][counter]).flatten()}}
            counter += 1

        # Reformat files into CDF format.
        for key, val in iterAble.items():
            # put .mat data into data_dict
            if 'Epoch' in key:
                data_dict = {**data_dict, **{key: [val, deepcopy(exampleEpoch)]}}
            else:
                data_dict = {**data_dict, **{key: [val, deepcopy(exampleVar) ]}}

            # Insert the special modifications
            if key in special_mods:

                for attrNam,attrVal in special_mods[key].items():

                    if attrNam in data_dict[key][1]: # if attribute already exists, overwrite it
                        data_dict[key][1][attrNam] = attrVal

                    else: # if attribute doesn't exist yet
                        data_dict[key][1] = {**data_dict[key][1],**{attrNam:attrVal}}


        # --- Handle the Time variable ---
        # If "Time" variable is actually "Time since launch" do the conversion from SECONDS
        if data_dict['Epoch'][0][0] < 100000000000000000:
            data_dict['Epoch'][0] = np.array([data_dict['Epoch'][0][i]*(10**(9)) + rocketAttrs.Launch_Times[wRocket - 4] for i in range(len(data_dict['Epoch'][0]))])

        # Calculate the Bmag value
        BmagAttrs = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                        'FORMAT': 'I5', 'UNITS': 'nT', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                        'SCALETYP': 'linear', 'LABLAXIS': 'Bmag'}

        # Real data
        Bmag = np.array(
            [
              np.linalg.norm([data_dict['Bx'][0][i],data_dict['By'][0][i],data_dict['Bz'][0][i]])
                for i in range(len(data_dict['Bx'][0]))
            ])

        data_dict = {**data_dict,**{'Bmag':[Bmag,BmagAttrs]}}

        # Model data
        # Bmag_Model = np.array(
        #     [
        #         np.linalg.norm([data_dict['Bx_Model'][0][i], data_dict['By_Model'][0][i], data_dict['Bz_Model'][0][i]])
        #         for i in range(len(data_dict['Bx_Model'][0]))
        # #     ])
        #
        # BmagAttrs['LABLAXIS'] = 'Bmag_Model'
        # data_dict = {**data_dict, **{'Bmag_Model': [Bmag_Model, BmagAttrs]}}

        Done(start_time)

        if flipYAxis:
            data_dict['By'][0] = np.array([-1*data_dict['By'][0][i] for i in range(len(data_dict['By'][0]))])

        # --- --- --- --- --- ---
        # --- ROTATE THE DATA ---
        # --- --- --- --- --- ---
        if rotateIntoRingCoreFrame:

            # Robert put the mag data into the payload frame to make his life easier. I need to convert it back to it's original frame
            from ACESII_code.class_var_func import Rz

            Bvec = np.array([
                    [data_dict['Bx'][0][i],
                     data_dict['By'][0][i],
                     data_dict['Bz'][0][i]]
                 for i in range(len(data_dict['Bmag'][0]))]
            )

            # rotate the components back. The -1 is multiplied to account for the point of view perspective.
            # I like to think in terms of rotating the axes themselves, with is the OPPOSITE angle change as the vector changing.
            Bcomps_rotated =np.array([
                np.matmul(Rz(-1*-90),Bvec[i]) for i in range(len(data_dict['Bmag'][0]))
                ]
            )

            # rename Bx,By,Bz to what they really are: B in the payload coordinates
            data_dict['Bx_payload'] = data_dict.pop('Bx')
            data_dict['By_payload'] = data_dict.pop('By')
            data_dict['Bz_payload'] = data_dict.pop('Bz')

            # store new "RingCore" B. The coordinates are now aligned to the RingCore instrument frame
            data_dict = {**data_dict, **{'Bx': [np.array([Bcomps_rotated[i][0] for i in range(len(Bcomps_rotated))]), data_dict['Bx_payload'][1]]}}
            data_dict = {**data_dict, **{'By': [np.array([Bcomps_rotated[i][1] for i in range(len(Bcomps_rotated))]), data_dict['By_payload'][1]]}}
            data_dict = {**data_dict, **{'Bz': [np.array([Bcomps_rotated[i][2] for i in range(len(Bcomps_rotated))]), data_dict['Bz_payload'][1]]}}


        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        prgMsg('Creating output file')

        outputPath = f'{outputrocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

        # --- --- --- --- --- ---
        # --- WRITE OUT DATA ---
        # --- --- --- --- --- ---

        # --- delete output file if it already exists ---
        if os.path.exists(outputPath):
            os.remove(outputPath)

        with pycdf.CDF(outputPath, '') as outputFile:
            outputFile.readonly(False)

            globalAttrsMod['Descriptor'] = wInstr[1]

            # --- write out global attributes ---
            inputGlobDic = outputModelData.cdfFile.globalattsget()
            for key, val in inputGlobDic.items():
                if key in globalAttrsMod:
                    outputFile.attrs[key] = globalAttrsMod[key]
                else:
                    outputFile.attrs[key] = val

            for varKey, varVal in data_dict.items():
                if 'epoch' in varKey.lower():
                    outputFile.new(varKey, data=varVal[0], type=33)
                else:
                    outputFile.new(varKey, data=varVal[0])

                # --- Write out the attributes and variable info ---
                for attrKey, attrVal in data_dict[varKey][1].items():
                    if attrKey == 'VALIDMIN':
                        outputFile[varKey].attrs[attrKey] = varVal[0].min()
                    elif attrKey == 'VALIDMAX':
                        outputFile[varKey].attrs[attrKey] = varVal[0].max()
                    elif attrKey == 'LABLAXIS':
                        outputFile[varKey].attrs[attrKey] = varKey
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

if len(glob(f'{rocketFolderPath}\ACESII_matlab\{inputPath_modifier}\{fliers[wflyer]}\*.mat')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        main(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.mat')))):
            main(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            main(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)