# --- L0_to_L1.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert electrostatic analyzer data from diffNFlux to Distribution Function



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

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

IsElectron = True
wIon = 0


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, L2_ACES_Quick,L2_TRICE_Quick, prgMsg, cm_to_m,q0,IonMasses,m_e
from glob import glob
from os.path import getsize
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


def Distribution_Function(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    outputFolderPath = rocketFolderPath + r'science\\DistFunc\\'

    if wRocket in [0,1,4,5]:

        # --- ACES II Flight/Integration Data ---
        rocketAttrs,b,c = ACES_mission_dicts()
        globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
        globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'Distribution_Function'
        L2ModelData = L2_ACES_Quick(wflyer)

    # --- TRICE II ---
    elif wRocket in [2,3]:
        globalAttrsMod = {}
        rocketAttrs,b,c = TRICE_mission_dicts()
        L2ModelData = L2_TRICE_Quick(wflyer)


    # Set the paths for the file names
    L2Files = glob(f'{rocketFolderPath}L2\{fliers[wflyer]}\*.cdf')
    DistFuncFiles = glob(f'{outputFolderPath}DistFunc\{fliers[wflyer]}\*.cdf')

    L2_names = [file.replace(f'{rocketFolderPath}L2\{fliers[wflyer]}\\', '') for file in L2Files]
    DistFunc_names = [file.replace(f'{outputFolderPath}DistFunc\{fliers[wflyer]}\\', '') for file in DistFuncFiles]

    L2_names_searchable = [file.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace('l2_', '').replace('_v00', '').replace('__', '_') for file in L2_names]
    DistFunc_names_searchable = [file.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace('distFunc_', '').replace('_v00', '') for file in DistFunc_names]

    dataFile_name = L2_names[wFile].replace(f'{rocketFolderPath}\{fliers[wflyer]}\\', '')
    fileoutName = dataFile_name.replace('l2', 'distFunc')


    if justPrintFileNames:
            for i, file in enumerate(L2Files):
                anws = ["yes" if L2_names_searchable[i].replace('.cdf', "") in DistFunc_names_searchable else "no"]
                print('[{:.0f}] {:70s}{:5.1f} MB   Made L2: {:3s} '.format(i, L2_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Calculating Distribution Function for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(L2Files[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg('Loading data from L2Files')
        data_dict = {}
        with pycdf.CDF(L2Files[wFile]) as L2DataFile:
            for key, val in L2DataFile.items():
                data_dict = {**data_dict, **{key : [L2DataFile[key][...] , {key:val for key,val in L2DataFile[key].attrs.items()  }  ]  }  }

        Done(start_time)

        # --- --- --- --- --- --- --- --- ---
        # --- Calculate Instrument Data ---
        # --- --- --- --- --- --- --- --- ---

        prgMsg('Calculating the Distribution Function')

        # --- CALCULATE DISTRIBUTION FUNCTION ---
        diffNFlux = data_dict['Differential_Number_Flux'][0]
        pitchAngle = data_dict['Pitch_Angle'][0]
        Energies = data_dict['Energy'][0]

        # define empty numpy array
        sizes = [len(diffNFlux),len(diffNFlux[0]), len(diffNFlux[0][0])]
        ranges = [range(sizes[0]), range(sizes[1]), range(sizes[2])]
        distFunc = np.zeros(shape=(sizes[0], sizes[1], sizes[2]))

        if IsElectron:
            m = m_e
        else:
            m = IonMasses[wIon]

        # --- Calculate DistFunc in SI units ---
        for tme, ptch, engy in tqdm(itertools.product(*ranges)):
            if diffNFlux[tme][ptch][engy] == -1e30:
                distFunc[tme][ptch][engy] = -1e30
            else:
                distFunc[tme][ptch][engy] = (cm_to_m*cm_to_m/(q0*q0))*(((m**2)*diffNFlux[tme][ptch][engy]) / (2 * Energies[engy]))

        distFunc = np.array(distFunc)

        del data_dict['Differential_Number_Flux'],data_dict['Differential_Energy_Flux']

        Done(start_time)

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        prgMsg('Creating output file')

        outputPath = f'{outputFolderPath}{fliers[wflyer]}\\{fileoutName}'

        data_dict = {**data_dict, **{'Distribution_Function':
                                         [distFunc, {'LABLAXIS': 'Distribution_Function',
                                                   'DEPEND_0': 'Epoch_esa',
                                                   'DEPEND_1': 'Pitch_Angle',
                                                   'DEPEND_2': 'Energy',
                                                   'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                   'UNITS': '!m!U-6!N s!U3!',
                                                   'VALIDMIN': distFunc.min(), 'VALIDMAX': distFunc.max(),
                                                   'VAR_TYPE': 'data', 'SCALETYP': 'log'}]}}

        # --- delete output file if it already exists ---
        if os.path.exists(outputPath):
            os.remove(outputPath)

        # --- open the output file ---
        with pycdf.CDF(outputPath, '') as sciFile:
            sciFile .readonly(False)

            # --- write out global attributes ---
            inputGlobDic = L2ModelData.cdfFile.globalattsget()
            for key, val in inputGlobDic.items():
                if key in globalAttrsMod:
                    sciFile.attrs[key] = globalAttrsMod[key]
                else:
                    sciFile.attrs[key] = val

            # --- WRITE OUT DATA ---
            for varKey, varVal in data_dict.items():
                if 'Epoch' in varKey: # epoch data
                    sciFile.new(varKey, data=varVal[0], type=33)
                elif 'Function' in varKey:
                    sciFile.new(varKey, data=varVal[0], type=pycdf.const.CDF_REAL8)
                else: # other data
                    sciFile.new(varKey, data=varVal[0])

                # --- Write out the attributes and variable info ---
                for attrKey, attrVal in data_dict[varKey][1].items():
                    if attrKey == 'VALIDMIN':
                        sciFile[varKey].attrs[attrKey] = varVal[0].min()
                    elif attrKey == 'VALIDMAX':
                        sciFile[varKey].attrs[attrKey] = varVal[0].max()
                    elif attrVal != None:
                        sciFile[varKey].attrs[attrKey] = attrVal

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

if len(glob(f'{rocketFolderPath}L2\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        Distribution_Function(wRocket, 0, rocketFolderPath, justPrintFileNames, wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}L2\{fliers[wflyer]}\*.cdf')))):
            Distribution_Function(wRocket, fileNo, rocketFolderPath, justPrintFileNames, wflyer)
    else:
        for filesNo in wFiles:
            Distribution_Function(wRocket, filesNo, rocketFolderPath, justPrintFileNames, wflyer)