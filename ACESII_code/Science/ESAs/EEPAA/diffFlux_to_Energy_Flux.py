# --- diffFlux_to_Energy_Flux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: using the specs of the ACESII ESAs, convert from differential Energy Flux
# to just energy flux as described in EEPAA_Flux_Conversion.pdf document in Overleaf



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
justPrintFileNames = False
wRocket = 5
wFiles = [0]
modifier = ''
inputPath_modifier = 'l2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'l3\Energy_Flux' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# ---------------------------
outputData = True
# ---------------------------

# solid angle contributions of each pitch bin (starting at 0deg to -190 deg)
solidAngle =[
    0.189604, # -10deg
    0.0477741, # 0 deg
    0.189604, # 10 deg
    0.351294, # 20 deg
    0.54761568, # 30
    0.704001, # 40
    0.838996, # 50
    0.9485, # 60
    1.0292, # 70
    1.07859, # 80
    1.0952, # 90
    1.07859, # 100
    1.0292, # 110
    0.9485, # 120
    0.838996, # 130
    0.704001, # 140
    0.54761568, # 150 deg
    0.351294, # 160 deg
    0.189604, # 170deg
    0.0477741, # 180deg
    0.189604 # 190 deg
]



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from numpy import trapz


def diffFlux_to_Energy_Flux(wRocket, rocketFolderPath, justPrintFileNames, wflyer, wfile):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    ModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace(inputPath_modifier.lower() +'_', '') for ifile in input_names]

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in inputFiles[wfile]:
            descriptiorNam = ['EEPAA', 'LEESA', 'IEPAA', 'Langmuir']
            wInstr = [index, instr, descriptiorNam[index]]
            break
        else:
            wInstr = [0, 'attitude', '']

    fileoutName = f'ACESII_{rocketID}_{wInstr[1].lower()}_Energy_Flux.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')

        # --- --- --- --- --- -
        # --- LOAD THE DATA ---
        # --- --- --- --- --- -

        # --- get the data from the file ---
        prgMsg(f'Loading data from ESA Files')
        # load data here
        data_dict = loadDictFromFile(input_data_dict={},inputFilePath=inputFiles[wfile],targetTimes=[],reduceData=False,wKeys=[])
        Done(start_time)

        # --- --- --- --- --- --- --- -
        # --- INTEGRATE OVER ENERGY ---
        # --- --- --- --- --- --- --- -

        prgMsg('Calculating Energy Flux')

        Epoch = data_dict['Epoch'][0]
        Pitch = data_dict['Pitch_Angle'][0]
        Energy = data_dict['Energy'][0][::-1]
        diffEflux = data_dict['Differential_Energy_Flux'][0]

        Eflux = [ [[] for ptch in Pitch] for tme in Epoch]

        # perform the integration
        for tme in tqdm(range(len(Epoch))):
            for ptch in range(len(Pitch)):
                diffFLux_data = diffEflux[tme][ptch][::-1]

                if all(np.array(diffFLux_data) >= 0): # see if there's any fillvals or negative values in array. If not, then you can integrate
                    integratedValue = trapz(y=diffFLux_data, x=Energy)
                    EFluxVal = (q0 * cm_to_m ** 2) * integratedValue * solidAngle[ptch]
                else:
                    EFluxVal = rocketAttrs.epoch_fillVal


                # if integratedValue < 0:
                #     print(tme,ptch,integratedValue)
                #     print(Energy)
                #     print(diffFLux_data)

                Eflux[tme][ptch].append(EFluxVal)
                Eflux[tme][ptch] = Eflux[tme][ptch][0] # reduce each element from a list to a single value

        Eflux = np.array(Eflux)



        Done(start_time)


        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('Creating output file')

            data_dict_output = {'Energy_Flux': [np.array(Eflux), data_dict['Differential_Energy_Flux'][1]],
                                'Pitch_Angle': data_dict['Pitch_Angle'],
                                'Energy': data_dict['Energy'],
                                'Epoch': data_dict['Epoch']}
            data_dict_output['Energy_Flux'][1]['LABLAXIS'] = 'Energy_Flux'
            data_dict_output['Energy_Flux'][1]['UNITS'] = 'J m!A-2!B!Ns!A-1!B'

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict_output, ModelData, globalAttrsMod, wInstr[1])

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
        diffFlux_to_Energy_Flux(wRocket, rocketFolderPath, justPrintFileNames, wflyer, 0)
    elif wFiles == []:
        for i, wfile in enumerate(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')):
            diffFlux_to_Energy_Flux(wRocket, rocketFolderPath, justPrintFileNames, wflyer, wfile)
    else:
        for wfile in wFiles:
            diffFlux_to_Energy_Flux(wRocket, rocketFolderPath, justPrintFileNames, wflyer, wfile)
