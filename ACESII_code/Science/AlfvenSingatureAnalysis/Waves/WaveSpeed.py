# --- PoyntingFlux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Determine the PoyntingFLux of the data using E-Field and B-Field Measurements.
# For the low flyer, it ONLY accepts despun data, high flyer has its own case



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
printMagFiles = False
printElecFiles = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

modifier = ''
inputPath_modifier_elec = 'science/deltaE'
wMagFile = 1
inputPath_modifier_mag = 'science/deltaB' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
wEFIFile = 1
outputPath_modifier = 'science/AlfvenSignatureWaveSpeed' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder


# --- --- --- Which Data --- -- ---
useDelta_E_B = True # use the deltaB, deltaE data
# --- --- --- reduce data --- --- ---
reduceData = True
targetTimes = [pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 24, 30, 000)),
               pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 25, 15, 000))]

# --- --- --- PLOT --- --- ---
plotSPoynting = False

# --- --- --- OUTPUT --- --- ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---

from ACESII_code.class_var_func import u0, IonMasses, InterpolateDataDict

def WaveSpeed(wRocket, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles_elec = glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wflyer]}{modifier}\*.cdf*')
    inputFiles_mag = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}{modifier}\*RingCore*')

    fileoutName = f'ACESII_{rocketID}_Signature_WaveSpeed.cdf'


    if justPrintFileNames:
        if printMagFiles:
            print('--- B-FIELD FILES ---')
            for i, file in enumerate(inputFiles_mag):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_mag[i], round(getsize(file) / (10 ** 6), 1)))
            print('\n')

        if printElecFiles:
            print('--- E-FIELD FILES ---')
            for i, file in enumerate(inputFiles_elec):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_elec[i], round(getsize(file) / (10 ** 6), 1)))
            print('\n')
    else:
        print('\n')
        print(color.UNDERLINE + f'Calculating Wave Speed for ACESII {rocketID}' + color.END)

        # --- get the data from the mag file ---
        prgMsg(f'Loading data from mag Files')
        data_dict_mag = loadDictFromFile(inputFiles_mag[wMagFile], {})
        data_dict_mag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag['Epoch'][0]])

        # component names for the magnetic field
        compNames_mag = [ key for key, val in data_dict_mag.items() if key.lower() not in ['epoch','db_mag']]


        if reduceData:
            indicies = [np.abs(data_dict_mag['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_mag['Epoch'][0] - targetTimes[1]).argmin()]

            for key, val in data_dict_mag.items():
                data_dict_mag[key][0] = data_dict_mag[key][0][indicies[0]:indicies[1]]

        # create vector variable and convert to tesla
        B_Field = (1E-9)*np.array([np.array([data_dict_mag[compNames_mag[0]][0][i], data_dict_mag[compNames_mag[1]][0][i], data_dict_mag[compNames_mag[2]][0][i]]) for i in range(len(data_dict_mag['Epoch'][0]))])
        Done(start_time)

        if wRocket == 4:

            # collect the Magnitude of B from L1 spun data
            prgMsg('Getting Bmag')
            inputFileBmag = glob('C:\Data\ACESII\L1\high\*RingCore_rktFrm*')
            data_dict_Bmag = loadDictFromFile(inputFileBmag[0], {})
            data_dict_Bmag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_Bmag['Epoch'][0]])
            Done(start_time)

            prgMsg('Getting Plasma Density')
            inputFileBmag = glob('C:\Data\ACESII\science\Langmuir\high\*Temp&Density*')
            data_dict_density = loadDictFromFile(inputFileBmag[0], {})
            data_dict_density['fixed_Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_density['fixed_Epoch'][0]])
            Done(start_time)

            # reduce the datasets:
            prgMsg('Reducing Data')

            # Bmag
            indicies = [np.abs(data_dict_Bmag['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_Bmag['Epoch'][0] - targetTimes[1]).argmin()]
            for key, val in data_dict_Bmag.items():
                data_dict_Bmag[key][0] = deepcopy(data_dict_Bmag[key][0][indicies[0]:indicies[1]])

            # Density
            indicies = [np.abs(data_dict_density['fixed_Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_density['fixed_Epoch'][0] - targetTimes[1]).argmin()]
            for key, val in data_dict_density.items():
                data_dict_density[key][0] = deepcopy(data_dict_density[key][0][indicies[0]:indicies[1]])

            Done(start_time)

            #################################
            # --- CALCULATE POYNTING FLUX ---
            #################################
            prgMsg('Calculating Poynting Flux using EigenFunction')

            # down sample the density data onto the magnetometer data
            indiciesDownsampled = [np.abs(data_dict_density['fixed_Epoch'][0] - data_dict_mag['Epoch'][0][i]).argmin() for i in range(len(data_dict_mag['Epoch'][0]))]
            plasmaDensity = np.array([data_dict_density['fixed_ni_density'][0][index] for index in indiciesDownsampled])

            # calculate Alfven Velocity
            AlfvenVelocity = np.array([((1E-9)*data_dict_Bmag['Bmag'][0][i])/np.sqrt(u0*plasmaDensity[i]*IonMasses[0]) for i in range(len(plasmaDensity))])


            # calculate Alfven eigenfunction E

            Done(start_time)

        elif wRocket == 5:

            # --- get the data from the electric file ---
            prgMsg(f'Loading data from Electric Field Files')
            data_dict_elec = loadDictFromFile(inputFiles_elec[wEFIFile], {})
            data_dict_elec['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_elec['Epoch'][0]])
            compNames_elec = [key for key, val in data_dict_elec.items() if key.lower() not in ['epoch', 'de_mag']]

            # --- reduce each dataset ---
            if reduceData:
                indicies = [np.abs(data_dict_elec['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_elec['Epoch'][0] - targetTimes[1]).argmin()]
                for key, val in data_dict_elec.items():
                    data_dict_elec[key][0] = deepcopy(data_dict_elec[key][0][indicies[0]:indicies[1]])

            # get the electric field and convert it to V/m
            E_Field = (1/1000)*np.array([ [data_dict_elec[compNames_elec[0]][0][i], data_dict_elec[compNames_elec[1]][0][i],data_dict_elec[compNames_elec[2]][0][i]] for i in range(len(data_dict_elec['Epoch'][0]))])
            Done(start_time)

            ###########################################
            # --- INTERPOLATE EFI ONTO MAG TIMEBASE ---
            ###########################################

            prgMsg('Downsampling via Interpolate the EFI Data')
            data_dict_elec_interp = InterpolateDataDict(InputDataDict=data_dict_elec,
                                                        InputEpochArray=data_dict_elec['Epoch'][0],
                                                        wKeys=[],
                                                        targetEpochArray=data_dict_mag['Epoch'][0])

            E_Field = np.array( [  [data_dict_elec_interp[compNames_elec[0]][0][i],data_dict_elec_interp[compNames_elec[1]][0][i],data_dict_elec_interp[compNames_elec[2]][0][i]] for i in range(len(data_dict_mag['Epoch'][0]))])
            Done(start_time)

            ##############################
            # --- CALCULATE WAVE SPEED ---
            ##############################
            prgMsg('Calculating Wave Speed')
            vErB = np.array([np.linalg.norm(E_Field[i])/np.linalg.norm(B_Field[i]) for i in range(len(E_Field))])
            Done(start_time)

        Done(start_time)

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('creating output file')

            data_dict_output = {}

            data = vErB
            varAttrs = {'LABLAXIS': 'Wave Speed', 'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None,
                        'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'm/s',
                        'VALIDMIN': data.min(), 'VALIDMAX': data.max(),
                        'VAR_TYPE': 'data', 'SCALETYP': 'linear'}

            data_dict_output = {**data_dict_output, **{'Wave Speed':[data,varAttrs]}}

            if wRocket == 4:
                Epoch_output = deepcopy(data_dict_mag['Epoch'])
            elif wRocket == 5:
                Epoch_output = deepcopy(data_dict_mag['Epoch'])

            Epoch_output[1]['VAR_TYPE'] = 'support_data'

            data_dict_output = {**data_dict_output, **{'Epoch': Epoch_output}}

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict_output, outputModelData, globalAttrsMod, 'ScienceData')

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

if len(glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no electric field .cdf files in the specified directory' + color.END)
elif len(glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no B-field .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        WaveSpeed(wRocket, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        WaveSpeed(wRocket, rocketFolderPath, justPrintFileNames, wflyer)
