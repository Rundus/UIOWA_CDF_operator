# --- AlfvenSpeed_rkt.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Determine the Alfven speed over the flight of the payloads using
# the density data and magnetometer data



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt

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
wFiles = []

modifier = ''
inputPath_modifier_density = 'l3/Langmuir' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
wFile_ni = 0
inputPath_modifier_mag = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
wFile_mag = 0
inputPath_modifier_deltaE = 'l3/deltaE'
wFile_deltaE =0
inputPath_modifier_deltaB = 'l3/deltaB'
wFile_deltaB =0
outputPath_modifier = 'science/AlfvenSpeed_rkt' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
outputData = True


# --- reduce data ---
reduceData = True
targetTimes = [dt.datetime(2022,11,20,17,24,30,00),dt.datetime(2022,11,20,17,25,30,00)]

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import u0,IonMasses
from numpy.fft import rfft, fftfreq

def AlfvenSpeed_rkt(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)


    inputFiles_mag = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}{modifier}\*RingCore_rktFrm*')
    inputFiles_density = glob(f'{rocketFolderPath}{inputPath_modifier_density}\{fliers[wflyer]}{modifier}\*fixed*.cdf')
    inputFiles_deltaE = glob(f'{rocketFolderPath}{inputPath_modifier_deltaE}\{fliers[wflyer]}{modifier}\*Field_Aligned*')
    inputFiles_deltaB = glob(f'{rocketFolderPath}{inputPath_modifier_deltaB}\{fliers[wflyer]}{modifier}\*Field_Aligned*')


    fileoutName = f'ACESII_{rocketID}_AlfvenSpeed_flight.cdf'

    if justPrintFileNames:
        print('--- MAG ---')
        for i, file in enumerate(inputFiles_mag):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_mag[i], round(getsize(file) / (10 ** 6), 1)))
        print('\n--- Density ---')
        for i, file in enumerate(inputFiles_density):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_density[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Calculating Alfven Speed' + color.END)

        # --- get the data from the Bmag file ---
        prgMsg(f'Loading data from mag File')
        data_dict_mag = loadDictFromFile(inputFiles_mag[wFile_mag],reduceData=reduceData,targetTimes=targetTimes)
        Done(start_time)

        # --- get the data from the density file ---
        prgMsg(f'Loading data from density File')
        data_dict_density = loadDictFromFile(inputFiles_density[wFile_ni],reduceData=reduceData,targetTimes=targetTimes)
        Done(start_time)

        # --- get the data from the deltaB file ---
        prgMsg(f'Loading data from deltaB File')
        data_dict_deltaB = loadDictFromFile(inputFiles_deltaB[wFile_deltaB], reduceData=reduceData, targetTimes=targetTimes)
        Done(start_time)

        # --- get the data from the deltaE file ---
        prgMsg(f'Loading data from deltaE File')
        data_dict_deltaE = loadDictFromFile(inputFiles_deltaE[wFile_deltaE], reduceData=reduceData, targetTimes=targetTimes)
        Done(start_time)


        # reduce data
        lowCut, highCut = np.abs(data_dict_mag['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_mag['Epoch'][0] - targetTimes[1]).argmin()

        for key, val in data_dict_mag.items():
            data_dict_mag[key][0] = data_dict_mag[key][0][lowCut:highCut]

        lowCut, highCut = np.abs(data_dict_density['Epoch'][0] - targetTimes[0]).argmin(), np.abs( data_dict_density['Epoch'][0] - targetTimes[1]).argmin()
        for key, val in data_dict_density.items():
            data_dict_density[key][0] = data_dict_density[key][0][lowCut:highCut]


        ###################################
        # --- DOWNSAMPLE DENSITY TO MAG ---
        ###################################
        prgMsg('Downsampling Density')
        from ACESII_code.class_var_func import InterpolateDataDict,dateTimetoTT2000

        data_dict_densityInterp = InterpolateDataDict(InputDataDict=data_dict_density,
                                                      InputEpochArray=dateTimetoTT2000(data_dict_density['Epoch'][0], False),
                                                      wKeys=['ni', 'Epoch'],
                                                      targetEpochArray=dateTimetoTT2000(data_dict_mag['Epoch'][0],False))
        Done(start_time)

        if wRocket == 5:
            prgMsg('Downsampling deltaE on deltaB')

            data_dict_deltaE = InterpolateDataDict(InputDataDict=data_dict_deltaE,
                                                          InputEpochArray=dateTimetoTT2000(data_dict_deltaE['Epoch'][0], False),
                                                          wKeys=['E_e','E_p','E_r', 'Epoch'],
                                                          targetEpochArray=dateTimetoTT2000(data_dict_deltaB['Epoch'][0], False))


            Done(start_time)

        #############################
        # --- SMOOTH DENSITY DATA ---
        #############################

        from scipy.signal import savgol_filter
        # densityVar =  savgol_filter(x=data_dict_densityInterp['ni'][0],
        #                             window_length=1000,
        #                             polyorder=5)

        data_dict_densityInterp['ni'][0] = savgol_filter(x=data_dict_densityInterp['ni'][0], window_length=1000, polyorder=5)




        ####################################
        # --- CALCULATE THE ALFVEN SPEED ---
        ####################################

        # MHD Alfven Speed - Plasma Density and Magnetic Field
        AlfvenSpeed_MHD = np.array([
            ((1E-9)*data_dict_mag['Bmag'][0][i])/np.sqrt(u0*(100**3)*data_dict_densityInterp['ni'][0][i]*IonMasses[0]) for i in range(len(data_dict_mag['Epoch'][0]))
        ])

        if wRocket == 5:
            # Alfven Speed - dEr/dBe
            AlfvenSpeed_ErBe = (1/1E-6)*np.array(data_dict_deltaE['E_r'][0])/np.array(data_dict_deltaB['B_e'][0])

            # Alfven Speed - -dEe/dBr
            AlfvenSpeed_EeBr = -1*(1/1E-6)*np.array(data_dict_deltaE['E_e'][0])/np.array(data_dict_deltaB['B_r'][0])

            # FFT of Br and -Ee
            N, T = len(data_dict_deltaB['Epoch'][0]), 1 / 128
            FFT_freqs = fftfreq(N, T)[:N // 2]

            yf_Br = rfft(data_dict_deltaB['B_r'][0])
            FFT_Br = 2.0 / N * np.abs(yf_Br[0:N // 2])

            yf_Ee = rfft(data_dict_deltaE['E_e'][0])
            FFT_Ee = 2.0 / N * np.abs(yf_Ee[0:N // 2])

            AlfvenSpeed_FFT = 1E6*FFT_Ee/FFT_Br


        # fig, ax = plt.subplots(2)
        # ax[0].plot(FFT_freqs, FFT_Br,color='blue')
        # ax[0].plot(FFT_freqs, FFT_Ee, color='red')
        # ax[1].plot(FFT_freqs, 1E6*AlfvenSpeed_FFT)
        # ax[1].set_ylim(1E5,3E8)
        # for i in range(2):
        #     ax[i].set_xlim(0, 20)
        #     ax[i].set_yscale('log')
        # plt.show()

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('creating output file')

            # create output data_dict
            data_dict_output = {}

            #prepare the alfven data
            AlfvenData = [AlfvenSpeed_MHD, AlfvenSpeed_ErBe, AlfvenSpeed_EeBr, AlfvenSpeed_FFT, FFT_freqs, FFT_Br, FFT_Ee]
            epochNames = ['Epoch', 'Epoch_delta', 'Epoch_delta', 'FFT_freqs', None,'FFT_freqs','FFT_freqs']
            varNames = ['MHD', 'ErBe', 'mEeBr', 'FFT_division','FFT_freqs','FFT_Br','FFT_Ee']
            for i, data in enumerate(AlfvenData):

                if i==1 and wRocket == 4:
                    break

                if i == 4:
                    lablaxis = 'Frequency'
                    units = 'Hz'
                elif i == 5:
                    lablaxis = 'B_r'
                    units = 'nT'
                elif i == 6:
                    lablaxis = 'E_e'
                    units = 'mV/m'
                else:
                    lablaxis = 'Alfven Speed'
                    units = 'm/s'

                varAttrs = {'LABLAXIS': lablaxis, 'DEPEND_0': epochNames[i], 'DEPEND_1': None, 'DEPEND_2': None,
                            'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': units,
                            'VALIDMIN': data.min(), 'VALIDMAX': data.max(),
                            'VAR_TYPE': 'data', 'SCALETYP': 'linear'}

                data_dict_output = {**data_dict_output, **{varNames[i]: [data, varAttrs]}}

            # - prepare the epoch variables -
            Epoch_output = deepcopy(data_dict_mag['Epoch'])
            Epoch_output[1]['VAR_TYPE'] = 'support_data'
            data_dict_output = {**data_dict_output, **{'Epoch': Epoch_output}}

            if wRocket == 5:
                Epoch_output = deepcopy(data_dict_deltaB['Epoch'])
                Epoch_output[1]['VAR_TYPE'] = 'support_data'
                data_dict_output = {**data_dict_output, **{'Epoch_delta': Epoch_output}}

            #output the data
            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'
            outputCDFdata(outputPath, data_dict_output, outputModelData, globalAttrsMod, 'Attitude')
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

if len(glob(f'{rocketFolderPath}{inputPath_modifier_density}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        AlfvenSpeed_rkt(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        AlfvenSpeed_rkt(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)