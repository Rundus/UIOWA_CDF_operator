# --- template.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:



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
justPrintFileNames = True

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

modifier = ''
inputPath_modifier = 'mag' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science/magSpectrograms' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# if using ENU coordniates
useENUcoordinates = False

##################
# FILTER TOGGLES #
##################
filtOrder = 1 # order of filter
cutoff_Freq = 1 # cut off frequency where gain has reached -3dB
dataSampleFreq = 128 # sample rate of the data

# Plot the filtered data
plotFilteredData = True
plotFreqSpectrogram = True

# output the data
outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.signal import butter, filtfilt, spectrogram
def butter_highpass(cutoff, fs, order):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a
def butter_highpass_filter(data, cutoff, fs, order):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def magSpectrograms(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'Science'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['Tesseract', 'RingCore']):
        if instr.lower() in dataFile_name.lower():
            wInstr = [index, instr, instr]

    fileoutName = f'ACESII_{rocketID}_{wInstr[1]}_filtered.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict = loadDictFromFile(inputFiles[wFile],{})
        data_dict['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][i]) for i in (range(len(data_dict['Epoch'][0])))])
        Done(start_time)

        ##################################
        # --- FITLER AND PLOT MAG DATA ---
        ##################################

        # format: [Bx,By,Bz]
        if useENUcoordinates:
            magAxes = ['B_east', 'B_north', 'B_up']
        else:
            magAxes = ['Bx', 'By', 'Bz']

        rawData = np.array([data_dict[axes][0] for axes in magAxes])
        conditionedData = np.array([butter_highpass_filter(data_dict[axes][0], cutoff_Freq, dataSampleFreq, filtOrder) for axes in magAxes])

        time = data_dict['Epoch'][0]

        if plotFreqSpectrogram:
            wAxes = 0
            fig, ax = plt.subplots()
            f, t, Sxx = spectrogram(rawData[wAxes], dataSampleFreq, return_onesided=True)
            cmap = ax.pcolormesh(t, f, Sxx, cmap='turbo', vmin=1E-3, vmax=1E3, norm='log')
            ax.set_ylabel('Frequency [Hz]')
            ax.set_xlabel('Time')
            ax.set_title(f'{magAxes[wAxes]}_filtered Spectrogram\nCutoff Freq: {cutoff_Freq} Hz \nOrder: {filtOrder}')
            cbar = fig.colorbar(cmap)
            cbar.set_label('Spectral Density')
            plt.show()


        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            # Add filtered data to data dict

            for i,data in enumerate(conditionedData):
                data_dict = {**data_dict, **{f'{magAxes[i]}_filtered':
                                                 [data, {'LABLAXIS': f'{magAxes[i]}_filtered',
                                                                        'DEPEND_0': 'Epoch',
                                                                        'DEPEND_1': None,
                                                                        'DEPEND_2': None,
                                                                        'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                                        'UNITS': 'nT',
                                                                        'VALIDMIN': data.min(),
                                                                        'VALIDMAX': data.max(),
                                                                        'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
            # --- --- --- --- --- --- ---
            # --- WRITE OUT THE DATA ---
            # --- --- --- --- --- --- ---
            prgMsg('Creating output file')

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            globalAttrsMod['Descriptor'] = wInstr[1]
            outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, wInstr[1])

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
        magSpectrograms(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            magSpectrograms(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            magSpectrograms(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)