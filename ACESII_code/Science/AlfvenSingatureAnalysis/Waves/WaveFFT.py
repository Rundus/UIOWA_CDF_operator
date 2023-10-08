# --- WaveFFT.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Make a nice plot of the FFT power in my waves



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
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
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5
modifier = ''
inputPath_modifier_elec = 'science\deltaE'
inputPath_modifier_mag = 'science\deltaB' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
plotBFFT = True
plotEFFT = True # only works for wRocket == 5


reduceData = True
targetTimes = [dt.datetime(2022, 11, 20, 17, 24, 50, 000000), dt.datetime(2022, 11, 20, 17, 25, 11, 000000)]


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from numpy.fft import rfft, fftfreq

def reduceData(targetTimes,data_dict):
    lowCutoff, highCutoff = np.abs(data_dict['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict['Epoch'][0] - targetTimes[1]).argmin()
    for key, val in data_dict.items():
        data_dict[key][0] = np.array(data_dict[key][0][lowCutoff:highCutoff])
    return data_dict


def WaveFFT(wRocket, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'

    inputFiles_elec = glob(f'{rocketFolderPath}{inputPath_modifier_elec}\low{modifier}\*Alfven*')
    inputFiles_mag_high = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\high{modifier}\*Alfven*')
    inputFiles_mag_low = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\low{modifier}\*Alfven*')

    fitResults = {
        'Bx': {'Spin Amp': 25.42873940404161, 'Spin Freq': 0.6463295881639182, 'Spin Phase': 91.9759995936283,
               'Cone Amp': 625.8772357084948, 'Cone Freq': 0.05294818121871208, 'Cone Phase': -138.77308595997619,
               'Offset': -44919.748937299344},
        'By': {'Spin Amp': 7.378420193701481, 'Spin Freq': 0.6442248190622027, 'Spin Phase': 109.20255873087793,
               'Cone Amp': 1380.5616077430786, 'Cone Freq': 0.02700105226961604, 'Cone Phase': 109.87799606103452,
               'Offset': -139.74554466082876},
        'Bz': {'Spin Amp': 8.095746809541962, 'Spin Freq': 0.6442537451458561, 'Spin Phase': 19.11852573798773,
               'Cone Amp': 1257.0313161879794, 'Cone Freq': 0.026874206798816504, 'Cone Phase': -69.78175516947503,
               'Offset': 32.456720919269245}
    }
    spinFreq = sum([fitResults['Bz']['Spin Freq'], fitResults['By']['Spin Freq'],
                    fitResults['Bz']['Spin Freq']]) / 3 if wRocket == 4 else 0.55


    print('\n')
    print(color.UNDERLINE + f'Calculating Poynting flux for ACESII {rocketID}' + color.END)

    # --- get the data from the mag file ---
    prgMsg(f'Loading data from mag Files')
    data_dict_mag_high = loadDictFromFile(inputFiles_mag_high[0], {})
    data_dict_mag_low = loadDictFromFile(inputFiles_mag_low[0], {})
    Done(start_time)

    # --- get the data from the electric file ---
    prgMsg(f'Loading data from Electric Field Files')
    data_dict_elec = loadDictFromFile(inputFiles_elec[0], {})
    Done(start_time)

    if reduceData:
        prgMsg('Reducing Data')
        data_dict_mag_high = reduceData(targetTimes, data_dict_mag_high)
        data_dict_mag_low = reduceData(targetTimes, data_dict_mag_low)
        data_dict_elec = reduceData(targetTimes, data_dict_elec)
        Done(start_time)

    B_Field_high = np.array([[data_dict_mag_high['dB_East'][0][i], data_dict_mag_high['dB_North'][0][i], data_dict_mag_high['dB_Up'][0][i]] for i in range(len(data_dict_mag_high['Epoch'][0]))])
    B_Field_low = np.array([[data_dict_mag_high['dB_East'][0][i], data_dict_mag_high['dB_North'][0][i], data_dict_mag_high['dB_Up'][0][i]] for i in range(len(data_dict_mag_high['Epoch'][0]))])
    E_Field = np.array([[data_dict_elec['dE_East'][0][i], data_dict_elec['dE_North'][0][i], data_dict_elec['dE_Up'][0][i]] for i in range(len(data_dict_elec['Epoch'][0]))])

    # --- --- --- --- --- -
    # --- PLOT THE DATA ---
    # --- --- --- --- --- -

    if plotBFFT:
        B_Field = B_Field_high if wRocket == 4 else B_Field_low
        Epoch = data_dict_mag_high['Epoch'][0] if wRocket == 4 else data_dict_mag_low['Epoch'][0]

        fig, ax = plt.subplots(nrows=3, ncols=2)

        fig.suptitle(f'ACESII {rocketID} RingCore')

        fs = 128
        N, T = len(B_Field), 1 / fs
        compNames = ['B_East', 'B_North', 'B_Up']
        for i in range(3):
            yf = rfft(B_Field[:, i])
            xf = fftfreq(N, T)[:N // 2]
            FFT = 2.0 / N * np.abs(yf[0:N // 2])

            # plot the data
            ax[i, 0].plot(Epoch, B_Field[:,i])
            ax[i, 0].set_ylabel(compNames[i] + ' [nT]')

            # plot the FFT
            ax[i, 1].set_ylabel('FFT')
            ax[i, 1].set_xlabel('Frequency [Hz]')
            ax[i, 1].plot(xf, FFT)
            ax[i, 1].vlines([spinFreq * (i + 1) for i in range(50)], ymin=min(FFT), ymax=max(FFT), alpha=0.25, color='red')
            ax[i, 1].set_xlim(0, 9)

        plt.show()

    if plotEFFT:

        fig, ax = plt.subplots(nrows=3, ncols=2)
        fig.suptitle(f'ACESII {rocketID} EFI')

        fs = 4000
        N, T = len(E_Field), 1 / fs
        compNames = ['E_East', 'E_North', 'E_Up']
        for i in range(3):
            yf = rfft(E_Field[:, i])
            xf = fftfreq(N, T)[:N // 2]
            FFT = 2.0 / N * np.abs(yf[0:N // 2])

            # plot the data
            ax[i, 0].plot(data_dict_elec['Epoch'][0], E_Field[:,i])
            ax[i, 0].set_ylabel(compNames[i] + ' [mV/m]')

            # plot the FFT
            ax[i, 1].set_ylabel('FFT')
            ax[i, 1].set_xlabel('Frequency [Hz]')
            ax[i, 1].plot(xf, FFT)
            ax[i, 1].vlines([spinFreq * (i + 1) for i in range(50)], ymin=min(FFT), ymax=max(FFT), alpha=0.25, color='red')
            ax[i, 1].set_xlim(0, 9)

        plt.show()






# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1


WaveFFT(wRocket, rocketFolderPath, justPrintFileNames, wflyer)
