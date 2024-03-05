# --- fitDistributionFunction.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Combine the High Flyer EEPAA and LEESA data into a single, large electron dataset



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
pitches = [0,1,2,3] # indicies of pitch bins to use in plot. Nominally [0,1,2,3,4] --> -10deg to 30deg
plotDataRegion = True
# ---------------------------
outputData = True
# ---------------------------

# Primary Auroral region
targetTimes = [dt.datetime(2022,11,20,17,25,30,000000), dt.datetime(2022, 11, 20, 17, 25, 35, 500000)]

# Quiet Region
# targetTimes = [dt.datetime(2022,11,20,17,22,15,000000), dt.datetime(2022, 11, 20, 17, 22, 35, 000000)]

# before dispersive region (First STEB occurs at dt.datetime(2022, 11, 20, 17, 24, 56, 000000))
# targetTimes = [dt.datetime(2022,11,20,17,24,55,000000), dt.datetime(2022, 11, 20, 17, 24, 55, 500000)]

def fitFunct(x,n,T)


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def fitDistributionFunction(wRocket, rocketFolderPath):
    wflyer = wRocket-4

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    ModelData = L2_TRICE_Quick(wflyer)
    inputFile_superSetDist = glob(f'{rocketFolderPath}L3\distFunc\{fliers[wflyer]}\*eSuperSet*')[0]
    fileoutName_Dist = f'ACESII_{rocketID}_'

    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    prgMsg(f'Loading Electron Files Files')
    # load data here
    data_dict_eSuperSet = loadDictFromFile(inputFile_superSetDist)
    Done(start_time)

    # --- Collect the Data ---
    prgMsg('Collecting Data')
    # in the region of interest, collect all the x,y points of the data [Energy, PSD]

    timeMin, timeMax = np.abs(data_dict_eSuperSet['Epoch'][0] - targetTimes[0]).argmin(),np.abs(data_dict_eSuperSet['Epoch'][0] - targetTimes[1]).argmin()
    distFunc = data_dict_eSuperSet['Distribution_Function'][0][timeMin:timeMax+1]
    oneCountLevel = data_dict_eSuperSet['oneCountLevel'][0][timeMin:timeMax+1]
    Energy = data_dict_eSuperSet['Energy'][0]
    Pitch = data_dict_eSuperSet['Pitch_Angle'][0]
    xData = []
    yData = []
    xData_singleCount = []
    yData_singleCount = []
    ranges = [range(len(distFunc)), pitches, range(len(distFunc[0][0]))]

    for tme, ptch, engy in itertools.product(*ranges):
        xData_singleCount.append(Energy[engy])
        yData_singleCount.append(oneCountLevel[tme][ptch][engy])
        if distFunc[tme][ptch][engy] != 0 and distFunc[tme][ptch][engy] != rocketAttrs.epoch_fillVal:
            xData.append(Energy[engy])
            yData.append(distFunc[tme][ptch][engy])


    # prepare the data
    xData = np.array(xData)
    yData = np.array(yData)
    xData_singleCount = np.array(xData_singleCount)
    yData_singleCount = np.array(yData_singleCount)
    Done(start_time)

    if plotDataRegion:
        prgMsg('Making Plot')

        fig, ax = plt.subplots()
        ax.scatter(xData, yData, label='EEPAA/LEESA PSD')
        # ax.plot(xData_singleCount,yData_singleCount)
        ax.set_title('Electron PSD [Observation]')
        # ax.set_title(rf'PSD for $\alpha=$ {Pitch[pitches[0]]}' +r'$^{\circ}$ to ' + f'{Pitch[pitches[-1]]}' +'$^{\circ}$')
        ax.text(x=100,y=0.8*yData.max(),
                s=f'pa < {Pitch[pitches[-1]]} deg'+ f'{targetTimes[0].strftime("%H:%M:%S")} UTC - {targetTimes[1].strftime("%H:%M:%S")}',
                ha='center',
                color='black',
                fontsize=12)
        ax.set_yscale('log')
        ax.set_ylabel('PSD [$s^{3}/m^{6}$]')
        ax.set_xscale('log')
        ax.set_xlabel('Energy [eV]')
        ax.legend()
        plt.show()
        Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
fitDistributionFunction(4, rocketFolderPath)
