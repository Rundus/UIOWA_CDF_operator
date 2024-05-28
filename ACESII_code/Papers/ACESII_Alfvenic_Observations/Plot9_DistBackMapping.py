# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing field-aligned
# particle data along with electric and magnetic signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools
import math

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from ACESII_code.class_var_func import EpochTo_T0_Rocket, q0,m_e,Re
import matplotlib.gridspec as gridspec
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---

# TODO: Need to implement an averaging for the Plasma Sheet/Inverted-V collection where you divide by the
# length of the number of times you have a non-zero distribution function value at specific energy all over a range of times.


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
print(color.UNDERLINE + f'Plot8_Conjugacy' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
# Physics Toggles
disR_TargetTime = [dt.datetime(2022,11, 20, 17, 24, 38, 000000), dt.datetime(2022,11,20,17,25,3,000000)]
targetVar = [disR_TargetTime, 'Epoch']
plasmaSheet_targetTimes = [dt.datetime(2022, 11, 20, 17, 24, 40, 000000), dt.datetime(2022,11,20,17,24,44,000000)]
PS_additionalTime = [dt.datetime(2022, 11, 20, 17, 24, 0, 000000), dt.datetime(2022,11,20,17,24,30,000000)]
invertedV_targetTimes = [dt.datetime(2022,11, 20, 17, 25, 1, 396000), dt.datetime(2022,11,20,17,25,1,596000)]
altsMapped = [400, 2000, 5000] # altitude [in kilometers] to determine the distribution functions at higher altitudes
countNoiseLevel = 1


# --- Plot - Slice in Pitch toggles ---
figure_width = 10 # in inches
figure_height =14 # in inches
Label_FontSize = 10
Label_Padding = 8
Tick_FontSize = 12
Tick_Length = 1
Tick_Width = 1
Tick_FontSize_minor = 10
Tick_Length_minor = 1
Tick_Width_minor = 1
Plot_LineWidth = 0.5
plot_MarkerSize = 14
legend_fontSize = 10
dpi = 800
distLimits = [1E-19, 1E-12]
xAxis_energyLimits = [10, 1E4]
cbarMin, cbarMax = 1E-18, 1E-14
cbarTickLabelSize = 14
cmap = apl_rainbow_black0_cmap()
cmap.set_bad(color=(1, 1, 1))

# Pitch Slice Toggles
Plot_sliceInPitch = True
wPtch_slice = 4

# FullDistribution Toggles
Plot_fullDistribution = True
ptchRange = [i for i in range(11)]
X_Velocity_limits, Y_Velocity_limit = [-0.5, 1.8], [-0.5, 2]



##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
prgMsg('Loading Data')

rocketAttrs, b, c = ACES_mission_dicts()

# EEPAA Distribution Data
inputEEPAA_high = glob('C:\Data\ACESII\L1\high\*eepaa_fullCal*')[0]
data_dict_dist_high = loadDictFromFile(inputFilePath=inputEEPAA_high, targetVar=targetVar, wKeys_Reduce=['eepaa', 'Epoch'])

# Define the data
counts = deepcopy(data_dict_dist_high['eepaa'][0])
Epoch = deepcopy(data_dict_dist_high['Epoch'][0])
Energy = deepcopy(data_dict_dist_high['Energy'][0])
Pitch = deepcopy(data_dict_dist_high['Pitch_Angle'][0])

Done(start_time)

##############################
# --- --- --- --- --- --- ----
# --- CALCULATED VARIABLES ---
# --- --- --- --- --- --- ----
##############################
prgMsg('Mapping Data')
# Time since launch
LaunchDateTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 000000))
Epoch_timeSince = EpochTo_T0_Rocket(InputEpoch=Epoch, T0=LaunchDateTime)
plasmaSheet_timeSince = EpochTo_T0_Rocket(InputEpoch=plasmaSheet_targetTimes, T0=LaunchDateTime)
invertedV_timeSince = EpochTo_T0_Rocket(InputEpoch=invertedV_targetTimes, T0=LaunchDateTime)

regionTimes = [plasmaSheet_targetTimes,invertedV_targetTimes]
mappedCounts_output = [[], []]


# loop over altitudes
for alt in altsMapped:

    # loop over plasma sheet or inverted-V altitudes
    for regionIdx in range(len(regionTimes)):

        rTime = regionTimes[regionIdx]
        low_idx, high_idx = np.abs(Epoch - rTime[0]).argmin(), np.abs(Epoch - rTime[1]).argmin()
        countsData = counts[low_idx:high_idx]
        mappedCounts = [[[[] for engy in range(len(Energy))] for ptch in range(len(Pitch))] for tme in range(len(counts[low_idx:high_idx]))]

        # loop over all values in the data and re-bin them to a new pitch angle position based on altitude
        sizes = [len(countsData), len(countsData[0]), len(countsData[0][0])]
        ranges = [range(sizes[0]), range(sizes[1]), range(sizes[2])]

        for tme, ptch,engy in itertools.product(*ranges):

            currPitch = Pitch[ptch]
            mappedPitch = np.degrees(np.arcsin(  np.sin(np.radians(currPitch))*np.power((Re + 400)/(Re + alt), 1.5)))

            # determine where new mapped pitch should go
            newPtchIdx = np.abs(Pitch - mappedPitch).argmin()
            mappedCounts[tme][newPtchIdx][engy].append(countsData[tme][ptch][engy])


        # reduce the lists by summing the total counts at each point
        for tme, ptch, engy in itertools.product(*ranges):

            val = np.array(mappedCounts[tme][ptch][engy])

            # remove fillvals from data
            val = val[np.where(val > 0)]

            if len(val) == 0:
                mappedCounts[tme][ptch][engy] = 0
            else:
                mappedCounts[tme][ptch][engy] = sum(val)

        # Calculate the Distribution Function for each point
        geo_factor = rocketAttrs.geometric_factor[0]
        count_interval = 0.8992E-3

        for tme, ptch, engy in itertools.product(*ranges):
            deltaT = (count_interval) - (mappedCounts[tme][ptch][engy] * rocketAttrs.deadtime[0])
            fluxVal = (mappedCounts[tme][ptch][engy]) / (geo_factor[0] * deltaT)
            mappedCounts[tme][ptch][engy] = ((cm_to_m * m_e / (q0 * Energy[engy])) ** 2) * (fluxVal / 2)


        mappedCounts_output[regionIdx].append(np.array(mappedCounts))

Done(start_time)




# --- Determine the 1-count threshold line ---
distFunc_NoiseCount = np.zeros(shape=(len(Energy)))
geo_factor = rocketAttrs.geometric_factor[0]
count_interval = 0.8992E-3
for engy in range(len(Energy)):
    deltaT = (count_interval) - (countNoiseLevel * rocketAttrs.deadtime[0])
    fluxVal = (countNoiseLevel) / (geo_factor[0] * deltaT)
    distFunc_NoiseCount[engy] = ((cm_to_m*m_e / (q0*Energy[engy]))**2) * (fluxVal/2)


# --- Calculate the average distribution function value for each pitch at each energy over a range of times
paraDist = [[[], []] for alt in altsMapped ]
paraEngy = [[[], []] for alt in altsMapped]
paraPtch = [[[], []] for alt in altsMapped]

for altidx, alt in enumerate(altsMapped):

    for regionIdx in range(len(paraDist[0])):
        rTime = regionTimes[regionIdx]
        low_idx, high_idx = np.abs(Epoch - rTime[0]).argmin(), np.abs(Epoch - rTime[1]).argmin()
        tempDist = []
        tempEngy = []
        tempPtch = []

        for ptch in ptchRange:

            data = mappedCounts_output[regionIdx][altidx][:, ptch, :].T

            average = np.array([np.average(arr[np.nonzero(arr)]) for arr in data]) # average over time for each energy bin
            average = np.array([val if np.isnan(val) == False else 0 for val in average ]) # remove the points where there was no distribution function for a specific enregy

            # remove all the zero values from average and from energy
            # zeroFinder = np.where(average!=0)
            # tempDist.append(average[zeroFinder])
            # tempEngy.append(Energy[zeroFinder])
            tempDist.append(average)
            tempEngy.append(Energy)
            tempPtch.append(ptch)

        paraDist[altidx][regionIdx] = tempDist
        paraEngy[altidx][regionIdx] = tempEngy
        paraPtch[altidx][regionIdx] = tempPtch

Done(start_time)


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

if Plot_sliceInPitch:

    # --- INVERTED-V parallel overlap ---
    fig, ax = plt.subplots(len(altsMapped))
    fig.set_size_inches(7, 10)
    fig.suptitle(f'Pitch {Pitch[wPtch_slice]}' + '$^{\circ}$')

    for z, alt in enumerate(altsMapped):
        ax[z].set_title(f'Alt {alt} [km]')
        ax[z].plot(paraEngy[z][0][wPtch_slice], paraDist[z][0][wPtch_slice], label='Plasma Sheet Ambient', color='blue', marker='.', ms=plot_MarkerSize-7)
        ax[z].plot(paraEngy[z][1][wPtch_slice], paraDist[z][1][wPtch_slice], label='Inverted-V', color='red',marker='.', ms=plot_MarkerSize-7)
        ax[z].plot(Energy, distFunc_NoiseCount, label=f'{countNoiseLevel}-count level', color='black')

        # ax_Compare.plot(xData_fit,yData_fit,label='Kappa Fit',color='green',alpha=0.5)
        ax[z].set_yscale('log')
        ax[z].set_ylim(distLimits[0], distLimits[1])
        ax[z].set_xlim(xAxis_energyLimits[0], xAxis_energyLimits[1])
        ax[z].set_xscale('symlog')
        ax[z].grid()
        ax[z].set_ylabel(r'Dist Func.', fontsize=Label_FontSize)
        ax[z].legend(prop={'size': legend_fontSize})

        if z == len(altsMapped):
            ax[z].set_xlabel('Energy [eV]')

    plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot8_distribution_slice_{wPtch_slice}_base.png',dpi=dpi)







if Plot_fullDistribution:

    prgMsg('Calculating Vperp and Vparallel')
    countsTemp = deepcopy(counts[0][0:len(ptchRange)])

    for regionIdx, tmes in enumerate(regionTimes):

        Vperp = deepcopy(countsTemp)
        Vpara = deepcopy(countsTemp)

        for ptch in ptchRange:
            for engy in range(len(Energy)):
                Emag = np.sqrt(2 * q0 * Energy[engy] / (m_e))
                Vperp[ptch][engy] = np.sin(np.radians(Pitch[ptch])) * Emag
                Vpara[ptch][engy] = np.cos(np.radians(Pitch[ptch])) * Emag

    Vpara, Vperp = np.array(Vpara) / (10000 * 1000), np.array(Vperp) / (10000 * 1000)
    Done(start_time)


    # --- INVERTED-V parallel overlap ---
    fig, ax = plt.subplots(len(altsMapped),2)
    fig.set_size_inches(7, 10)
    fig.suptitle('Plasma Sheet vs Inverted-V')
    for z, alt in enumerate(altsMapped):

        # ax[z, 0].set_title(f'Alt {alt} [km]')
        # ax[z, 1].set_title(f'Alt {alt} [km]')

        ax[z, 0].pcolormesh(Vperp, Vpara, paraDist[z][0], cmap='turbo', shading='nearest', norm='log', vmin=cbarMin, vmax=cbarMax)
        ax[z,0].scatter([0],[0],alpha=0, label=f'P.S. Alt {alt} [km]')
        ax[z, 1].pcolormesh(Vperp, Vpara, paraDist[z][1], cmap='turbo', shading='nearest', norm='log', vmin=cbarMin, vmax=cbarMax, label=f'InV  Alt {alt} [km]')
        ax[z, 1].scatter([0], [0], alpha=0, label=f'InV Alt {alt} [km]')

        for i in range(2):
            ax[z,i].set_xlim(X_Velocity_limits[0], X_Velocity_limits[1])
            ax[z,i].set_ylim(Y_Velocity_limit[0], Y_Velocity_limit[1])
            ax[z,i].invert_yaxis()
            ax[z,i].legend()

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
    plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot9_Mapping_Distribution_base.png', dpi=dpi)

