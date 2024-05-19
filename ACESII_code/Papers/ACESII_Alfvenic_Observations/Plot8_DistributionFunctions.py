# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing field-aligned
# particle data along with electric and magnetic signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from ACESII_code.myImports import *
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from ACESII_code.class_var_func import EpochTo_T0_Rocket, q0,m_e
import matplotlib.gridspec as gridspec
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
print(color.UNDERLINE + f'Plot8_Conjugacy' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
# Physics Toggles
disR_TargetTime = [dt.datetime(2022,11, 20, 17, 24, 54, 000000), dt.datetime(2022,11,20,17,25,3,000000)]
targetVar = [disR_TargetTime, 'Epoch']
# plasmaSheet_targetTimes = [dt.datetime(2022,11, 20, 17, 24, 54, 762000), dt.datetime(2022,11,20,17,24,56,212000)]
plasmaSheet_targetTimes = [dt.datetime(2022,11, 20, 17, 24, 55, 000000), dt.datetime(2022,11,20,17,24,55,900000)]
invertedV_targetTimes = [dt.datetime(2022,11, 20, 17, 25, 1, 512000), dt.datetime(2022,11,20,17,25,1,570000)]
# ptchRange = [2 + i for i in range(17)] # determines which pitchs to include in the parallel distribution functions
ptchRange = [2,17] # determines which pitchs to include in the parallel distribution functions


# --- Cbar ---
cbarMin, cbarMax = 1E-20, 1E-14
cbarTickLabelSize = 14
cmap = apl_rainbow_black0_cmap()
cmap.set_bad(color=(1, 1, 1))

# --- Plot toggles ---
normalizeByTime = True
figure_width = 10 # in inches
figure_height =14# in inches
Label_FontSize = 20
Tick_FontSize = 5
Tick_Length = 1
Tick_Width = 1
Tick_FontSize_minor = 2
Tick_Length_minor = 1
Tick_Width_minor = 1
Plot_LineWidth = 0.5
Label_Padding = 8
dpi = 800
distLimits = [1E-20, 1E-14]
xAxis_energyLimits = [28, 3E3]



##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
prgMsg('Loading Data')

rocketAttrs, b, c = ACES_mission_dicts()

# EEPAA Distribution Data
inputEEPAA_low = glob('C:\Data\ACESII\L3\DistFunc\low\*ACESII_36364_distFunc_eepaa_fullCal*')[0]
data_dict_dist_low = loadDictFromFile(inputFilePath=inputEEPAA_low, targetVar=targetVar, wKeys_Reduce=['Distribution_Function', 'Epoch', 'Alt'])
inputEEPAA_high = glob('C:\Data\ACESII\L3\DistFunc\high\*ACESII_36359_distFunc_eepaa_fullCal*')[0]
data_dict_dist_high = loadDictFromFile(inputFilePath=inputEEPAA_high, targetVar=targetVar, wKeys_Reduce=['Distribution_Function', 'Epoch', 'Alt'])

# Define the data
distFunc = deepcopy(data_dict_dist_high['Distribution_Function'][0])
Epoch = deepcopy(data_dict_dist_high['Epoch'][0])
Energy = deepcopy(data_dict_dist_high['Energy'][0])
Pitch = deepcopy(data_dict_dist_high['Pitch_Angle'][0])

Done(start_time)

##############################
# --- --- --- --- --- --- ----
# --- CALCULATED VARIABLES ---
# --- --- --- --- --- --- ----
##############################
prgMsg('Preparing Data')
# Time since launch
LaunchDateTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 000000))
Epoch_timeSince = EpochTo_T0_Rocket(InputEpoch=Epoch, T0=LaunchDateTime)
plasmaSheet_timeSince = EpochTo_T0_Rocket(InputEpoch=plasmaSheet_targetTimes, T0=LaunchDateTime)
invertedV_timeSince = EpochTo_T0_Rocket(InputEpoch=invertedV_targetTimes, T0=LaunchDateTime)

# averaged distribution function over pitch angle
omniDist_high = np.zeros(shape=(len(distFunc), len(Energy)))
for tme in range(len(Epoch)):
    for engy in range(len(Energy)):
        sumVal = 0

        for ptch in ptchRange:
            val = distFunc[tme, ptch, engy]
            if val > 0:
                sumVal += val

        # Average the Omni-flux by the number of bins. ONLY include bins 10deg - 170 since they have full coverage
        omniDist_high[tme][engy] = sumVal / len(range(2, 18 + 1))



# --- PLASMASHEET - parallel distribution function over pitch angle ---
low_idx, high_idx = np.abs(Epoch-plasmaSheet_targetTimes[0]).argmin(), np.abs(Epoch-plasmaSheet_targetTimes[1]).argmin()
parallel_Dist_plasmaSheet = []
parallel_energy_plasmaSheet = []
for tme in range(len(Epoch[low_idx:high_idx])):
    for engy in range(len(Energy)):
        for ptchIdx in ptchRange:
            val = distFunc[tme][ptchIdx][engy]*np.abs(np.cos(np.radians(Pitch[ptchIdx])))
            modifier = -1 if Pitch[ptchIdx] >= 90 else 1
            parEngy = modifier*Energy[engy] * (np.cos(np.radians(Pitch[ptchIdx])))**2

            if val > 0:
                parallel_Dist_plasmaSheet.append(val)
                parallel_energy_plasmaSheet.append(parEngy)


# Sort the data based on velocity
paraEngy_PS,paraDist_PS = zip(*sorted(zip(parallel_energy_plasmaSheet,parallel_Dist_plasmaSheet)))
paraEngy_PS, paraDist_PS = np.array(paraEngy_PS),np.array(paraDist_PS)

# sum all points at the same velocity
uniqueEneries = set(paraEngy_PS)
uniqueDist = []
for engyIdx,uEngy in enumerate(uniqueEneries):
    indicies = np.where(paraEngy_PS==uEngy)[0]

    if normalizeByTime:
        uniqueDist.append(sum(paraDist_PS[indicies]) / len(Epoch[low_idx:high_idx]))
    else:
        uniqueDist.append(sum(paraDist_PS[indicies]))

paraEngy_PS = np.array(list(uniqueEneries))
paraDist_PS = np.array(uniqueDist)





# --- INVERTED-V parallel distribution function over pitch angle ---
low_idx, high_idx = np.abs(Epoch - invertedV_targetTimes[0]).argmin(), np.abs(Epoch - invertedV_targetTimes[1]).argmin()
parDist_InV = []
parEngy_InV = []
for tme in range(len(Epoch[low_idx:high_idx])):
    for engy in range(len(Energy)):
        for ptchIdx in ptchRange:
            val = distFunc[tme][ptchIdx][engy] * np.abs(np.cos(np.radians(Pitch[ptchIdx])))
            modifier = -1 if Pitch[ptchIdx] >= 90 else 1
            parEngy = modifier * Energy[engy] * (np.cos(np.radians(Pitch[ptchIdx]))) ** 2

            if val > 0:
                parDist_InV.append(val)
                parEngy_InV.append(parEngy)

# Sort the data based on velocity
parEngy_InV, parDist_InV = zip(*sorted(zip(parEngy_InV,parDist_InV)))
parEngy_InV, parDist_InV = np.array(parEngy_InV),np.array(parDist_InV)
# sum all points at the same velocity
uniqueEnergies = set(parEngy_InV)
uniqueDist = []
runningSum = 0

for engyIdx, uEngy in enumerate(uniqueEnergies):
    indicies = np.where(parEngy_InV == uEngy)[0]
    runningSum += len(indicies)
    if normalizeByTime:
        uniqueDist.append(sum(parDist_InV[indicies])/len(Epoch[low_idx:high_idx]))
    else:
        uniqueDist.append(sum(parDist_InV[indicies]))

parEngy_InV = np.array(list(uniqueEnergies))
parDist_InV = np.array(uniqueDist)
Done(start_time)

############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################
prgMsg('Plotting Figure')

# DETERMINE THE PLOT DIMESNIONS
fig = plt.figure()
fig.set_size_inches(figure_width, figure_height)
gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[0.25, 0.75], hspace=0.12)

# Top Pllot
gs00 = gs0[0].subgridspec(1, 1)
axDistOverview = fig.add_subplot(gs00[0])

# Subplots
gs01 = gs0[1].subgridspec(2, 2)
ax_PS_ParDist = fig.add_subplot(gs01[0, 0])
ax_InV_ParDist = fig.add_subplot(gs01[0, 1])
axS3 = fig.add_subplot(gs01[1, 0])
axS4 = fig.add_subplot(gs01[1, 1])
subAxes = [ax_PS_ParDist, ax_InV_ParDist, axS3, axS4]



# ---Distribution Data Overview---
cmap = axDistOverview.pcolormesh(Epoch_timeSince, Energy, omniDist_high.T, cmap=cmap,vmin=cbarMin,vmax=cbarMax, norm='log')
axDistOverview.set_ylabel('Energy [eV]', fontsize=Label_FontSize)
axDistOverview.tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
axDistOverview.tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize_minor, length=Tick_Length_minor, width=Tick_Width_minor)
axDistOverview.set_yscale('log')
axDistOverview.set_title(r'Distribution Function, $\alpha = 10^{\circ}$')

axDistOverview.axvspan(invertedV_timeSince[0], invertedV_timeSince[1], color='gray', alpha=0.25,zorder=1)
axDistOverview.axvspan(plasmaSheet_timeSince[0], plasmaSheet_timeSince[1], color='gray', alpha=0.25,zorder=1)


# --- PLASMA SHEET - PARALLEL DISTRIBUTION ---
ax_PS_ParDist.set_title('"Plasma Sheet" Distribution',fontsize=Label_FontSize)
ax_PS_ParDist.scatter(paraEngy_PS,paraDist_PS)
ax_PS_ParDist.set_yscale('log')
ax_PS_ParDist.set_ylim(distLimits[0],distLimits[1])
ax_PS_ParDist.set_xlim(xAxis_energyLimits[0],xAxis_energyLimits[1])
ax_PS_ParDist.set_xscale('symlog')
ax_PS_ParDist.grid()

# --- INVERTEDV - PARALLEL DISTRIBUTION ---
ax_InV_ParDist.set_title('Inverted-V Distribution',fontsize=Label_FontSize)
ax_InV_ParDist.scatter(parEngy_InV, parDist_InV)
ax_InV_ParDist.set_yscale('log')
ax_InV_ParDist.set_ylim(distLimits[0],distLimits[1])
ax_InV_ParDist.set_xlim(xAxis_energyLimits[0],xAxis_energyLimits[1])
ax_InV_ParDist.set_xscale('symlog')
ax_InV_ParDist.grid()

# --- INVERTED-V parallel overlap ---
axS4.scatter(paraEngy_PS,paraDist_PS,label='Plasma Sheet Ambient',color='blue')
axS4.scatter(parEngy_InV,parDist_InV,label='Inverted-V',color='red')
axS4.set_yscale('log')
axS4.set_ylim(distLimits[0],distLimits[1])
axS4.set_xlim(xAxis_energyLimits[0],xAxis_energyLimits[1])
axS4.set_xscale('symlog')
axS4.grid()
axS4.set_title('Comparison',fontsize=Label_FontSize)
axS4.set_ylabel(r'Distribution Function, $\alpha = 10^{\circ}$',fontsize=Label_FontSize)


# --- cbar ---
cax = fig.add_axes([0.91, 0.288, 0.02, 0.592])
cbar = plt.colorbar(cmap, cax=cax)
# cbar.set_label('Omni-Dir. diff E. Flux \n' + '[cm$^{-2}$str$^{-1}$eV/eV]', rotation=-90, labelpad=20, fontsize=General_LabelFontSize)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbarTickLabelSize + 5)


fig.subplots_adjust(left=None, bottom=0.3, right=None, top=0.99, wspace=None, hspace=None)  # remove the space between plots
plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot8_DistributionFunc_base.png', dpi=dpi)
Done(start_time)

