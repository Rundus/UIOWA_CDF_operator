# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing field-aligned
# particle data along with electric and magnetic signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
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
plasmaSheet_targetTimes = [dt.datetime(2022,11, 20, 17, 24, 55, 446000), dt.datetime(2022,11,20,17,24,56,000000)]
invertedV_targetTimes = [dt.datetime(2022,11, 20, 17, 25, 1, 396000), dt.datetime(2022,11,20,17,25,1,596000)]
ptchRange = [2] # determines which pitchs to include in the parallel distribution functions


# --- Cbar ---
cbarMin, cbarMax = 1E-18, 1E-14
cbarTickLabelSize = 14
cmap = apl_rainbow_black0_cmap()
cmap.set_bad(color=(1, 1, 1))

# --- Plot toggles ---
normalizeByTime = True
figure_width = 10 # in inches
figure_height =14# in inches
Label_FontSize = 20
Label_Padding = 8
Tick_FontSize = 12
Tick_Length = 1
Tick_Width = 1
Tick_FontSize_minor = 10
Tick_Length_minor = 1
Tick_Width_minor = 1
Plot_LineWidth = 0.5
legend_fontSize = 15
dpi = 800
distLimits = [1E-20, 1E-14]
xAxis_energyLimits = [28, 1E4]



##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
prgMsg('Loading Data')

rocketAttrs, b, c = ACES_mission_dicts()

# EEPAA Distribution Data
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


# --- Get the Distribution Function data as Parallel Distribution vs Energy [eV] ---
paraDist = [[],[]]
paraEngy = [[],[]]
regionTimes = [plasmaSheet_targetTimes,invertedV_targetTimes]


for regionIdx in range(len(paraDist)):

    rTime = regionTimes[regionIdx]
    low_idx, high_idx = np.abs(Epoch - rTime[0]).argmin(), np.abs(Epoch - rTime[1]).argmin()
    tempDist = []
    tempEngy = []

    # collect all the data in the pitch angle range here and store in in temp variables
    for tme in range(low_idx,high_idx):
        for engy in range(len(Energy)):
            for ptchIdx in ptchRange:
                val = distFunc[tme][ptchIdx][engy] * np.abs(np.cos(np.radians(Pitch[ptchIdx])))
                modifier = -1 if Pitch[ptchIdx] > 90 else 1
                parallel_value = modifier * Energy[engy] * (np.cos(np.radians(Pitch[ptchIdx]))) ** 2

                if val > 0:
                    tempDist.append(val)
                    tempEngy.append(parallel_value)

    # Sort the data based on velocity
    tempEngy_sorted, tempDist_sorted = zip(*sorted(zip(tempEngy, tempDist)))
    tempEngy_sorted, tempDist_sorted = np.array(tempEngy_sorted), np.array(tempDist_sorted)


    # sum all points at the same velocity
    uEngies = set(tempEngy_sorted)
    uDist = []
    for engyIdx, uVal in enumerate(uEngies):
        indicies = np.where(tempEngy_sorted==uVal)[0]

        if normalizeByTime:
            uDist.append(sum(tempDist_sorted[indicies]) / len(Epoch[low_idx:high_idx]))
        else:
            uDist.append(sum(tempDist_sorted[indicies]))

    paraDist[regionIdx] = np.array(uDist)
    paraEngy[regionIdx] = np.array(list(uEngies))

Done(start_time)



###############################
# --- --- --- --- --- --- --- -
# --- FIT THE DISTRIBUTIONS ---
# --- --- --- --- --- --- --- -
###############################


def MaxwellianDist(x, n, T):  # Fits the NATURAL LOG of a Maxwellian for a Uniform/Thermalized Plasma
    return 4 * pi * n * ((m_e / (2 * pi * q0 * T)) ** 1.5) * (2 * x / m_e) * exp(-x / (T))


def KappaDist(x, kappa, n, T):
    return n * gamma(kappa + 1) * (
                (1 + 2 * q0 * x / (m_e * kappa * ((sqrt(q0 * T * (2 * kappa - 3) / (kappa * m_e))) ** 2))) ** (
                    -kappa - 1)) / (((sqrt(pi) * (sqrt(q0 * T * (2 * kappa - 3) / (kappa * m_e)))) ** 3) * (
                (kappa ** 1.5) * gamma(kappa - 0.5)))




############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################
prgMsg('Plotting Figure')

# DETERMINE THE PLOT DIMESNIONS
fig = plt.figure()
fig.set_size_inches(figure_width, figure_height)
gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[0.5, 0.5], hspace=0.12)

# Top Pllot
gs00 = gs0[0].subgridspec(1, 1)
axDistOverview = fig.add_subplot(gs00[0])

# Subplots
gs01 = gs0[1].subgridspec(1, 1)
ax_Compare = fig.add_subplot(gs01[0,0])
subAxes = [ax_Compare]


# ---Distribution Data Overview---
cmap = axDistOverview.pcolormesh(Epoch_timeSince, Energy, omniDist_high.T, cmap=cmap,vmin=cbarMin,vmax=cbarMax, norm='log')
axDistOverview.set_ylabel('Energy [eV]', fontsize=Label_FontSize)
axDistOverview.tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
axDistOverview.tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize_minor, length=Tick_Length_minor, width=Tick_Width_minor)
axDistOverview.set_yscale('log')
axDistOverview.set_title(r'Distribution Function, $\alpha = 10^{\circ}$')

axDistOverview.axvspan(invertedV_timeSince[0], invertedV_timeSince[1], color='gray', alpha=0.25,zorder=1)
axDistOverview.axvspan(plasmaSheet_timeSince[0], plasmaSheet_timeSince[1], color='gray', alpha=0.25,zorder=1)

# --- INVERTED-V parallel overlap ---
ax_Compare.set_title('Plasma Sheet and Inverted-V',fontsize=Label_FontSize)
ax_Compare.scatter(paraEngy[0], paraDist[0], label='Plasma Sheet Ambient', color='blue')
ax_Compare.scatter(paraEngy[1], paraDist[1], label='Inverted-V', color='red')
ax_Compare.set_yscale('log')
ax_Compare.set_ylim(distLimits[0],distLimits[1])
ax_Compare.set_xlim(xAxis_energyLimits[0],xAxis_energyLimits[1])
ax_Compare.set_xscale('symlog')
ax_Compare.grid()
ax_Compare.set_ylabel(r'Parallel Distribution Function',fontsize=Label_FontSize)
ax_Compare.legend(prop={'size': legend_fontSize})

# --- cbar ---
# cax = fig.add_axes([0.91, 0.288, 0.02, 0.592])
cbar = plt.colorbar(cmap)
# cbar.set_label('Omni-Dir. diff E. Flux \n' + '[cm$^{-2}$str$^{-1}$eV/eV]', rotation=-90, labelpad=20, fontsize=General_LabelFontSize)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbarTickLabelSize + 5)


fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)  # remove the space between plots
plt.tight_layout()
plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot8_DistributionFunc_base.png', dpi=dpi)
Done(start_time)

