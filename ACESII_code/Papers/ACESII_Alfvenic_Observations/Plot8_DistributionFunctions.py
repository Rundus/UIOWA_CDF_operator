# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing field-aligned
# particle data along with electric and magnetic signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import math

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from ACESII_code.class_var_func import EpochTo_T0_Rocket, q0,m_e
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
disR_TargetTime = [dt.datetime(2022,11, 20, 17, 24, 22, 000000), dt.datetime(2022,11,20,17,25,3,000000)]
targetVar = [disR_TargetTime, 'Epoch']
plasmaSheet_targetTimes = [dt.datetime(2022, 11, 20, 17, 24, 24, 000000), dt.datetime(2022,11,20,17,24,30,000000)]
# plasmaSheet_targetTimes = [dt.datetime(2022, 11, 20, 17, 24, 40, 000000), dt.datetime(2022,11,20,17,24,44,000000)]
invertedV_targetTimes = [dt.datetime(2022,11, 20, 17, 25, 1, 396000), dt.datetime(2022,11,20,17,25,1,596000)]
ptchRange = [2] # determines which pitchs to include in the parallel distribution functions
countNoiseLevel = 1

# --- Cbar ---
cbarMin, cbarMax = 1E-18, 1E-14
cbarTickLabelSize = 14
cmap = apl_rainbow_black0_cmap()
cmap.set_bad(color=(1, 1, 1))

# --- Plot toggles ---
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
distLimits = [1E-19, 1E-13]
xAxis_energyLimits = [10, 1E4]



##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
prgMsg('Loading Data')

rocketAttrs, b, c = ACES_mission_dicts()

# EEPAA Distribution Data
inputEEPAA_high = glob('C:\Data\ACESII\L3\DistFunc\high\*ACESII_36359_distFunc_eepaa_fullCal*')[0]
data_dict_dist_high = loadDictFromFile(inputFilePath=inputEEPAA_high, targetVar=targetVar, wKeys_Reduce=['Distribution_Function', 'Epoch'])

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

# --- Determine the 1-count threshold line ---
distFunc_NoiseCount = np.zeros(shape=(len(Energy)))
geo_factor = rocketAttrs.geometric_factor[0]
count_interval = 0.8992E-3

for engy in range(len(Energy)):
    deltaT = (count_interval) - (countNoiseLevel * rocketAttrs.deadtime[0])
    fluxVal = (countNoiseLevel) / (geo_factor[0] * deltaT)
    distFunc_NoiseCount[engy] = ((cm_to_m*m_e / (q0*Energy[engy]))**2) * (fluxVal/2)





# --- averaged distribution function over pitch angle ---
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
paraDist = [[], []]
paraEngy = [[], []]
regionTimes = [plasmaSheet_targetTimes,invertedV_targetTimes]

for regionIdx in range(len(paraDist)):
    rTime = regionTimes[regionIdx]
    low_idx, high_idx = np.abs(Epoch - rTime[0]).argmin(), np.abs(Epoch - rTime[1]).argmin()
    tempDist = []
    tempEngy = []
    tempPtch = []

    for ptch in ptchRange:
        data = distFunc[low_idx:high_idx, ptch, :].T
        data[np.where(data == rocketAttrs.epoch_fillVal)] = 0

        average = np.array([np.average(arr[np.nonzero(arr)])    for arr in data]) # average over time for each energy bin
        average = np.array([val if np.isnan(val) == False else 0 for val in average ]) # remove the points where there was no distribution function for a specific enregy

        # remove all the zero values from average and from energy
        zeroFinder = np.where(average!=0)


        tempDist.append(average[zeroFinder])
        tempEngy.append(Energy[zeroFinder])
        tempPtch.append(ptch)

    if len(tempPtch) > 1:
        print('GOTTA WRITE SOMETHING')

    else:
        paraDist[regionIdx] = tempDist[0]
        paraEngy[regionIdx] = tempEngy[0]

Done(start_time)


# for i in range(len(paraDist[0])):
#     print(i,paraEngy[0][i],paraDist[0][i])
#     # print(i, Energy[i],distFunc_NoiseCount[i], paraDist[0][i],paraDist[0][i] - distFunc_NoiseCount[i])







###############################
# --- --- --- --- --- --- --- -
# --- FIT THE DISTRIBUTIONS ---
# --- --- --- --- --- --- --- -
###############################
from numpy import pi, exp, sqrt, power
from math import gamma

def MaxwellianDist(x, n, T):  # Fits the NATURAL LOG of a Maxwellian for a Uniform/Thermalized Plasma
    return 4 * pi * n * ((m_e / (2 * pi * q0 * T)) ** 1.5) * (2 * x / m_e) * exp(-x / (T))


def KappaDist(x, kappa, n, T):
    w = np.sqrt( (2*q0*T/m_e) * ((kappa - 3/2)/kappa))
    # term1 = n / (power(pi, 3/2)*power(w,3))
    # term2 = gamma(kappa+1)/(power(kappa,3/2)*gamma(kappa-1/2))
    # term3 = (1 + power(x,2)/(kappa*power(w,2)) )**(-kappa-1)

    term1 = n * gamma(kappa + 1)
    term2 = (1/(((sqrt(pi) * (sqrt(q0 * T * (2 * kappa - 3) / (kappa * m_e)))) ** 3) * ((kappa ** 3/2) * gamma(kappa - 0.5))))
    term3 = ((1 + 2 * q0 * x / (m_e * kappa * ((sqrt(q0 * T * (2 * kappa - 3) / (kappa * m_e))) ** 2))) ** (-kappa - 1))
    return term1*term2*term3


# guess = [1.509 1.3E6, 500]
# NOTE: 1 cm^3 = 1E6 m^3. Magnetospheric particles from plasma sheet are around 1 cm^-3 at 500eV
guess = [1.51, 1E6, 200] # observed plasma at dispersive region is 0.5E5 cm^-3 BUT this doesn't make sense to use as the kappa fit since the kappa fit comes from MUCH less dense populations above
boundVals = [[1.5000000001, 2], [1E4, 1E14], [0.001, 1E4]]  # kappa, n, T
bounds = tuple([[boundVals[i][0] for i in range(len(boundVals))], [boundVals[i][1] for i in range(len(boundVals))]])
# params, cov = curve_fit(KappaDist, paraEngy[0], paraDist[0], p0=guess, maxfev=int(1E9), bounds=bounds)
# xData_fit, yData_fit = np.array(paraEngy[0]), np.array([KappaDist(val, *params) for val in paraEngy[0]])



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
gs01 = gs0[1].subgridspec(2, 1)
ax_Compare = fig.add_subplot(gs01[0,0])
ax_Fit = fig.add_subplot(gs01[1,0])
subAxes = [ax_Compare]


# ---Distribution Data Overview---
cmap = axDistOverview.pcolormesh(Epoch_timeSince, Energy, omniDist_high.T, cmap=cmap,vmin=cbarMin,vmax=cbarMax, norm='log')
axDistOverview.set_ylabel('Energy [eV]', fontsize=Label_FontSize)
axDistOverview.tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
axDistOverview.tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize_minor, length=Tick_Length_minor, width=Tick_Width_minor)
axDistOverview.set_yscale('log')

axDistOverview.axvspan(invertedV_timeSince[0], invertedV_timeSince[1], color='gray', alpha=0.25,zorder=1)
axDistOverview.axvspan(plasmaSheet_timeSince[0], plasmaSheet_timeSince[1], color='gray', alpha=0.25,zorder=1)

# --- INVERTED-V parallel overlap ---
# ax_Compare.set_title('Plasma Sheet and Inverted-V', fontsize=Label_FontSize)
ax_Compare.scatter(paraEngy[0], paraDist[0], label='Plasma Sheet Ambient', color='blue')
ax_Compare.scatter(paraEngy[1], paraDist[1], label='Inverted-V', color='red')
ax_Compare.plot(Energy,distFunc_NoiseCount,label=f'{countNoiseLevel}-count level',color='black')

# ax_Compare.plot(xData_fit,yData_fit,label='Kappa Fit',color='green',alpha=0.5)
ax_Compare.set_yscale('log')
ax_Compare.set_ylim(distLimits[0], distLimits[1])
ax_Compare.set_xlim(xAxis_energyLimits[0], xAxis_energyLimits[1])
ax_Compare.set_xscale('symlog')
ax_Compare.grid()
ax_Compare.set_ylabel(r'Parallel Distribution Function',fontsize=Label_FontSize)
ax_Compare.legend(prop={'size': legend_fontSize})

# --- Fit a Kappa Distribution ---
# ax_Fit.scatter(paraEngy[0], paraDist[0], label='Plasma Sheet Ambient', color='blue')
# legendLabel = '$\kappa$=' + f'{params[0]}\n' + f'n={params[1]}\n' + f'T={params[2]}'
# ax_Fit.plot(xData_fit,yData_fit,label=legendLabel,color='green',alpha=1)
ax_Fit.set_yscale('log')
ax_Fit.set_xscale('symlog')
ax_Fit.set_xlim(xAxis_energyLimits[0],xAxis_energyLimits[1])
ax_Fit.grid()
ax_Fit.legend(prop={'size': legend_fontSize})
ax_Fit.axis('off')

# --- cbar ---
# cax = fig.add_axes([0.91, 0.288, 0.02, 0.592])
# cbar = plt.colorbar(cmap, cax= cax)
# cbar.set_label('Omni-Dir. diff E. Flux \n' + '[cm$^{-2}$str$^{-1}$eV/eV]', rotation=-90, labelpad=20, fontsize=General_LabelFontSize)
# cbar.ax.minorticks_on()
# cbar.ax.tick_params(labelsize=cbarTickLabelSize + 5)

plt.tight_layout

fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)  # remove the space between plots
plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot8_DistributionFunc_base.png', dpi=dpi)
Done(start_time)

