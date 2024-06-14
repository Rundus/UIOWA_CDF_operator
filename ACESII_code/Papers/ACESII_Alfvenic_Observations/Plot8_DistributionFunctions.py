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
from ACESII_code.class_var_func import EpochTo_T0_Rocket, q0,m_e, Re
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

# --- Physics Toggles ---
datasetReduction_TargetTime = [dt.datetime(2022,11, 20, 17, 24, 50, 000000), dt.datetime(2022,11,20,17,25,15,000000)]
targetVar = [datasetReduction_TargetTime, 'Epoch']
invertedV_TargetTimes_data = dt.datetime(2022, 11, 20, 17, 25, 1, 212000)
model_n = 7.42533896
model_T = 86.9711
model_V0 = 126.731-15
countNoiseLevel = 4

# --- Plot toggles - General ---
figure_width = 10 # in inches
figure_height =18 # in inches
Title_FontSize = 20
Label_FontSize = 20
Label_Padding = 8
Tick_FontSize = 12
Tick_Length = 1
Tick_Width = 1
Tick_FontSize_minor = 10
Tick_Length_minor = 1
Tick_Width_minor = 1
Plot_LineWidth = 0.5
plot_MarkerSize = 14
legend_fontSize = 15
dpi = 200

# --- Cbar ---
mycmap = apl_rainbow_black0_cmap()
cbarMin, cbarMax = 1E-18, 1E-14
cbarTickLabelSize = 14


##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
prgMsg('Loading Data')

rocketAttrs, b, c = ACES_mission_dicts()

# EEPAA Distribution Data
inputEEPAA_dist_high = glob('C:\Data\ACESII\L3\DistFunc\high\*ACESII_36359_distFunc_eepaa_fullCal*')[0]
data_dict_dist = loadDictFromFile(inputFilePath=inputEEPAA_dist_high, targetVar=targetVar, wKeys_Reduce=['Distribution_Function', 'Epoch'])
inputEEPAA_counts_high = glob('C:\Data\ACESII\L1\high\*eepaa_fullCal*')[0]
data_dict_counts = loadDictFromFile(inputFilePath=inputEEPAA_counts_high, targetVar=targetVar, wKeys_Reduce=['eepaa', 'Epoch'])
inputEEPAA_diffFlux_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0]
data_dict_diffFlux = loadDictFromFile(inputEEPAA_diffFlux_high, targetVar=targetVar, wKeys_Reduce=['Differential_Energy_Flux','Differential_Number_Flux', 'Epoch'])

# Define the data
counts = deepcopy(data_dict_counts['eepaa'][0])
distFunc = deepcopy(data_dict_dist['Distribution_Function'][0])
diffEFlux = deepcopy(data_dict_diffFlux['Differential_Energy_Flux'][0])
diffNFlux = deepcopy(data_dict_diffFlux['Differential_Number_Flux'][0])
Epoch = deepcopy(data_dict_dist['Epoch'][0])
Energy = deepcopy(data_dict_dist['Energy'][0])
Pitch = deepcopy(data_dict_dist['Pitch_Angle'][0])
Done(start_time)

##############################
# --- --- --- --- --- --- ----
# --- CALCULATED VARIABLES ---
# --- --- --- --- --- --- ----
##############################
prgMsg('Preparing Data')

###### 1-Count THRESHOLD LINE ######
distFunc_NoiseCount = np.zeros(shape=(len(Energy)))
geo_factor = rocketAttrs.geometric_factor[0]
count_interval = 0.8992E-3

for engy in range(len(Energy)):
    deltaT = (count_interval) - (countNoiseLevel * rocketAttrs.deadtime[0])
    fluxVal = (countNoiseLevel) / (geo_factor[0] * deltaT)
    distFunc_NoiseCount[engy] = ((cm_to_m * m_e / (q0 * Energy[engy])) ** 2) * (fluxVal / 2)

# Time since launch
LaunchDateTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 000000))
Epoch_timeSince = EpochTo_T0_Rocket(InputEpoch=Epoch, T0=LaunchDateTime)

Done(start_time)

################################
# --- --- --- --- --- --- --- --
# --- PLOT THE DISTRIBUTIONS ---
# --- --- --- --- --- --- --- --
################################

# --- plot parameters ---
fig = plt.figure()
gs0 = gridspec.GridSpec(2, 1, figure=fig)
gs00 = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=gs0[1, :])
fig.suptitle()


# --- diffNFlux fit plot ---
ax_diffNFluxFit = fig.add_subplot(gs0[0, :])


# --- distribution func plot - HIGH ALT ---
ax_atBeta_PS = fig.add_subplot(gs00[0, 0])
ax_atBeta_invertedV = fig.add_subplot(gs00[0, 1])

# --- distribution func plot - MID ALT ---
ax_atMidAlt_PS = fig.add_subplot(gs00[1, 0])
ax_atMidAlt_invertedV = fig.add_subplot(gs00[1, 1])

# --- distribution func plot - LOW ALT ---
ax_atLowAlt_PS = fig.add_subplot(gs00[2, 0])
ax_atLowAlt_invertedV = fig.add_subplot(gs00[2, 1])


plt.tight_layout()
plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot8_DistributionFunc_Base.png')
