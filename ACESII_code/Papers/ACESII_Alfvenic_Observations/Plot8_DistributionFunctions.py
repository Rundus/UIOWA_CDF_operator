# --- Plot8_DistributionFunctions.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Plot the mapped distribution function at various altitudes

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
from ACESII_code.myImports import *
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
print(color.UNDERLINE + f'Plot8_Conjugacy' + color.END)
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from ACESII_code.class_var_func import EpochTo_T0_Rocket, q0,m_e, Re
import matplotlib.gridspec as gridspec
from functools import partial
from ACESII_code.Science.InvertedV.Evans_class_var_funcs import diffNFlux_for_mappedMaxwellian

#########################
# --- Physics Toggles ---
#########################
datasetReduction_TargetTime = [dt.datetime(2022, 11, 20, 17, 24, 50, 000000), dt.datetime(2022,11,20,17,25,15,000000)]
targetVar = [datasetReduction_TargetTime, 'Epoch']
invertedV_TargetTimes_data = dt.datetime(2022, 11, 20, 17, 25, 1, 212000)

# diffNFlux fit parameters
wPitchParameters = 0
fitParameters = [[7.42533896,  86.9711336,  1, 126.73113657], # 0deg
                 [7.34108025, 111.85470989, 1, 147.92852777], # 10deg
                 [2.94001376, 104.81292613, 1, 201.52473787] # 20deg
                 ] # values taken from ''2022, 11, 20, 17, 25, 1, 212000''. Format: density [cm^-3], Temp [eV], beta [doens't matter],Accel Potential [eV]
diffNFlux_Fit_threshEngy = 100
countNoiseLevel = 4

# backscatter toggles
BackScatter_EnergyThreshold = 100 # in eV
BackScatter_wPitchAngles = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

################################
# --- Plot toggles - General ---
################################
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
Plot_MarkerSize = 14
Legend_fontSize = 15
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
diffFlux_NoiseCount = np.zeros(shape=(len(Energy)))
geo_factor = rocketAttrs.geometric_factor[0]
count_interval = 0.8992E-3

for engy in range(len(Energy)):
    deltaT = (count_interval - countNoiseLevel * rocketAttrs.deadtime[0])
    diffFlux_NoiseCount[engy] = countNoiseLevel / (geo_factor[0] * deltaT*Energy[engy])

# Time since launch
LaunchDateTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 000000))
Epoch_timeSince = EpochTo_T0_Rocket(InputEpoch=Epoch, T0=LaunchDateTime)
diffNFlux_ChoiceIdx = np.abs(Epoch - invertedV_TargetTimes_data).argmin()
diffNFluxEpochTimeStamp = Epoch_timeSince[diffNFlux_ChoiceIdx]

Done(start_time)

################################
# --- --- --- --- --- --- --- --
# --- PLOT THE DISTRIBUTIONS ---
# --- --- --- --- --- --- --- --
################################

# --- plot parameters ---
fig = plt.figure()
fig.set_size_inches(figure_width, figure_height)
gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1/4,3/4])
gs00 = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=gs0[1, :])
fig.suptitle(f'Time Since Launch {diffNFluxEpochTimeStamp}', fontsize=Title_FontSize)


# --- diffNFlux fit plot ---
ax_diffNFluxFit = fig.add_subplot(gs0[0, :])

# collect data
fitData = diffNFlux[diffNFlux_ChoiceIdx][wPitchParameters+1]
EngyIdx = np.abs(Energy - diffNFlux_Fit_threshEngy).argmin()
peakDiffNVal = fitData[:EngyIdx].max()
peakDiffNVal_index = np.argmax(fitData[:EngyIdx])
xData_fit = np.array(Energy[:peakDiffNVal_index+1])
yData_fit = np.array(fitData[:peakDiffNVal_index+1])
nonZeroIndicies = np.where(yData_fit!=0)[0]
xData_fit = xData_fit[nonZeroIndicies]
yData_fit = yData_fit[nonZeroIndicies]
fittedX = np.linspace(xData_fit.min(), xData_fit.max(), 100)
fitFuncAtPitch = partial(diffNFlux_for_mappedMaxwellian, alpha=Pitch[wPitchParameters+1])
fittedY = fitFuncAtPitch(fittedX, *fitParameters[wPitchParameters])

# plot it
ax_diffNFluxFit.plot(Energy, fitData, '-o')
ax_diffNFluxFit.plot(fittedX, fittedY, color='red', label=f'n = {fitParameters[wPitchParameters][0]}' +' [cm$^{-3}$]' +
                                                          f'\n T = {fitParameters[wPitchParameters][1]} [eV]\n'
                                                          + f'V = {fitParameters[wPitchParameters][3]} [eV]\n')
ax_diffNFluxFit.axvline(Energy[peakDiffNVal_index],color='red')
ax_diffNFluxFit.set_yscale('log')
ax_diffNFluxFit.set_xscale('log')
ax_diffNFluxFit.set_xlabel('Energy [eV]', fontsize=Label_FontSize)
ax_diffNFluxFit.set_ylabel('diffNFlux [cm$^{-2}$s$^{-1}$str$^{-1}$ eV/eV]', fontsize=Label_FontSize-4)
ax_diffNFluxFit.set_xlim(28, 1E3)
ax_diffNFluxFit.set_ylim(1E4, 1E8)

ax_diffNFluxFit.plot(Energy, diffFlux_NoiseCount, color='black', label=f'{countNoiseLevel}-count level')
ax_diffNFluxFit.legend(fontsize=Legend_fontSize)



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
