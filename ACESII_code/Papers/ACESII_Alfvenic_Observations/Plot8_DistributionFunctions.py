# --- Plot8_DistributionFunctions.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Plot the mapped distribution function at various altitudes

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import os

import numpy as np

from ACESII_code.myImports import *
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
print(color.UNDERLINE + f'Plot8_DistributionFunctions' + color.END)
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from ACESII_code.class_var_func import EpochTo_T0_Rocket, q0,m_e, Re
import matplotlib.gridspec as gridspec
from functools import partial
from ACESII_code.Science.InvertedV.Evans_class_var_funcs import diffNFlux_for_mappedMaxwellian,dist_Maxwellian, calc_DistributionMapping

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

# Model Primary Beam toggles
N = 1000
PS_BeamThreshEnergy = 500 # in eV
PS_BeamThreshPitch = 45
InV_BeamThreshEnergy = 500 # in eV
InV_BeamThreshPitch = 45


# backscatter toggles
BackScatter_AverageData = False # if == False just use invertdV_TargetTimes_data
BackScatter_betaAltitudes = [1, 3, 6.5]
BackScatter_EnergyThreshold = 100 # in eV
BackScatter_wPitchAngles = np.array([10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])

################################
# --- Plot toggles - General ---
################################
figure_width = 14 # in inches
figure_height = 18 # in inches
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

# --- Velocity Space ---
VelSpace_Norm = 1E7
VelSpace_Max = 2

# --- Cbar ---
mycmap = apl_rainbow_black0_cmap()
cbarMin, cbarMax = 1E3, 5E7
cbarTickLabelSize = 14
cbarFontSize = 15


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

###########################################
# --- --- --- --- --- --- --- --- --- --- -
# --- CALCULATED VARIABLES/COLLECT DATA ---
# --- --- --- --- --- --- --- --- --- --- -
###########################################
prgMsg('Preparing Data')

####################################
###### noise THRESHOLD LINE ######
####################################
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

###########################
# --- diffNFlux fitting ---
###########################
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
Done(start_time)

############################
# --- Primary Beam Model ---
############################

# define the Velocity Space grids
Emax_vel = np.sqrt(2*q0*2000/m_e)
Vperp_Vals = np.linspace(-1*Emax_vel, Emax_vel, N)
Vpara_Vals = np.linspace(0, Emax_vel, N)
Vperp_grid, Vpara_grid = np.meshgrid(Vperp_Vals, Vpara_Vals)
distFunc_Vals = dist_Maxwellian(Vperp=Vperp_grid,
                                Vpara=Vpara_grid,
                                n=fitParameters[wPitchParameters][0],
                                T=fitParameters[wPitchParameters][1])


# determine the mapped distribution function
Vperp_mappedGrids_Accel = []
Vpara_mappedGrids_Accel = []
Vperp_mappedGrids_noAccel = []
Vpara_mappedGrids_noAccel = []
diffNFlux_mappedGrids_Accel = []
diffNFlux_mappedGrids_noAccel = []

for betaVal in BackScatter_betaAltitudes:
    distGrid, VperpGrid, VparaGrid, diffNFluxGrid, VperpGrid_Accel, VparaGrid_Accel, diffNFluxGrid_Accel, VperpGrid_iono, VparaGrid_iono, diffNFluxGrid_iono = calc_DistributionMapping(Vperp_gridVals=Vperp_Vals,
                                                                                                                                                                                   Vpara_gridVals=Vpara_Vals,
                                                                                                                                                                                   model_T=fitParameters[wPitchParameters][1],
                                                                                                                                                                                   model_n=fitParameters[wPitchParameters][0],
                                                                                                                                                                                   model_V0=0,
                                                                                                                                                                                   beta=betaVal,
                                                                                                                                                                                   modifyInitialBeam=True,
                                                                                                                                                                                   beamPitchThreshold=PS_BeamThreshPitch,
                                                                                                                                                                                   beamEnergyThreshold=PS_BeamThreshEnergy)
    Vperp_mappedGrids_noAccel.append(VperpGrid_iono/VelSpace_Norm)
    Vpara_mappedGrids_noAccel.append(VparaGrid_iono/VelSpace_Norm)
    diffNFlux_mappedGrids_noAccel.append(diffNFluxGrid_iono)

    distGrid, VperpGrid, VparaGrid, diffNFluxGrid, VperpGrid_Accel, VparaGrid_Accel, diffNFluxGrid_Accel, VperpGrid_iono, VparaGrid_iono, diffNFluxGrid_iono = calc_DistributionMapping(
        Vperp_gridVals=Vperp_Vals,
        Vpara_gridVals=Vpara_Vals,
        model_T=fitParameters[wPitchParameters][1],
        model_n=fitParameters[wPitchParameters][0],
        model_V0=fitParameters[wPitchParameters][3],
        beta=betaVal,
        modifyInitialBeam=True,
        beamPitchThreshold=InV_BeamThreshPitch,
        beamEnergyThreshold=InV_BeamThreshEnergy)
    Vperp_mappedGrids_Accel.append(VperpGrid_iono/VelSpace_Norm)
    Vpara_mappedGrids_Accel.append(VparaGrid_iono/VelSpace_Norm)
    diffNFlux_mappedGrids_Accel.append(diffNFluxGrid_iono)


###########################
# --- BackScatter Data ---
###########################

# Determine the Average distribution Function across all energies up to BackScatter_EnergyThreshold
thresholdIdx = np.abs(Energy - BackScatter_EnergyThreshold).argmin()
BackScatter_Energies = Energy[thresholdIdx:]
if BackScatter_AverageData:
    lowIdx, highIdx = np.abs(Epoch - datasetReduction_TargetTime[0]).argmin(), np.abs(Epoch - datasetReduction_TargetTime[1]).argmin()
    backScatterRangeData = diffNFlux[lowIdx:highIdx][BackScatter_wPitchAngles][thresholdIdx:]

    # handle fillvals and average
    placeHolder = []
    for tme in range(len(backScatterRangeData)):
        timeSliceData = np.array(backScatterRangeData[tme]).T
        avgDataOverOneTime = []

        for engyRow in timeSliceData:
            temp = [val for val in engyRow if np.abs(val) <= 1E11]
            avgDataOverOneTime.append(sum(temp) / len(temp))

        placeHolder.append(avgDataOverOneTime)

    # Average over all times
    backScatterData = np.sum(placeHolder,axis=0)
else:
    # handle fillvals and Average the lower energies over all pitch angles
    tempData = np.array(distFunc[diffNFlux_ChoiceIdx, BackScatter_wPitchAngles, thresholdIdx:]).T
    backScatterData = []
    for engyRow in tempData:
        temp = [val for val in engyRow if np.abs(val) <= 1E11]
        backScatterData.append(sum(temp)/len(temp))

# --- Use the average BackScatter curves to create a Velocity Space mapping ---

# interpolate the curve onto the Primary Beam's Energy range
from scipy.interpolate import _cubi











################################
# --- --- --- --- --- --- --- --
# --- PLOT THE DISTRIBUTIONS ---
# --- --- --- --- --- --- --- --
################################

# --- plot parameters ---
fig = plt.figure()
fig.set_size_inches(figure_width, figure_height)
gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1/4, 3/4])
gs00 = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=gs0[1, :])
fig.suptitle(f'Time Since Launch {diffNFluxEpochTimeStamp}', fontsize=Title_FontSize)

# --- diffNFlux fit plot ---
ax_diffNFluxFit = fig.add_subplot(gs0[0, :])

# plot it
ax_diffNFluxFit.plot(Energy, fitData, '-o')
ax_diffNFluxFit.plot(fittedX, fittedY, color='red', label=f'n = {fitParameters[wPitchParameters][0]}' +' [cm$^{-3}$]' +
                                                          f'\n T = {fitParameters[wPitchParameters][1]} [eV]\n'
                                                          + f'V = {fitParameters[wPitchParameters][3]} [eV]\n')
ax_diffNFluxFit.axvline(Energy[peakDiffNVal_index], color='red')
ax_diffNFluxFit.set_yscale('log')
ax_diffNFluxFit.set_xscale('log')
ax_diffNFluxFit.set_xlabel('Energy [eV]', fontsize=Label_FontSize)
ax_diffNFluxFit.set_ylabel('diffNFlux \n [cm$^{-2}$s$^{-1}$str$^{-1}$ eV$^{-1}$]', fontsize=Label_FontSize)
ax_diffNFluxFit.set_xlim(28, 1E3)
ax_diffNFluxFit.set_ylim(1E4, 1E8)
ax_diffNFluxFit.tick_params(axis='y', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax_diffNFluxFit.tick_params(axis='y', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length )
ax_diffNFluxFit.tick_params(axis='x', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax_diffNFluxFit.tick_params(axis='x', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax_diffNFluxFit.plot(Energy, diffFlux_NoiseCount, color='black', label=f'{countNoiseLevel}-count level')
ax_diffNFluxFit.plot(BackScatter_Energies, backScatterData,marker='o')
ax_diffNFluxFit.legend(fontsize=Legend_fontSize)



# --- distribution func plot - HIGH ALT ---

# PS at beta
IDX = 0
ax_atBeta_PS = fig.add_subplot(gs00[0, 0])
ax_atBeta_PS.pcolormesh(Vperp_mappedGrids_noAccel[IDX], Vpara_mappedGrids_noAccel[IDX], diffNFlux_mappedGrids_noAccel[IDX],vmin=cbarMin,vmax=cbarMax,cmap=mycmap, norm='log')
ax_atBeta_PS.set_title(rf'$\beta$ = {BackScatter_betaAltitudes[IDX]}',fontsize=Label_FontSize)
ax_atBeta_PS.set_ylabel('$V_{\parallel}$ '+ f'[{VelSpace_Norm/1000} km/s]', fontsize=Label_FontSize)

# inverted-V at beta
ax_atBeta_invertedV = fig.add_subplot(gs00[0, 1])
ax_atBeta_invertedV.pcolormesh(Vperp_mappedGrids_Accel[IDX], Vpara_mappedGrids_Accel[IDX], diffNFlux_mappedGrids_Accel[IDX], vmin=cbarMin,vmax=cbarMax,cmap=mycmap, norm='log')
ax_atBeta_invertedV.axhline(np.sqrt(2*fitParameters[wPitchParameters][3]*q0/m_e)/(VelSpace_Norm),color='red', label='$V_{0}$'+f'= {fitParameters[wPitchParameters][3]} eV')
ax_atBeta_invertedV.set_title(rf'$\beta$ = {BackScatter_betaAltitudes[IDX]}',fontsize=Label_FontSize)
ax_atBeta_invertedV.legend(fontsize=Legend_fontSize)


# --- distribution func plot - MID ALT ---

# PS at midAltbeta
IDX = 1
ax_atMidAlt_PS = fig.add_subplot(gs00[1, 0])
ax_atMidAlt_PS.pcolormesh(Vperp_mappedGrids_noAccel[IDX], Vpara_mappedGrids_noAccel[IDX], diffNFlux_mappedGrids_noAccel[IDX], vmin=cbarMin,vmax=cbarMax,cmap=mycmap, norm='log')
ax_atMidAlt_PS.set_title(rf'$\beta$ = {BackScatter_betaAltitudes[IDX]}',fontsize=Label_FontSize)
ax_atMidAlt_PS.set_ylabel('$V_{\parallel}$ '+ f'[{VelSpace_Norm/1000} km/s]', fontsize=Label_FontSize)

# inverted-V at midAltbeta
ax_atMidAlt_invertedV = fig.add_subplot(gs00[1, 1])
ax_atMidAlt_invertedV.pcolormesh(Vperp_mappedGrids_Accel[IDX], Vpara_mappedGrids_Accel[IDX], diffNFlux_mappedGrids_Accel[IDX], vmin=cbarMin,vmax=cbarMax,cmap=mycmap, norm='log')
ax_atMidAlt_invertedV.set_title(rf'$\beta$ = {BackScatter_betaAltitudes[IDX]}',fontsize=Label_FontSize)


# --- distribution func plot - LOW ALT ---
IDX = 2
# PS at lowAltbeta
ax_atLowAlt_PS = fig.add_subplot(gs00[2, 0])
ax_atLowAlt_PS.pcolormesh(Vperp_mappedGrids_noAccel[IDX], Vpara_mappedGrids_noAccel[IDX], diffNFlux_mappedGrids_noAccel[IDX], vmin=cbarMin,vmax=cbarMax,cmap=mycmap, norm='log')
ax_atLowAlt_PS.set_title(rf'$\beta$ = {BackScatter_betaAltitudes[IDX]}',fontsize=Label_FontSize)
ax_atLowAlt_PS.set_ylabel('$V_{\parallel}$ '+ f'[{VelSpace_Norm/1000} km/s]', fontsize=Label_FontSize)
ax_atLowAlt_PS.set_xlabel('$V_{\perp}$ '+ f'[{VelSpace_Norm/1000} km/s]', fontsize=Label_FontSize)

# invertedV at lowAltbeta
ax_atLowAlt_invertedV = fig.add_subplot(gs00[2, 1])
cmapObj = ax_atLowAlt_invertedV.pcolormesh(Vperp_mappedGrids_Accel[IDX], Vpara_mappedGrids_Accel[IDX], diffNFlux_mappedGrids_Accel[IDX], vmin=cbarMin,vmax=cbarMax,cmap=mycmap, norm='log')
ax_atLowAlt_invertedV.set_title(rf'$\beta$ = {BackScatter_betaAltitudes[IDX]}',fontsize=Label_FontSize)
ax_atLowAlt_invertedV.set_xlabel('$V_{\perp}$ '+ f'[{VelSpace_Norm/1000} km/s]', fontsize=Label_FontSize)

# adjust the mapped distribution function labels, xlimits
mappedAxes = [ax_atBeta_PS, ax_atBeta_invertedV, ax_atMidAlt_PS, ax_atMidAlt_invertedV, ax_atLowAlt_PS, ax_atLowAlt_invertedV]
for ax in mappedAxes:
    ax.set_ylim(0, VelSpace_Max)
    ax.set_xlim(-VelSpace_Max, VelSpace_Max)
    ax.invert_yaxis()


# Cbar
cax = fig.add_axes([0.89, 0.05, 0.025, 0.64])
cbar = plt.colorbar(cmapObj,cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbarTickLabelSize)
cbar.ax.get_yaxis().labelpad = 17
cbar.set_label(r'diffNFlux [cm$^{-2}$str$^{-1}$s$^{-1}$eV$^{-1}$]', fontsize=Label_FontSize)
for l in cbar.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(cbarFontSize)

plt.subplots_adjust(left=0.08, bottom=0.05, right=1, top=0.99, wspace=None, hspace=None)
try:
    os.remove(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot8_DistributionFunc_Base.png')
except:
    print('')
plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot8_DistributionFunc_Base.png')
