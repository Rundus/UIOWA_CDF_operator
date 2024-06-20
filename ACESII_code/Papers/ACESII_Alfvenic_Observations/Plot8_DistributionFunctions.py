# --- Plot8_DistributionFunctions.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Plot the mapped distribution function at various altitudes


# TODO: Handle the betavalue in the backscatter mapping better because beta=1 should correspond to AT the inverted-V altitude
# TODO: NOT at the ionosphere (as it is currently programmed)

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from ACESII_code.myImports import *
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
print(color.UNDERLINE + f'Plot8_DistributionFunctions' + color.END)

from myspaceToolsLib.time import EpochTo_T0_Rocket
import matplotlib.gridspec as gridspec
from functools import partial
from ACESII_code.Science.InvertedV.Evans_class_var_funcs import *
from scipy.interpolate import CubicSpline

#########################
# --- Physics Toggles ---
#########################

# --- Collect the Data ---
datasetReduction_TargetTime = [dt.datetime(2022, 11, 20, 17, 24, 50, 000000), dt.datetime(2022,11,20,17,25,15,000000)]
targetVar = [datasetReduction_TargetTime, 'Epoch']
mappingAltitudes = [6270, 4850, 2252] # These correspond to 2252 km, 4850 km and  6270 km (respectively), which correspond to the STEB TOFs


# --- PRIMARY BEAM diffNFlux fit  ---
primaryBeam_TargetTimes_data = dt.datetime(2022, 11, 20, 17, 25, 1, 212210)
AverageFitData = True # use an average set of data across the SAME pitch angle
NSweeps = 1 # the WIDTH of How many Sweeps to Average Over. E.g. if NSweeps = 2 --> +2 sweeps on either side of the center point
wPitchToPlot = 2
fitParameters = [5.2,  108.9, 185.5] #  average values taken from 0deg bin fitting. Format: density [cm^-3], Temp [eV], beta [doens't matter],Accel Potential [eV]
diffNFlux_Fit_threshEngy = 170
countNoiseLevel = 4


# --- Model Primary Beam toggles ---
N = 101
PS_BeamThreshEnergy = 500 # in eV
PS_BeamThreshPitch = 90
InV_BeamThreshEnergy = 500 # in eV
InV_BeamThreshPitch = 90

# backscatter toggles
BackScatter_TargetTimes_AverageData = [dt.datetime(2022, 11, 20, 17, 25, 1, 312207),dt.datetime(2022, 11, 20, 17, 25, 1, 512208)]
BackScatter_EnergyThreshold = [28, fitParameters[2]] # in eV
BackScatter_wPitchAngles = np.array([10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])

Plot_ModelSlice = False # show the Model Slice
Plot_useDistributionFunction = False # use distribution Function instead of diffNFlux

Plot_Alfven_Resonance_Bands = True # # Model Inertial Alfven Wave Parameters and plot them
wavePhi = 50 # the electric potential


################################
# --- Plot toggles - General ---
################################
figure_width = 13 # in inches
figure_height = 18 # in inches
Title_FontSize = 25
Label_FontSize = 20
Label_Padding = 8
Line_LineWidth = 2.2
Text_Fontsize = 15
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
VelSpace_Max = 2.55

# --- Cbar ---
mycmap = apl_rainbow_black0_cmap()
if Plot_useDistributionFunction:
    cbarMin, cbarMax = 1E-22, 1E-12
else:
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

#####################################
# --- TOP PLOT: diffNFlux fitting ---
#####################################
diffNFlux_ChoiceIdx = np.abs(Epoch - primaryBeam_TargetTimes_data).argmin()
EngyIdx = np.abs(Energy - diffNFlux_Fit_threshEngy).argmin()+1

# --- Perform the fit ---
fitData = diffNFlux[diffNFlux_ChoiceIdx][wPitchToPlot]
peakDiffNVal = fitData[:EngyIdx].max()
peakDiffNVal_index = np.argmax(fitData[:EngyIdx])
xData_fit = np.array(Energy[:peakDiffNVal_index+1])
yData_fit = np.array(fitData[:peakDiffNVal_index+1])
nonZeroIndicies = np.where(yData_fit!=0)[0]
xData_fit = xData_fit[nonZeroIndicies]
yData_fit = yData_fit[nonZeroIndicies]
fittedX = np.linspace(xData_fit.min(), xData_fit.max(), 100)
fitFuncAtPitch = partial(diffNFlux_for_mappedMaxwellian, alpha=Pitch[wPitchToPlot])
fittedY = fitFuncAtPitch(fittedX, *fitParameters)

# --- STD dev ---
dataToAverage = np.transpose(diffNFlux[(diffNFlux_ChoiceIdx-NSweeps):(diffNFlux_ChoiceIdx+NSweeps),wPitchToPlot]) # Collect the data
data_stdDevs = np.array([np.std(engyArr[np.where(engyArr > 0)]) for engyArr in dataToAverage])

# --- ChiSquare ---
stdDev_fit = data_stdDevs[EngyIdx:]# get the error bars only for the data in the fit range
print(stdDev_fit)
ChiSquare = (1/(3-1))*sum([(fitFuncAtPitch(xData_fit[i],*fitParameters) - yData_fit[i])**2 / (stdDev_fit[i]**2) for i in range(len(xData_fit))])

Done(start_time)

##################################################
# --- VELOCITY SPACE GRID VALUES USED IN MODEL ---
##################################################
Emax_vel = np.sqrt(2*q0*2000/m_e)
Vperp_Vals = np.linspace(-1*Emax_vel, Emax_vel, N)
Vpara_Vals = np.linspace(0, Emax_vel, N)
Vperp_Grid, Vpara_Grid = np.meshgrid(Vperp_Vals, Vpara_Vals)
DistFunc_Grid = dist_Maxwellian(Vperp=Vperp_Grid,Vpara=Vpara_Grid,model_Params=fitParameters)


###########################
# --- BackScatter Data ---
###########################
prgMsg('Calculating BackScatter')

# Determine the Average distribution Function across all energies up to BackScatter_EnergyThreshold
thresholdIdx = np.abs(Energy - BackScatter_EnergyThreshold[1]).argmin()
BackScatter_Energies = Energy[thresholdIdx:]

# handle fillvals and Average the lower energies over all pitch angles
tempData = np.array(diffNFlux[diffNFlux_ChoiceIdx, BackScatter_wPitchAngles, thresholdIdx:]).T
backScatterData = []
for engyRow in tempData:
    temp = [val for val in engyRow if np.abs(val) <= 1E11]
    backScatterData.append(sum(temp)/len(temp))

# Sort the data
BackScatter_Energies, backScatterData = zip(*sorted(zip(BackScatter_Energies,backScatterData)))

# interpolate the curve onto the Primary Beam's Energy range
BackScatter_diffNFlux_spline = CubicSpline(np.array(BackScatter_Energies), np.array(backScatterData))
BackScatter_diffNFlux_Grid = calc_BackScatter_onto_Velspace(VperpGrid=Vperp_Grid,
                                                            VparaGrid=Vpara_Grid,
                                                            BackScatterSpline=BackScatter_diffNFlux_spline,
                                                            EngyLimit=BackScatter_EnergyThreshold)
# print(BackScatter_diffNFlux_Grid)
BackScatter_distFunc_Grid = diffNFlux_to_distFunc(Vperp=Vperp_Grid,
                                                  Vpara=Vpara_Grid,
                                                  diffNFlux=BackScatter_diffNFlux_Grid)

# Map the backscatter FROM the ionosphere to the Magnetosphere at various beta values
Data_BackScatter_diffNFluxAtBeta = []
Data_BackScatter_VperpAtBeta = []
Data_BackScatter_VparaAtBeta = []
for mapAlt in mappingAltitudes:
    VperpGrid_mapped, VparaGrid_mapped, diffNFlux_mapped = mapping_VelSpace_magMirror(VperpGrid=Vperp_Grid,
                                                                                       VparaGrid=Vpara_Grid,
                                                                                       distFuncGrid=BackScatter_distFunc_Grid,
                                                                                       targetAlt=400,
                                                                                       startingAlt=mapAlt,
                                                                                       mapToMagSph=True)
    Data_BackScatter_VperpAtBeta.append(VperpGrid_mapped)
    Data_BackScatter_VparaAtBeta.append(VparaGrid_mapped)
    Data_BackScatter_diffNFluxAtBeta.append(diffNFlux_mapped)


Done(start_time)


##########################################
# --- Primary Beam + BackScatter Model ---
##########################################
prgMsg('Calculating Primary Beam')

# --- INVERTED V AT DIFFERENT BETA ---
Data_InvertedV_VperpAtBeta = []
Data_InvertedV_VparaAtBeta = []
Data_InvertedV_diffNFluxAtBeta = []

# Generate the new ACCELERATED distribution
VperpGrid_inV, VparaGrid_inV, distGrid_inV, diffNFluxGrid_inV = calc_velSpace_DistFuncDiffNFluxGrid(Vperp_gridVals=Vperp_Vals,
                                                                                    Vpara_gridVals=Vpara_Vals,
                                                                                    model_Params=fitParameters,
                                                                                    initalBeamParams=[InV_BeamThreshPitch,InV_BeamThreshEnergy])

for mapAlt in mappingAltitudes:

    # Map it to a different beta altitude

    VperpGrid_mapped, VparaGrid_mapped, diffNFlux_mapped = mapping_VelSpace_magMirror(
        VperpGrid=VperpGrid_inV,
        VparaGrid=VparaGrid_inV,
        distFuncGrid=distGrid_inV,
        targetAlt=mapAlt,
        startingAlt=max(mappingAltitudes),
        mapToMagSph=False)
    Data_InvertedV_VperpAtBeta.append(VperpGrid_mapped)
    Data_InvertedV_VparaAtBeta.append(VparaGrid_mapped)
    Data_InvertedV_diffNFluxAtBeta.append(diffNFlux_mapped)



# --- PLASMA SHEET AT DIFFERENT BETAs ---
Data_PS_VperpAtBeta = []
Data_PS_VparaAtBeta = []
Data_PS_diffNFluxAtBeta = []
fitParamsNoAccel = deepcopy(fitParameters)
fitParamsNoAccel[-1] = 0
# Generate the new ACCELERATED distribution
VperpGrid_PS, VparaGrid_PS, distGrid_PS, diffNFluxGrid_PS = calc_velSpace_DistFuncDiffNFluxGrid(
    Vperp_gridVals=Vperp_Vals,
    Vpara_gridVals=Vpara_Vals,
    model_Params=fitParamsNoAccel,
    initalBeamParams=[PS_BeamThreshPitch, PS_BeamThreshEnergy])

for mapAlt in mappingAltitudes:

    # Map it to a different beta altitude
    VperpGrid_mapped, VparaGrid_mapped, diffNFlux_mapped = mapping_VelSpace_magMirror(
        VperpGrid=VperpGrid_PS,
        VparaGrid=VparaGrid_PS,
        distFuncGrid=distGrid_PS,
        targetAlt=mapAlt,
        startingAlt=max(mappingAltitudes),
        mapToMagSph=False)
    Data_PS_VperpAtBeta.append(VperpGrid_mapped)
    Data_PS_VparaAtBeta.append(VparaGrid_mapped)
    Data_PS_diffNFluxAtBeta.append(diffNFlux_mapped)

Done(start_time)

if Plot_ModelSlice:
    ##########################################
    # --- Collect the Slice at V_perp = 0 ---
    ##########################################
    MinVperpIdx = np.abs(Data_BackScatter_VperpAtBeta[-1][0] - 0).argmin()
    MinVperpEnergies = [0.5*m_e*(val**2) for val in Vpara_Vals]
    Data_BackScatter_Energies_Slice = np.array([0.5*m_e*(val**2)/q0 for val in np.transpose(Data_BackScatter_VparaAtBeta[-1])[MinVperpIdx]])
    Data_BackScatter_diffNflux_Slice = np.transpose(Data_BackScatter_diffNFluxAtBeta[-1])[MinVperpIdx]
    Data_InvertedV_Energies_Slice = np.array([0.5*m_e*(val**2)/q0 for val in np.transpose(Data_InvertedV_VparaAtBeta[-1])[MinVperpIdx] ])
    Data_InvertedV_diffNFlux_Slice = np.transpose(Data_InvertedV_diffNFluxAtBeta[-1])[MinVperpIdx]


##################################################
# --- COLLECT INERTIAL ALFVEN SPEED MODEL DATA ---
##################################################
if Plot_Alfven_Resonance_Bands:
    # Use my own model data and my chosen potential to determine some resonance Velocities
    data_dict_alfModel = loadDictFromFile(r'C:\Data\ACESII\science\simulations\TestParticle\plasmaEnvironment\plasmaEnvironment.cdf')

    ResonanceLimits = []
    for mapAlt in mappingAltitudes:
        Altitude = mapAlt*1000
        altIdx = np.abs(data_dict_alfModel['simAlt'][0] - Altitude).argmin()
        localAlfSpdInertial = data_dict_alfModel['alfSpdInertial'][0][altIdx]
        ResonanceLimits.append([localAlfSpdInertial - np.sqrt(2*q0*wavePhi/m_e), localAlfSpdInertial])



################################
# --- --- --- --- --- --- --- --
# --- PLOT THE DISTRIBUTIONS ---
# --- --- --- --- --- --- --- --
################################
prgMsg('Plotting Fig')

# --- plot parameters ---
fig = plt.figure()
fig.set_size_inches(figure_width, figure_height)

if Plot_ModelSlice:
    gs0 = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[1/4, 2/4,1/4])
else:
    gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1/5, 4/5])

gs00 = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=gs0[1, :])

# --- diffNFlux fit plot ---
ax_diffNFluxFit = fig.add_subplot(gs0[0, :])

# plot it
print(data_stdDevs)

ax_diffNFluxFit.errorbar(x=Energy, y=fitData, yerr=data_stdDevs,xerr=None, color='tab:blue', linewidth=Line_LineWidth, label= f'{Pitch[wPitchToPlot]}'+'$^{\circ}$ Pitch Data (T = 301.21 sec)', capsize=7)
ax_diffNFluxFit.plot(Energy, fitData, color='tab:blue', marker='o', markersize=4)
ax_diffNFluxFit.plot(BackScatter_Energies,backScatterData,marker='o',color='tab:purple',label='Avg. BackScatter (Up Going)', linewidth=Line_LineWidth)
ax_diffNFluxFit.plot(fittedX, fittedY, color='tab:red',label='Maxwellian Fit', linewidth=Line_LineWidth)
ax_diffNFluxFit.text(x=550,y=3E5,s=f'  $n_0$ = {round(fitParameters[0],1)}' +' [cm$^{-3}$]' +
                                                          f'\n  T = {round(fitParameters[1],1)} [eV]\n'
                                                          + f'$V_0$ = {round(fitParameters[2],1)} [eV]\n'
                                                          + r'$\chi^{2}_{\nu}$'+f'{ChiSquare}', fontsize=Text_Fontsize, color='tab:red')

ax_diffNFluxFit.text(x=100, y=7E4, s='Secondary/BackScatter Flux', fontsize=Text_Fontsize, ha='right')
ax_diffNFluxFit.text(x=200, y=7E4, s='Primary Beam Flux', fontsize=Text_Fontsize,ha='left')
ax_diffNFluxFit.axvline(Energy[peakDiffNVal_index], color='black',linestyle='--')
ax_diffNFluxFit.set_yscale('log')
ax_diffNFluxFit.set_xscale('log')
ax_diffNFluxFit.set_xlabel('Energy [eV]', fontsize=Label_FontSize)
ax_diffNFluxFit.set_ylabel('diffNFlux \n [cm$^{-2}$s$^{-1}$str$^{-1}$ eV$^{-1}$]', fontsize=Label_FontSize)
ax_diffNFluxFit.set_xlim(28, 1E3)
ax_diffNFluxFit.set_ylim(6E4, 1E8)
ax_diffNFluxFit.tick_params(axis='y', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax_diffNFluxFit.tick_params(axis='y', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length )
ax_diffNFluxFit.tick_params(axis='x', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax_diffNFluxFit.tick_params(axis='x', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax_diffNFluxFit.plot(Energy, diffFlux_NoiseCount, color='black', label=f'{countNoiseLevel}-count level')
ax_diffNFluxFit.legend(fontsize=Legend_fontSize, loc='upper right')


# --- distribution func plot - HIGH ALT ---
VperpData = [Data_PS_VperpAtBeta, Data_InvertedV_VperpAtBeta]
VparaData = [Data_PS_VparaAtBeta, Data_InvertedV_VparaAtBeta]
ZData = [Data_PS_diffNFluxAtBeta, Data_InvertedV_diffNFluxAtBeta] if not Plot_useDistributionFunction else [[DistFunc_Grid for i in range(len(mappingAltitudes))], [DistFunc_Grid for i in range(len(mappingAltitudes))]]


for i in range(3): # Which Beta I'm considering
    for j in range(2): # Plasma Sheet vs InvertedV
        ax = fig.add_subplot(gs00[i, j])

        cmapObj = ax.pcolormesh(VperpData[j][i]/VelSpace_Norm, VparaData[j][i]/VelSpace_Norm, ZData[j][i], vmin=cbarMin, vmax=cbarMax, cmap=mycmap, norm='log')

        # General Labels

        ax.set_ylim(0, VelSpace_Max)
        ax.set_xlim(-VelSpace_Max, VelSpace_Max)
        ax.invert_yaxis()
        ax.tick_params(axis='y', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
        ax.tick_params(axis='y', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length / 2)
        ax.tick_params(axis='x', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
        ax.tick_params(axis='x', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length / 2)
        ax.text(x=-VelSpace_Max+0.2,y=0.1,s=rf'{mappingAltitudes[i]} km', fontsize=Label_FontSize,bbox=dict(facecolor='white', edgecolor='black',boxstyle='round'),va='top',ha='left')

        # Specific things
        if i == 0 and j == 0:
            ax.set_title('Plasma Sheet', fontsize=Title_FontSize)
        if i == 0 and j == 1:
            ax.set_title('Inverted-V',fontsize=Title_FontSize)

        if j==1: # plot the backscatter
            ax.pcolormesh(Data_BackScatter_VperpAtBeta[i] / VelSpace_Norm, Data_BackScatter_VparaAtBeta[i] / VelSpace_Norm, Data_BackScatter_diffNFluxAtBeta[i], vmin=cbarMin, vmax=cbarMax, cmap=mycmap, norm='log')
        if j == 1 and i == 0: # plot the accelerated potential line
            ax.axhline(np.sqrt(2 * fitParameters[2] * q0 / m_e) / VelSpace_Norm, color='tab:red')
            ax.text(x=1, y=(np.sqrt(2 * fitParameters[2] * q0 / m_e) / VelSpace_Norm)-0.075,s='$V_{0}$' + f'= {round(fitParameters[2],1)} eV',fontsize=Text_Fontsize, color='black'  )
        if j == 0: # plot Vparallel labels
            ax.set_ylabel('$V_{\parallel}$ ' + f'[{VelSpace_Norm / 1000} km/s]', fontsize=Label_FontSize)
        if i == 2:# plot Vperp labels
            ax.set_xlabel('$V_{\perp}$ ' + f'[{VelSpace_Norm / 1000} km/s]', fontsize=Label_FontSize)
        if i == 2 and j == 1 and Plot_ModelSlice: # plot the slice of model data
            ax.axvline(Data_InvertedV_VperpAtBeta[-1][0][MinVperpIdx]/VelSpace_Norm,color='black')

        if Plot_Alfven_Resonance_Bands:
            ax.axhline(ResonanceLimits[i][0]/VelSpace_Norm, linestyle='--', color='black',linewidth=Line_LineWidth)
            ax.text(x=0.65, y= sum(ResonanceLimits[i])/(2* VelSpace_Norm)+0.075, s=r'$|\frac{\omega}{k} - v_{\parallel} |< \sqrt{\frac{2e\phi}{m_{e}}}$', fontsize=Text_Fontsize,bbox=dict(facecolor='none', edgecolor='none',boxstyle='round'))
            ax.axhline(ResonanceLimits[i][1]/VelSpace_Norm, color='black',linewidth=Line_LineWidth)
            ax.fill_between([-1*VelSpace_Max,VelSpace_Max],ResonanceLimits[i][0]/VelSpace_Norm,ResonanceLimits[i][1]/VelSpace_Norm,color='grey',alpha=0.15)



if Plot_ModelSlice:
    # Vperp == 0 Slice
    ax_diffNFluxSlice = fig.add_subplot(gs0[2, :])
    ax_diffNFluxSlice.plot(Data_BackScatter_Energies_Slice,Data_BackScatter_diffNflux_Slice,marker='o', color='tab:purple')
    ax_diffNFluxSlice.plot(Data_InvertedV_Energies_Slice,Data_InvertedV_diffNFlux_Slice,marker='o', color='tab:red')
    ax_diffNFluxSlice.set_yscale('log')
    ax_diffNFluxSlice.set_xscale('log')
    ax_diffNFluxSlice.set_xlim(28,1E3)
    ax_diffNFluxSlice.set_ylim(1E4, 1E8)



# Cbar
cax = fig.add_axes([0.90, 0.05, 0.025, 0.64])
cbar = plt.colorbar(cmapObj,cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbarTickLabelSize)
cbar.ax.get_yaxis().labelpad = 17
cbar.set_label(r'DiffNFlux [cm$^{-2}$str$^{-1}$s$^{-1}$eV$^{-1}$]', fontsize=Label_FontSize)
for l in cbar.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(cbarFontSize)

plt.subplots_adjust(left=0.08, bottom=0.05, right=0.89, top=0.98, wspace=None, hspace=None)
try:
    os.remove(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot8_DistributionFunc_Base.png')
except:
    print('')
plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot8_DistributionFunc_Base.png',dpi=dpi)
Done(start_time)