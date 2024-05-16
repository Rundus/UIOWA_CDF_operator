# --- Plots4_pitchAnglePlots.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Recreation of DiffEFlux pitch angle Plots, focusing on
# a few particle signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from ACESII_code.myImports import *
import matplotlib.gridspec as gridspec
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import math
import matplotlib as mpl
from ACESII_code.class_var_func import EpochTo_T0_Rocket

print(color.UNDERLINE + f'Plot7_keyObservations' + color.END)

#################
# --- TOGGLES ---
#################
figure_height = (15)
figure_width = (11)

cmap = 'turbo'
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
cmap = apl_rainbow_black0_cmap()


plot_LineWidth = 3
plot_textFontSize = 20
plot_MarkerSize = 20
plot_Colors = ['tab:blue', 'tab:red', 'tab:orange', 'tab:green', 'tab:purple', 'tab:olive', 'black']
title_FontSize = 14
labels_FontSize = 20
labels_Padding = 1
legend_FontSize = 15
tick_LabelSize = 15
tick_SubplotLabelSize = 10
tick_Width = 2
tick_Length = 4
cbar_FontSize = 15
dpi = 1200

# dispersiveRegionTargetTime = [dt.datetime(2022,11,20,17,24,57,600000),
#                               dt.datetime(2022,11,20,17,25,2,000000)]
dispersiveRegionTargetTime = [dt.datetime(2022,11,20,17,24,55,900000),
                              dt.datetime(2022,11,20,17,25,2,000000)]

# plot toggles - Show STEB itself ----------
cbarLow_counts, cbarHigh_counts = 1, 100
cbarLow_diffEFlux, cbarHigh_diffEFlux = 5E6, 1E9
wPitch_Engy_vs_Time = 2 # the pitch angle index to plot for the Energy vs time plot
Energy_yLimit = 1350

# plot toggles - Histogram ---------------------------
Histogram_countsthresh = 2
# consider only the pitch angles between -10 and 90
# [  0 1  2  3  4  5  6  7  8  9 10 ...]
# [-10 0 10 20 30 40 50 60 70 80 90 ...]
pitchAngleWidthAcceptance_lowerlimit = 2
pitchAngleWidthAcceptance_upperlimit = 10 +1 # NEEDS +1 for array indexing
tempCmap = plt.cm.turbo_r # define new colormap
cmaplist = [tempCmap(i) for i in range(tempCmap.N)] # extract all colors from the colormap
cmap_hist = mpl.colors.LinearSegmentedColormap.from_list('turboCustom',cmaplist,tempCmap.N) # create the new map
bounds = 10*np.array([-1.5+i for i in range(11)])
histNorm = mpl.colors.BoundaryNorm(bounds,tempCmap.N)

# --- plot Poynting Flux Toggles ---
wSTEBtoPlot = [1,2, 3, 4, 5] # STEB number NOT index. Don't -1
PoyntingScale = 1000# convert from W/m^2 to ergs/cm^2

# --- plot COUNTS and ENERGY  toggles ---
wPitchs_to_plot = [2, 3, 4, 5] # decide which pitch angles to get the peak energy for
# countsMask = 3
fluxMask = 5E7
Energy_yLimit_Idx = 15 # ONLY consider energies above this index -> energies BELOW ~ 1345 eV



##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################

prgMsg('Loading Data')
targetVarName = 'Epoch'
targetVar = dispersiveRegionTargetTime

# Attitude data (for geomagnetic lat/long info)
inputFiles_Attitude = [glob(r'C:\Data\ACESII\attitude\high\*Attitude_Solution*')[0], glob(r'C:\Data\ACESII\attitude\low\*Attitude_Solution*')[0]]
data_dict_attitude_high = loadDictFromFile(inputFilePath=inputFiles_Attitude[0])
data_dict_attitude_low = loadDictFromFile(inputFilePath=inputFiles_Attitude[1])

# Magnetometer Data
inputFile_B = 'C:\Data\ACESII\L2\high\ACESII_36359_l2_RingCore_ENU.cdf'  # get the B data
data_dict_B = loadDictFromFile(inputFile_B,targetVar=[targetVar, targetVarName])

inputFile_deltaB = glob('C:\Data\ACESII\L3\deltaB\high\ACESII_36359_RingCore_Field_Aligned_WL250_stitchedFlight.cdf')[0] # get the deltaB data
data_dict_deltaB = deepcopy(loadDictFromFile(inputFile_deltaB,targetVar=[targetVar, targetVarName]))

# Langmuir Data
inputFile_Langmuir = 'C:\Data\ACESII\L3\Langmuir\high\ACESII_36359_langmuir_fixed_LowPass_low0.3_high0.3.cdf'
data_dict_langmuir = deepcopy(loadDictFromFile(inputFile_Langmuir,targetVar=[targetVar, targetVarName],wKeys_Load=['ni', 'Epoch', 'ILat']))
indexVals = [np.abs(data_dict_langmuir['Epoch'][0] - tme).argmin() for i,tme in enumerate(data_dict_B['Epoch'][0])]
data_dict_langmuir['ni'][0] = deepcopy(data_dict_langmuir['ni'][0][indexVals])
data_dict_langmuir['Epoch'][0] = deepcopy(data_dict_langmuir['Epoch'][0][indexVals])

# EISCAT Data (# Up-sample the EISCAT data (it's fine since the EISCAT variables are very slowly varying))
inputFile_EISCAT = 'C:\Data\ACESII\science\EISCAT_ACESII_Slice\high\ACESII_36359_EISCAT_Tromso_rktSlice.cdf'
data_dict_EISCAT = deepcopy(loadDictFromFile(inputFile_EISCAT, targetVar=[targetVar, targetVarName],wKeys_Load=['Ion_Comp', 'Op_Comp','ILat','Epoch']))
data_dict_EISCAT_interp = InterpolateDataDict(InputDataDict=data_dict_EISCAT,InputEpochArray=data_dict_EISCAT['Epoch'][0],targetEpochArray=data_dict_B['Epoch'][0],wKeys=[])
data_dict_EISCAT = deepcopy(data_dict_EISCAT_interp)

# EEPAA Particle Data
inputFiles_eepaa = [glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0], glob('C:\Data\ACESII\L2\low\*eepaa_fullCal*')[0]]
data_dict_eepaa_high = loadDictFromFile(inputFilePath=inputFiles_eepaa[0], targetVar=[targetVar, targetVarName])

inputFiles_eepaa_counts = [glob('C:\Data\ACESII\L1\high\*eepaa_fullCal*')[0], glob('C:\Data\ACESII\L1\low\*eepaa_fullCal*')[0]]
data_dict_counts_high = loadDictFromFile(inputFilePath=inputFiles_eepaa_counts[0], targetVar=[targetVar, targetVarName],wKeys_Reduce=['eepaa','Epoch'])
Done(start_time)


##########################
# --- --- --- --- --- ---
# --- PREPARE THE DATA ---
# --- --- --- --- --- ---
##########################
# --- time ---
from ACESII_code.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes
STEBtimes = [dispersionAttributes.keyDispersionDeltaT[idx-1] for idx in wSTEBtoPlot]
LaunchDateTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022,11,20,17,20,00,000000))
STEBtimes_rkt = [[(pycdf.lib.datetime_to_tt2000(tme[0]) - LaunchDateTime)/1E9, (pycdf.lib.datetime_to_tt2000(tme[1]) - LaunchDateTime)/1E9] for tme in STEBtimes]

rktTime_deltaB = EpochTo_T0_Rocket(InputEpoch=data_dict_deltaB['Epoch'][0], T0=LaunchDateTime)
rktTime_counts = EpochTo_T0_Rocket(InputEpoch=data_dict_counts_high['Epoch'][0], T0=LaunchDateTime)

# --- particles ---
Pitch = data_dict_counts_high['Pitch_Angle'][0]
Energy = data_dict_counts_high['Energy'][0]
counts = data_dict_counts_high['eepaa'][0]
diffEFlux = data_dict_eepaa_high['Differential_Energy_Flux'][0]

# --- Flux ---
B0 = 1E-9 * data_dict_B['Bmag'][0]
dB_e = data_dict_deltaB['B_e'][0] # in nanotesla
dB_r = data_dict_deltaB['B_r'][0]
ni = (cm_to_m ** 3) * data_dict_langmuir['ni'][0]
rhoCalc = (ni*data_dict_EISCAT['Ion_Comp'][0]*((IonMasses[4] + IonMasses[5] + IonMasses[7])/3) + IonMasses[1]*ni*data_dict_EISCAT['Op_Comp'][0])
VA_t = B0 / np.sqrt(u0 * rhoCalc)
VA_avg = sum(VA_t)/len(VA_t)
dBperp = np.array([np.sqrt((dB_e[i]*1E-9)**2 + (dB_r[i]*1E-9)**2) for i in range(len(dB_e))])
S_est = VA_t*(dBperp**2)/ (u0) # calculate Estimated Poynting Flux


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

prgMsg('Beginning Plot')
fig = plt.figure()
fig.set_size_inches(figure_width, figure_height)
gs0 = gridspec.GridSpec(2,1,figure=fig, height_ratios=[2/3,1/3],hspace=0.1)

gs00 = gs0[0].subgridspec(3, 1)
axPeakE = fig.add_subplot(gs00[2])
axPitchHist = fig.add_subplot(gs00[0],sharex=axPeakE)
axPoynting = fig.add_subplot(gs00[1],sharex=axPeakE)

mainThreeAxes = [axPitchHist, axPoynting, axPeakE]

gs01 = gs0[1].subgridspec(2, 2)
axS1 = fig.add_subplot(gs01[0,0])
axS2 = fig.add_subplot(gs01[0,1])
axS3 = fig.add_subplot(gs01[1,0])
axS4 = fig.add_subplot(gs01[1,1])
subAxes = [axS1,axS2,axS3,axS4]

# --- --- --- --- --- --- --- ---
# --- MEDIAN PITCH ANGLE PLOT ---
# --- --- --- --- --- --- --- ---
X, Y = np.meshgrid(rktTime_counts, Energy)
Z = deepcopy(X)

# populate the new data
for tme in range(len(rktTime_counts)):
    for engy in range(len(Energy)):

        # get the data across all pitch angles here
        # pitchData = dataArray_counts[tme, :, engy]
        pitchData = counts[tme, :, engy]

        # Find errors in data and eliminate them from the median calculation
        for k, val in enumerate(pitchData):
            if val < Histogram_countsthresh: # theshold away small COUNT values
                pitchData[k] = 0
            elif val > 1E14: # for extremely large outlier values
                pitchData[k] = 0
            elif math.isnan(val): # if is a nan
                pitchData[k] = 0

        pitchData_useThis = pitchData[pitchAngleWidthAcceptance_lowerlimit:pitchAngleWidthAcceptance_upperlimit]
        medianVal_pitchVal = np.nan

        halfDataVal = sum(pitchData_useThis) / 2
        for h in range(len(pitchData_useThis)):
            if halfDataVal - sum(pitchData_useThis[:h + 1]) < 0:
                medianVal_pitchVal = Pitch[h+pitchAngleWidthAcceptance_lowerlimit]
                break

        Z[engy][tme] = medianVal_pitchVal

# adjust the plot
cmapHist = axPitchHist.pcolormesh(X, Y, Z, cmap=cmap_hist, norm=histNorm)
axPitchHist.set_yscale('log')
axPitchHist.set_ylim(Energy[-1], Energy_yLimit)
axPitchHist.set_ylabel('Energy [eV]', fontsize=labels_FontSize)
axPitchHist.tick_params(axis='x', which='major', labelbottom=False)

# --- --- --- --- --- --- --
# --- FLUX ESTIMATE PLOT ---
# --- --- --- --- --- --- --
axPoynting.plot(rktTime_deltaB, PoyntingScale*S_est, plot_Colors[2],linewidth=plot_LineWidth,zorder=2)
axPoynting.set_ylabel('$\delta$ S$_{p}$ [ergs/cm$^{2}$s]',fontsize=labels_FontSize)
axPoynting.set_ylim(-0.1E-2, 5.3E-2)
axPoynting.set_xmargin(0)
axPoynting.tick_params(axis='both', which='major', labelsize=tick_LabelSize)
axPoynting.tick_params(axis='both', which='minor', labelsize=tick_LabelSize-2)
axPoynting.tick_params(axis='x', which='major', labelbottom=False)

# plot the color regions
for k,tme in enumerate(STEBtimes_rkt):
    lowIdx = np.abs(rktTime_deltaB - tme[0]).argmin()
    highIdx = np.abs(rktTime_deltaB - tme[1]).argmin()
    axPoynting.axvspan(rktTime_deltaB[lowIdx], rktTime_deltaB[highIdx], color=plot_Colors[k], alpha=0.25,zorder=1)
    textPlacement = (STEBtimes_rkt[k][0]+STEBtimes_rkt[k][1])/2
    axPoynting.text(textPlacement, 0.048, f'S{wSTEBtoPlot[k]}', ha='center', weight='bold',fontsize=plot_textFontSize)

# axPoynting.set_xticklabels([])



# --- --- --- --- --- --- --- --- ---
# --- PEAK COUNTS AND ENERGY PLOTS ---
# --- --- --- --- --- --- --- --- ---

outside_peakFlux = []
outside_peakEnergy_atpeakFlux = []
outside_avgtimes = []

for idx, ptchVal in enumerate(wPitchs_to_plot):

    peakFlux = []
    peakEnergy = []
    peakEnergy_atPeakFlux = []
    avgTimes = []

    for t,tme in enumerate(STEBtimes_rkt):
        # determine the x-value for the point (i.e. the average time)
        lowIdx = np.abs(rktTime_counts - tme[0]).argmin()
        highIdx = np.abs(rktTime_counts - tme[1]).argmin()
        avgTimes.append((rktTime_counts[lowIdx] + rktTime_counts[highIdx])/2)

        # Isolate the dispersion Feature
        wDis = wSTEBtoPlot[t]
        wDispersion_key = f's{wDis}'

        lowCut, highCut = np.abs(data_dict_eepaa_high['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis - 1][0]).argmin(), np.abs(data_dict_eepaa_high['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis - 1][1]).argmin()
        Epoch_dis = deepcopy(dateTimetoTT2000(data_dict_eepaa_high['Epoch'][0][lowCut:highCut + 1], inverse=False))
        eepaa_dis_pre = deepcopy(diffEFlux[lowCut:highCut + 1])
        eepaa_dis = deepcopy(np.array(dispersionAttributes.isolationFunctions[wDispersion_key](eepaa_dis_pre, Energy, Epoch_dis)))  # pply the isolation functions found in dispersionAttributes.py
        # prepare the data by removing fillvals and applying the mask
        STEBdata = deepcopy(eepaa_dis)
        STEBdata[STEBdata < 0] = 0  # set anything below 0 = 0
        STEBdata -= fluxMask
        STEBdata[STEBdata < 0] = 0  # set anything below 0 = 0

        # --- PEAK ENERGY ---
        peakEval = 0

        for timeIdx in range(len(STEBdata)):

            engySubSlice = STEBdata[timeIdx, ptchVal, Energy_yLimit_Idx:]
            peakE_Idx = next((i for i, x in enumerate(engySubSlice) if x), None)

            if peakE_Idx == None:
                EpVal = 0
            else:
                EpVal = Energy[Energy_yLimit_Idx+peakE_Idx]

            if EpVal > peakEval:
                peakEval = Energy[Energy_yLimit_Idx+peakE_Idx]

        peakEnergy.append(peakEval)

        # --- PEAK FLUX ---
        # get the peak Counts and Peak Energy
        pitchSlice = STEBdata[:, ptchVal, Energy_yLimit_Idx:]+fluxMask
        engyIdx, tmeIdx = np.where(pitchSlice == pitchSlice.max())
        peakEnergy_atPeakFlux.append(Energy[Energy_yLimit_Idx+engyIdx[0]])
        peakFlux.append(pitchSlice.max())
        outside_peakEnergy_atpeakFlux.append(pitchSlice.max())

    # plot the results
    axPeakE.plot(avgTimes, peakEnergy, color=plot_Colors[idx], label=rf'$\alpha = {Pitch[ptchVal]}^\circ$', marker='.', ms=plot_MarkerSize)
    # ax[3].plot(avgTimes, peakFlux, color=plot_Colors[idx], label=rf'$\alpha = {Pitch[ptchVal]}^\circ$', marker='.', ms=plot_MarkerSize)

    outside_peakFlux.append(peakFlux)
    outside_peakEnergy_atpeakFlux.append(peakEnergy_atPeakFlux)
    outside_avgtimes.append(avgTimes)

axPeakE.set_ylabel('Peak Energy [eV]',fontsize=labels_FontSize)
axPeakE.legend(fontsize=legend_FontSize)
axPeakE.set_ylim(0,Energy[Energy_yLimit_Idx])
axPeakE.set_xlabel('Time Since Launch [s]',fontsize=labels_FontSize)
# ax[3].set_ylabel('Peak Flux',fontsize=labels_FontSize)
# ax[3].legend(fontsize=legend_FontSize)
# ax[3].set_ylim(1E7, 2E9)
# ax[3].set_xlabel('Time Since Launch [s]', fontsize=labels_FontSize, labelpad=labels_Padding)
# ax[3].set_yscale('log')

# set the ticks for all plots
for ax in mainThreeAxes:
    ax.tick_params(axis='y', which='major', labelsize=tick_LabelSize, width=tick_Width, length=tick_Length)
    ax.tick_params(axis='y', which='minor', labelsize=tick_LabelSize, width=tick_Width, length=tick_Length / 2)
    ax.tick_params(axis='x', which='major', labelsize=tick_LabelSize, width=tick_Width, length=tick_Length)
    ax.tick_params(axis='x', which='minor', labelsize=tick_LabelSize, width=tick_Width, length=tick_Length / 2)
    # ax.margins(0)


# Histogram Colorbar
caxH = fig.add_axes([0.91,  0.81, 0.025, 0.22])
cbar_hist = plt.colorbar(mpl.cm.ScalarMappable(norm=histNorm, cmap=cmap_hist), cax=caxH, drawedges=False)
cbar_hist.set_label(f'Median Pitch Angle', fontsize=labels_FontSize, rotation=270)
cbar_hist.ax.get_yaxis().labelpad = 22
cbar_hist.set_ticks(ticks=[0, 20, 40, 60, 80], labels=[20*i for i in range(5)])
for l in cbar_hist.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(cbar_FontSize)


fig.align_ylabels(mainThreeAxes)

# --- Calculate Flux vs energy graphs ---

for t, tme in enumerate(STEBtimes[1:]):

    wDis = wSTEBtoPlot[t+1]
    wDispersion_key = f's{wDis}'
    lowCut, highCut = np.abs(data_dict_eepaa_high['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis - 1][0]).argmin(), np.abs(data_dict_eepaa_high['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis - 1][1]).argmin()
    Epoch_dis = deepcopy(dateTimetoTT2000(data_dict_eepaa_high['Epoch'][0][lowCut:highCut + 1], inverse=False))
    eepaa_dis_pre = deepcopy(diffEFlux[lowCut:highCut + 1])
    eepaa_dis = deepcopy(np.array(dispersionAttributes.isolationFunctions[wDispersion_key](eepaa_dis_pre, Energy, Epoch_dis)))  # pply the isolation functions found in dispersionAttributes.py
    # prepare the data by removing fillvals and applying the mask
    STEBdata = deepcopy(eepaa_dis)
    STEBdata[STEBdata < 0] = 0  # set anything below 0 = 0
    # STEBdata -= fluxMask
    # STEBdata[STEBdata < 0] = 0  # set anything below 0 = 0
    # STEBdata += fluxMask

    for idx, ptchVal in enumerate(wPitchs_to_plot):

        ptchSlice = STEBdata[:, ptchVal, Energy_yLimit_Idx:].T
        Energies = Energy[Energy_yLimit_Idx+1:]
        fluxSum = [sum(ptchSlice[engy])/len(ptchSlice[engy]) for engy in range(len(Energies))]
        subAxes[t].plot(Energies, fluxSum, color=plot_Colors[idx], label=rf'$\alpha = {Pitch[ptchVal]}^\circ$', marker='.', ms=plot_MarkerSize-7)

    # set the labels and such
    subAxes[t].legend()
    if t in [0,2]:
        subAxes[t].set_ylabel('Avg. Flux',fontsize=labels_FontSize,)
    if t in [2,3]:
        subAxes[t].set_xlabel('Energy',fontsize=labels_FontSize,)
    if t in [0,1]:
        subAxes[t].set_xticklabels([])
    subAxes[t].set_yscale('log')
    subAxes[t].set_ylim(1E6, 1E9)
    subAxes[t].tick_params(axis='y', which='major', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length)
    subAxes[t].tick_params(axis='y', which='minor', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length / 2)
    subAxes[t].tick_params(axis='x', which='major', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length)
    subAxes[t].tick_params(axis='x', which='minor', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length / 2)

# output the figure
plt.subplots_adjust(left=0.12, bottom=0.055, right=0.9, top=0.99, wspace=None, hspace=0.08)
fileOutName = rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot7\\Plot7_KeyObservations_base.png'
plt.savefig(fileOutName)
Done(start_time)
plt.close()


# fig, ax = plt.subplots()
# wSTEB = 2
# for i in range(len(outside_peakEnergy_atCounts)):
# # for i in range(1):
#     # ax.plot(outside_peakEnergy_atCounts[i][wSTEB],outside_peakCounts[i][wSTEB], color=plot_Colors[i], label=rf'$\alpha = {Pitch[i]}^\circ$', marker='.', ms=plot_MarkerSize)
#     ax.plot(outside_avgtimes[i], outside_peakEnergy_atCounts[i], color=plot_Colors[i], label=rf'$\alpha = {Pitch[i]}^\circ$', marker='.', ms=plot_MarkerSize)
# plt.legend()
# fig.suptitle(f'STEB {wSTEB+1}')
# ax.set_xlabel('time')
# ax.set_ylabel('Energy at Peak Counts')
# plt.show()