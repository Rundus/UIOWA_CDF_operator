# --- Plots4_pitchAnglePlots.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Recreation of DiffEFlux pitch angle Plots, focusing on
# a few particle signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import math
import matplotlib as mpl
from ACESII_code.class_var_func import EpochTo_T0_Rocket
from ACESII_code.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes



print(color.UNDERLINE + f'Plot4_pitchAnglePlots' + color.END)

#################
# --- TOGGLES ---
#################

# plot toggles - Overview -----------------------------
sliceEpochIndicies = {
    's1':[5934, 5940, 5946],
    's2':[5959, 5966, 5974],
    's3':[5987 - 3, 5990 - 1, 5995],
    's4':[6003 - 1, 6006 - 1, 6009 - 1],
    's5':[6014, 6018 + 1, 6021 + 3],
    's10':[6139, 6142, 6145]  # s10 The Big One on the poleward side of the aurora
}
# figure_height = (11.5)
# figure_width = (9.9318)
figure_height = (13)
figure_width = (10)

cmap = 'turbo'
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
cmap = apl_rainbow_black0_cmap()

labelPadding = -0.25
textFontSize = 10
titleFontSize = 18
labelsFontSize = 10
SidelabelAdjust = 2.6
tickFontSize = 10
tickWidth = 1
tickLength = 2
lineWidth = 2
cbarFont = 15
dpi = 1200


# plot toggles - Show STEB itself ----------
useCounts = False # if False --> use diffEFlux
cbarLow_counts, cbarHigh_counts = 1, 100
cbarLow_diffEFlux, cbarHigh_diffEFlux = 5E6, 1E9
wDispersions = np.array([2,3,4,5])-1 # [s1, s2, s3, s4, etc] <-- Index
wPitch_Engy_vs_Time = 2 # the pitch angle index to plot for the Energy vs time plot
# colors = ['red', 'green', 'black','red', 'green', 'black','red', 'green', 'black']
sliceLineColors = np.array([['black', 'black', 'black'] for i in range(len(wDispersions))]).flatten()
# sliceLineColors = ['red', 'green', 'black','red', 'green', 'black','red', 'green', 'black']
Energy_yLimit = 1350

# plot toggles - Slices pitch angle ------------------
X_Velocity_limits, Y_Velocity_limit = [-0.5, 1.6], [-1.6, 1.6]
NoOfSlices = 3

# plot toggles - Histogram ---------------------------
countsthresh = 2
# consider only the pitch angles between -10 and 90
# [  0 1  2  3  4  5  6  7  8  9 10 ...]
# [-10 0 10 20 30 40 50 60 70 80 90 ...]
pitchAngleWidthAcceptance_lowerlimit = 2
pitchAngleWidthAcceptance_upperlimit =10 +1 # NEEDS +1 for array indexing
useCustomHistrogramColorbar = True # create custom histogram Colorbar
if useCustomHistrogramColorbar:
    tempCmap = plt.cm.turbo_r # define new colormap
    cmaplist = [tempCmap(i) for i in range(tempCmap.N)] # extract all colors from the colormap
    cmap_hist = mpl.colors.LinearSegmentedColormap.from_list('turboCustom',cmaplist,tempCmap.N) # create the new map
    bounds = 10*np.array([-1.5+i for i in range(11)])
    histNorm = mpl.colors.BoundaryNorm(bounds,tempCmap.N)
else:
    cmap_hist = 'turbo_r'
    histNorm = mpl.colors.Normalize(vmin=-10, vmax=90)


# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---

prgMsg('Loading Data')

dataKey = 'eepaa'
labelNam = 'Counts'

# Attitude data (for geomagnetic lat/long info)
inputFiles_Attitude = [glob(r'C:\Data\ACESII\attitude\high\*Attitude_Solution*')[0], glob(r'C:\Data\ACESII\attitude\low\*Attitude_Solution*')[0]]
data_dict_attitude_high = loadDictFromFile(inputFilePath=inputFiles_Attitude[0])
data_dict_attitude_low = loadDictFromFile(inputFilePath=inputFiles_Attitude[1])

# EEPAA Particle Data
inputFiles_eepaa = [glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0], glob('C:\Data\ACESII\L2\low\*eepaa_fullCal*')[0]]
data_dict_eepaa_high = loadDictFromFile(inputFilePath=inputFiles_eepaa[0])

inputFiles_eepaa_counts = [glob('C:\Data\ACESII\L1\high\*eepaa_fullCal*')[0], glob('C:\Data\ACESII\L1\low\*eepaa_fullCal*')[0]]
data_dict_counts_high = loadDictFromFile(inputFilePath=inputFiles_eepaa_counts[0])

# Set the cbar limits
if useCounts:
    cbarLow, cbarHigh = cbarLow_counts, cbarHigh_counts
else:
    cbarLow, cbarHigh = cbarLow_diffEFlux, cbarHigh_diffEFlux

Done(start_time)


# --- --- --- --- --- --- --
# --- Calc Vpara & Vperp ---
# --- --- --- --- --- --- --
prgMsg('Calculating Vperp and Vparallel')
Energy = data_dict_counts_high['Energy'][0]
Pitch = data_dict_counts_high['Pitch_Angle'][0]
countsTemp = data_dict_counts_high['eepaa'][0][6000]
Vperp = deepcopy(countsTemp)
Vpara = deepcopy(countsTemp)

for ptch in range(len(countsTemp)):
    for engy in range(len(countsTemp[0])):
        Emag = np.sqrt(2*q0*Energy[engy]/(m_e))
        Vperp[ptch][engy] = np.sin(np.radians(Pitch[ptch]))*Emag
        Vpara[ptch][engy] = np.cos(np.radians(Pitch[ptch]))*Emag

Vpara, Vperp = np.array(Vpara)/(10000*1000), np.array(Vperp)/(10000*1000)
Done(start_time)


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

# define (in time) where the STEBS occur, this has been done already

dispersionTimes = dispersionAttributes.keyDispersionTimes
# dispersionTimes = dispersionAttributes.keyDispersionDeltaT

# the slices in time for each dispersion used
sliceTimes = {key:[data_dict_eepaa_high['Epoch'][0][val] for val in sliceEpochIndicies[key]] for key,val in sliceEpochIndicies.items()}

prgMsg('Beginning Plot')
fig, ax = plt.subplots(nrows=NoOfSlices+1+1, ncols=len(wDispersions))
fig.set_size_inches(figure_width,figure_height)

# fig.subplots_adjust(hspace=0) # remove the space between plots



for colIndx in range(len(wDispersions)):

    # get the data
    IndexLow, IndexHigh = np.abs(data_dict_counts_high['Epoch'][0] - dispersionTimes[wDispersions[colIndx]][0]).argmin(), np.abs(data_dict_counts_high['Epoch'][0] - dispersionTimes[wDispersions[colIndx]][1]).argmin()
    Epoch = EpochTo_T0_Rocket(InputEpoch=data_dict_counts_high['Epoch'][0][IndexLow:IndexHigh], T0=data_dict_counts_high['Epoch'][0][0])
    dataArray = np.array(data_dict_counts_high['eepaa'][0][IndexLow:IndexHigh]) if useCounts else np.array(data_dict_eepaa_high['Differential_Energy_Flux'][0][IndexLow:IndexHigh])


    # --- loop through the axis: (STEB Plot, Slice 1, slice 2, slice 3) ---
    for rowIndx in range(NoOfSlices+1+1):

        ax[rowIndx, colIndx].xaxis.set_tick_params(labelsize=tickFontSize)
        ax[rowIndx, colIndx].yaxis.set_tick_params(labelsize=tickFontSize)

        # --- STEB top Plot ---
        if rowIndx == 0:

            if useCounts:
                # get the formatted eepaa data - Energy vs Time
                dataToPlot = deepcopy(np.transpose(dataArray[:, wPitch_Engy_vs_Time, :]))

                # Set the background black by turning all 0 values into 1's (just for display purposes)
                for tme in range(len(dataToPlot)):
                    for engy in range(len(dataToPlot[0])):
                        if dataToPlot[tme][engy] == 0:
                            dataToPlot[tme][engy] = 1
            else:
                dataToPlot = deepcopy(np.transpose(dataArray[:, wPitch_Engy_vs_Time, :]))


            dispersionTitleTime = pycdf.lib.tt2000_to_datetime(pycdf.lib.datetime_to_tt2000(dispersionTimes[wDispersions[colIndx]][0]) + int((pycdf.lib.datetime_to_tt2000(dispersionTimes[wDispersions[colIndx]][1]) - pycdf.lib.datetime_to_tt2000(dispersionTimes[wDispersions[colIndx]][0]))/2)).strftime("%H:%M:%S.%f")[:-3]

            ax[rowIndx,colIndx].set_title(f'STEB {wDispersions[colIndx]+1}\n' + dispersionTitleTime  + ' UTC',fontsize=titleFontSize-4, weight='bold')
            alfSigPlot = ax[rowIndx,colIndx].pcolormesh(Epoch, Energy, dataToPlot, cmap=cmap, shading='nearest',norm='log', vmin=cbarLow, vmax=cbarHigh)

            # format the plot
            ax[rowIndx, colIndx].set_ylim(Energy[-1], Energy_yLimit)

            if colIndx == 0:
                ax[rowIndx, colIndx].set_ylabel('P.A. = 10$^{\circ}$ \n Energy [eV]',fontsize=labelsFontSize+SidelabelAdjust+2)

            ax[rowIndx, colIndx].set_xlabel('Time Since Launch [s]', fontsize=labelsFontSize+SidelabelAdjust-2, weight='bold',labelpad=labelPadding)
            ax[rowIndx, colIndx].set_yscale('log')
            ax[rowIndx, colIndx].tick_params(axis='y', which='major', labelsize=tickFontSize, width=tickWidth, length=tickLength)
            ax[rowIndx, colIndx].tick_params(axis='y', which='minor', labelsize=tickFontSize, width=tickWidth, length=tickLength/2)
            ax[rowIndx, colIndx].tick_params(axis='x', which='major', labelsize=tickFontSize, width=tickWidth, length=tickLength)
            ax[rowIndx, colIndx].tick_params(axis='x', which='minor', labelsize=tickFontSize, width=tickWidth, length=tickLength/2)

        # --- STEB Slices Plots ---
        elif rowIndx != 0 and rowIndx != NoOfSlices+1:
            # get the formatted eepaa data - Vpara vs Vperp

            # colored textbox
            timeTag = round(EpochTo_T0_Rocket(InputEpoch=[sliceTimes[f's{wDispersions[colIndx] + 1}'][rowIndx - 1]], T0=data_dict_eepaa_high['Epoch'][0][0])[0], 2)
            props = dict(boxstyle='round', facecolor='white', alpha=1)
            ax[rowIndx, colIndx].text(0.5, -1.3, f'$t_{rowIndx}$=' +f'{timeTag} s',fontsize=labelsFontSize, weight='bold', color=sliceLineColors[rowIndx-1], bbox=props, ha='center')

            # dataToPlot
            dataArray_Slice = data_dict_counts_high['eepaa'][0][np.abs(data_dict_eepaa_high['Epoch'][0] - sliceTimes[f's{wDispersions[colIndx]+1}'][rowIndx - 1]).argmin()] if useCounts else data_dict_eepaa_high['Differential_Energy_Flux'][0][np.abs(data_dict_eepaa_high['Epoch'][0] - sliceTimes[f's{wDispersions[colIndx]+1}'][rowIndx - 1]).argmin()]

            # Set the background black by turning all 0 values into 1's (just for display purposes)
            for tme in range(len(dataArray_Slice)):
                for engy in range(len(dataArray_Slice[0])):
                    if not dataArray_Slice[tme][engy] >= 1:
                        dataArray_Slice[tme][engy] = 1

            ax[rowIndx,colIndx].pcolormesh(Vperp, Vpara, dataArray_Slice, cmap=cmap, shading='nearest',norm='log', vmin=cbarLow, vmax=cbarHigh)
            ax[rowIndx,colIndx].set_xlim(X_Velocity_limits[0], X_Velocity_limits[1])
            ax[rowIndx,colIndx].set_ylim(Y_Velocity_limit[0], Y_Velocity_limit[1])
            if colIndx == 0:
                ax[rowIndx,colIndx].set_ylabel('V$_{\parallel}$ [10,000 km/s]', fontsize=labelsFontSize+SidelabelAdjust+2,labelpad=labelPadding+5)
            ax[rowIndx,colIndx].tick_params(axis='both', which='major', labelsize=tickFontSize, width=tickWidth, length=tickLength)
            ax[rowIndx,colIndx].tick_params(axis='both', which='minor', labelsize=tickFontSize, width=tickWidth, length=tickLength/2)

            # if j == NoOfSlices:
            ax[rowIndx,colIndx].set_xlabel('V$_{\perp}$ [10,000 km/s]', fontsize=labelsFontSize + SidelabelAdjust-2,labelpad=labelPadding)
            ax[rowIndx,colIndx].invert_yaxis()

            # add a verical line on the Alfvenic signature plot
            ax[0, colIndx].axvline(x=timeTag, color=sliceLineColors[rowIndx-1], linewidth=lineWidth, linestyle='--',alpha=0.6)


        else:
            # --- Pitch Angle COUNTS Histogram ---
            # construct the new dataset variable
            x = Epoch
            y = Energy
            X, Y = np.meshgrid(x, y)
            Z = deepcopy(X)
            countsArray = np.array(data_dict_counts_high['eepaa'][0][IndexLow:IndexHigh])

            # populate the new data
            for tme in range(len(Epoch)):
                for engy in range(len(Energy)):

                    # get the data across all pitch angles here
                    # pitchData = dataArray_counts[tme, :, engy]
                    pitchData = countsArray[tme, :, engy]

                    # Find errors in data and eliminate them from the median calculation
                    for k, val in enumerate(pitchData):
                        if val < countsthresh: # theshold away small COUNT values
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
            cmapObj = ax[rowIndx,colIndx].pcolormesh(X, Y, Z, cmap=cmap_hist, norm=histNorm, shading='nearest')
            ax[rowIndx,colIndx].set_yscale('log')
            ax[rowIndx,colIndx].set_ylim(Energy[-1], Energy_yLimit)
            if colIndx == 0:
                ax[rowIndx,colIndx].set_ylabel('Energy [eV]', fontsize=labelsFontSize+SidelabelAdjust+2)
            ax[rowIndx,colIndx].set_xlabel('Time Since Launch [s]',fontsize=labelsFontSize+SidelabelAdjust-2, weight='bold',labelpad=labelPadding)
            ax[rowIndx,colIndx].tick_params(axis='y', which='major', labelsize=labelsFontSize, width=tickWidth,length=tickLength)
            ax[rowIndx,colIndx].tick_params(axis='y', which='minor', labelsize=labelsFontSize, width=tickWidth,length=tickLength/2)
            ax[rowIndx,colIndx].tick_params(axis='x', which='major', labelsize=labelsFontSize, width=tickWidth, length=tickLength)
            ax[rowIndx,colIndx].tick_params(axis='x', which='minor', labelsize=labelsFontSize, width=tickWidth, length=tickLength/2)



plt.tight_layout(w_pad=-0.5,h_pad=0.15,rect=[0,0,0.91,1])

prgMsg('Creating Colorbar Plot')

# --- Counts Colorbar ---
cax = fig.add_axes([0.9, 0.23, 0.025, 0.723])
cbar = plt.colorbar(alfSigPlot, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=tickFontSize)
if useCounts:
    cbar.set_label('Counts', fontsize=labelsFontSize+10, rotation=270)
else:
    cbar.set_label('cm$^{-2}$str$^{-1}$s$^{-1}$eV/eV', fontsize=labelsFontSize + 8, rotation=270)

cbar.ax.get_yaxis().labelpad = 17
for l in cbar.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(cbarFont)

# Slices Colorbar
# axesToColorbar = ax.ravel().tolist()
# sliceNorm = mpl.colors.LogNorm(vmin=cbarLow,vmax=cbarHigh)
# cbar_slice = fig.colorbar(mpl.cm.ScalarMappable(norm=sliceNorm, cmap=cmap), ax=axesToColorbar[0:4])
# cbar_slice.ax.get_yaxis().labelpad = 70
# cbar_slice.set_label('Counts', fontsize=50, rotation=270)
# for l in cbar_slice.ax.yaxis.get_ticklabels():
#     l.set_weight("bold")
#     l.set_fontsize(40)

# Histogram Colorbar
caxH = fig.add_axes([0.9, 0.040, 0.025, 0.158])
cbar_hist = plt.colorbar(mpl.cm.ScalarMappable(norm=histNorm, cmap=cmap_hist), cax=caxH, drawedges=False)
cbar_hist.set_label(f'Median Pitch Angle', fontsize=labelsFontSize+8, rotation=270)
cbar_hist.ax.get_yaxis().labelpad = 22

cbar_hist.set_ticks(ticks=[0, 20, 40, 60, 80], labels=[20*i for i in range(5)])

for l in cbar_hist.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(cbarFont)


# output the figure
pitchThreshUsed = Pitch[pitchAngleWidthAcceptance_upperlimit-1]
fileOutName = rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot4\\Plot4_pitchAngle_median_{pitchThreshUsed}deg_{countsthresh}countsThresh.png'

plt.savefig(fileOutName)
Done(start_time)