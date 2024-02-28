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
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import math
import matplotlib.colors


print(color.UNDERLINE + f'Plot4_pitchAnglePlots' + color.END)

#################
# --- TOGGLES ---
#################

# plot toggles - Overview -----------------------------
isolateSignatures = False # applies functions that remove data outside of signature
sliceEpochIndicies = {
    's1':[5934, 5940, 5946],
    's2':[5959, 5966, 5974],
    's3':[5987 - 3, 5990 - 1, 5995],
    's4':[6003 - 1, 6006 - 1, 6009 - 1],
    's5':[6014, 6018 + 1, 6021 + 3],
    's10':[6139, 6142, 6145]  # s10 The Big One on the poleward side of the aurora
}
figure_height =40
figure_width = 8
cmap = 'turbo'
cbarLow_counts, cbarHigh_counts = 1, 100
labelsFontSize = 30
TimeSinceLaunchLabelFontSize = 32


# plot toggles - Show STEB itself ----------
wDispersions = np.array([1,2,3,4,5,10])-1 # [s1, s2, s3, s4, etc] <-- Index
wPitch_Engy_vs_Time = 2 # the pitch angle index to plot for the Energy vs time plot
# colors = ['red', 'green', 'black','red', 'green', 'black','red', 'green', 'black']
sliceLineColors = np.array([['red', 'red', 'red'] for i in range(len(wDispersions))]).flatten()
# sliceLineColors = ['red', 'green', 'black','red', 'green', 'black','red', 'green', 'black']
Energy_yLimit = 1350

# plot toggles - Slices pitch angle ------------------
X_Velocity_limits, Y_Velocity_limit = [-0.5, 1.6], [-1.6, 1.6]
NoOfSlices = 3

# plot toggles - Histogram ---------------------------
countsthresh = 5
# consider only the pitch angles between -10 and 90
# [  0 1  2  3  4  5  6  7  8  9 10 ...]
# [-10 0 10 20 30 40 50 60 70 80 90 ...]
pitchAngleWidthAcceptance_lowerlimit = 2
pitchAngleWidthAcceptance_upperlimit =10 +1 # NEEDS +1 for array indexing
useCustomHistrogramColorbar = True # create custom histogram Colorbar
if useCustomHistrogramColorbar:
    tempCmap = plt.cm.turbo_r # define new colormap
    cmaplist = [tempCmap(i) for i in range(tempCmap.N)] # extract all colors from the colormap
    cmap_hist = matplotlib.colors.LinearSegmentedColormap.from_list('turboCustom',cmaplist,tempCmap.N) # create the new map
    bounds = 10*np.array([-1.5+i for i in range(11)])
    histNorm = matplotlib.colors.BoundaryNorm(bounds,tempCmap.N)
else:
    cmap_hist = 'turbo_r'
    histNorm = matplotlib.colors.Normalize(vmin=-10, vmax=90)


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
data_dict_eepaa_low = loadDictFromFile(inputFilePath=inputFiles_eepaa[1])
cbarLow,cbarHigh = cbarLow_counts, cbarHigh_counts
inputFiles_eepaa_counts = [glob('C:\Data\ACESII\L1\high\*eepaa_fullCal*')[0], glob('C:\Data\ACESII\L1\low\*eepaa_fullCal*')[0]]
data_dict_counts_high = loadDictFromFile(inputFilePath=inputFiles_eepaa_counts[0])

Done(start_time)



# --- --- --- --- --- --- --
# --- Calc Vpara & Vperp ---
# --- --- --- --- --- --- --
prgMsg('Calculating Vperp and Vparallel')
Energy = data_dict_eepaa_high['Energy'][0]
Pitch = data_dict_eepaa_high['Pitch_Angle'][0]
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
from ACESII_code.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes
dispersionTimes = dispersionAttributes.keyDispersionTimes

# the slices in time for each dispersion used
sliceTimes = {key:[data_dict_eepaa_high['Epoch'][0][val] for val in sliceEpochIndicies[key]] for key,val in sliceEpochIndicies.items()}

###############################
# --- Make the Column Plots ---
###############################
for i in range(len(wDispersions)):
    prgMsg('Beginning Plot')

    fig, ax = plt.subplots(NoOfSlices+1+1)
    fig.set_figwidth(figure_width)
    fig.set_figheight(figure_height)

    # get the data
    IndexLow, IndexHigh = np.abs(data_dict_counts_high['Epoch'][0] - dispersionTimes[wDispersions[i]][0]).argmin(), np.abs(data_dict_counts_high['Epoch'][0] - dispersionTimes[wDispersions[i]][1]).argmin()
    Epoch = EpochTo_T0_Rocket(InputEpoch=data_dict_counts_high['Epoch'][0][IndexLow:IndexHigh], T0=data_dict_counts_high['Epoch'][0][0])
    dataArray = np.array(data_dict_counts_high['eepaa'][0][IndexLow:IndexHigh])

    if isolateSignatures:
        Epoch_for_one_dispersion = np.array(Epoch) - Epoch[0]
        dataArray = dispersionAttributes.isolationFunctions[f's{wDispersions[i] + 1}'](dataArray, Energy,  Epoch_for_one_dispersion)

    # --- loop through the axis: (STEB Plot, Slice 1, slice 2, slice 3) ---
    for j in range(NoOfSlices+1+1):

        ax[j].xaxis.set_tick_params(labelsize=labelsFontSize - 8)
        ax[j].yaxis.set_tick_params(labelsize=labelsFontSize - 8)

        # --- STEB top Plot ---
        if j == 0:
            # get the formatted eepaa data - Energy vs Time
            dataToPlot = deepcopy(np.transpose(dataArray[:, wPitch_Engy_vs_Time, :]))

            # Set the background black by turning all 0 values into 1's (just for display purposes)
            for tme in range(len(dataToPlot)):
                for engy in range(len(dataToPlot[0])):
                    if dataToPlot[tme][engy] == 0:
                        dataToPlot[tme][engy] = 1

            dispersionTitleTime = pycdf.lib.tt2000_to_datetime(pycdf.lib.datetime_to_tt2000(dispersionTimes[wDispersions[i]][0]) + int((pycdf.lib.datetime_to_tt2000(dispersionTimes[wDispersions[i]][0]) - pycdf.lib.datetime_to_tt2000(dispersionTimes[wDispersions[i]][1]))/2)).strftime("%H:%M:%S.%f")[:-3]

            ax[j].set_title(f'STEB {wDispersions[i]+1}\n' + dispersionTitleTime  + ' UTC',fontsize=labelsFontSize+5)
            alfSigPlot = ax[j].pcolormesh(Epoch, Energy, dataToPlot, cmap=cmap, shading='nearest',norm='log', vmin=cbarLow, vmax=cbarHigh)

            # format the plot
            ax[j].set_ylim(Energy[-1], Energy_yLimit)
            ax[j].set_ylabel('Energy [eV]',fontsize=labelsFontSize+8)
            ax[j].set_xlabel('Time Since Launch [s]', fontsize=TimeSinceLaunchLabelFontSize, weight='bold')
            ax[j].set_yscale('log')
            ax[j].tick_params(axis='y', which='major', labelsize=labelsFontSize+10,width=3,length=16)
            ax[j].tick_params(axis='y', which='minor', labelsize=labelsFontSize+5,width=3,length=16)
            ax[j].tick_params(axis='x', which='major', labelsize=labelsFontSize, width=3, length=16)
            ax[j].tick_params(axis='x', which='minor', labelsize=labelsFontSize, width=3, length=16)

        # --- STEB Slices Plots ---
        elif j != 0 and j != NoOfSlices+1:
            # get the formatted eepaa data - Vpara vs Vperp

            # colored textbox
            timeTag = round(EpochTo_T0_Rocket(InputEpoch=[sliceTimes[f's{wDispersions[i] + 1}'][j - 1]], T0=data_dict_eepaa_high['Epoch'][0][0])[0], 2)
            props = dict(boxstyle='round', facecolor='white', alpha=1)
            ax[j].text(0.5, -1.35, f'$t_{j}$=' +f'{timeTag} s',fontsize=labelsFontSize, weight='bold', color=sliceLineColors[j-1],bbox=props, ha='center')

            # dataToPlot
            dataArray_Slice = data_dict_counts_high['eepaa'][0][np.abs(data_dict_eepaa_high['Epoch'][0] - sliceTimes[f's{wDispersions[i]+1}'][j - 1]).argmin()]

            # Set the background black by turning all 0 values into 1's (just for display purposes)
            for tme in range(len(dataArray_Slice)):
                for engy in range(len(dataArray_Slice[0])):
                    if not dataArray_Slice[tme][engy] >= 1:
                        dataArray_Slice[tme][engy] = 1

            ax[j].pcolormesh(Vperp, Vpara, dataArray_Slice, cmap=cmap, shading='nearest',norm='log', vmin=cbarLow, vmax=cbarHigh)
            ax[j].set_xlim(X_Velocity_limits[0], X_Velocity_limits[1])
            ax[j].set_ylim(Y_Velocity_limit[0], Y_Velocity_limit[1])
            ax[j].set_ylabel('V$_{\parallel}$ [10,000 km/s]', fontsize=labelsFontSize+6)
            ax[j].tick_params(axis='both', which='major', labelsize=labelsFontSize-5, width=3, length=10)
            ax[j].tick_params(axis='both', which='minor', labelsize=labelsFontSize-5, width=3, length=16)

            # if j == NoOfSlices:
            ax[j].set_xlabel('V$_{\perp}$ [10,000 km/s]', fontsize=labelsFontSize+6)
            ax[j].invert_yaxis()

            # add a verical line on the Alfvenic signature plot
            ax[0].axvline(x=timeTag, color=sliceLineColors[j-1], linewidth=5)

        # --- Pitch Angle COUNTS Histogram ---
        else:
            # construct the new dataset variable
            x = Epoch
            y = Energy
            X, Y = np.meshgrid(x, y)
            Z = deepcopy(X)

            # populate the new data
            for tme in range(len(Epoch)):
                for engy in range(len(Energy)):

                    # get the data across all pitch angles here
                    # pitchData = dataArray_counts[tme, :, engy]
                    pitchData = dataArray[tme, :, engy]

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
            cmapObj = ax[j].pcolormesh(X, Y, Z, cmap=cmap_hist, norm=histNorm, shading='nearest')
            ax[j].set_yscale('log')
            ax[j].set_ylim(Energy[-1], Energy_yLimit)
            ax[j].set_ylabel('Energy [eV]', fontsize=labelsFontSize+8)
            ax[j].set_xlabel('Time Since Launch [s]',fontsize=TimeSinceLaunchLabelFontSize, weight='bold')
            ax[j].tick_params(axis='y', which='major', labelsize=labelsFontSize + 10,width=3,length=16)
            ax[j].tick_params(axis='y', which='minor', labelsize=labelsFontSize + 5,width=3,length=8)
            ax[j].tick_params(axis='x', which='major', labelsize=labelsFontSize, width=3, length=16)
            ax[j].tick_params(axis='x', which='minor', labelsize=labelsFontSize, width=3, length=16)

    # output the figure
    pitchThreshUsed = Pitch[pitchAngleWidthAcceptance_upperlimit-1]
    fileOutName = rf'C:\Users\cfelt\OneDrive\Desktop\Paper_Photos\Plot4\\Plot4_pitchAngle_STEB_s{wDispersions[i]+1}_median_{pitchThreshUsed}deg_{countsthresh}countsThresh.png'
    plt.tight_layout(pad=1.5)
    plt.savefig(fileOutName)
    plt.close()
    Done(start_time)

    prgMsg('Creating Colorbar Plot')
    # --- Colorbar Plot ---
    fig, ax = plt.subplots(NoOfSlices + 1 + 1)
    fig.set_figwidth(figure_width+15)
    fig.set_figheight(figure_height)
    import matplotlib as mpl

    # Slices Colorbar
    axesToColorbar = ax.ravel().tolist()
    sliceNorm = mpl.colors.LogNorm(vmin=cbarLow,vmax=cbarHigh)
    cbar_slice = fig.colorbar(mpl.cm.ScalarMappable(norm=sliceNorm, cmap=cmap), ax=axesToColorbar[0:4])
    cbar_slice.ax.get_yaxis().labelpad = 70
    cbar_slice.set_label('Counts', fontsize=50, rotation=270)
    for l in cbar_slice.ax.yaxis.get_ticklabels():
        l.set_weight("bold")
        l.set_fontsize(40)

    # Histogram Colorbar
    cbar_hist = plt.colorbar(mpl.cm.ScalarMappable(norm=histNorm,cmap=cmap_hist), ax=axesToColorbar[4])
    cbar_hist.set_label(f'Median Pitch Angle', fontsize=labelsFontSize+10, rotation=270)
    cbar_hist.ax.get_yaxis().labelpad = 40
    for l in cbar_hist.ax.yaxis.get_ticklabels():
        l.set_weight("bold")
        l.set_fontsize(25)

    plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Paper_Photos\Plot4\\Plot4_pitchAngle_STEB_colorbar.png')
    # output
    Done(start_time)