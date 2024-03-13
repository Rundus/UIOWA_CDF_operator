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

from ACESII_code.myImports import *
from scipy.signal import spectrogram
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from my_matplotlib_Assets.colorbars.blue_green_White_yellow_red import blue_green_white_yellow_red_cmap
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
print(color.UNDERLINE + f'Plot2_Conjugacy' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
dpi = 500
cbarMin, cbarMax = 5E7, 1E10
cbarTickLabelSize = 14
import matplotlib

my_cmap = apl_rainbow_black0_cmap()
my_cmap.set_bad((0, 0, 0))



spectrogramCmap = blue_green_white_yellow_red_cmap()

# --- HF/LF General ---
plot_General = False
General_figure_width = 8.5 # in inches
General_figure_height = 11/2 # in inches
General_targetILat = [71.20, 73.55]
General_LabelFontSize = 14
General_TickFontSize = 10
General_PlotLineWidth = 0.5
GeneralCmap = my_cmap

# --- HF/LF Dispersive Region ---
plot_Dispersive = True
Dispersive_figure_width = 8.5
Dispersive_figure_height = 11/2
Disp_targetILat= [71.90, 72.03]
Dispersive_wPitch = 2 # 10deg pitch
Dispersive_LabelFontSize = 14
Dispersive_TickFontSize = 13
Dispersive_TickLength = 5
Dispersive_TickWidth = 2
Dispersive_LabelPadding = 10
DispersiveCmap = my_cmap
specCbarMin, specCbarMax = 1E-2, 1E1
DispersiveFreqlimits = [0, 12]
Disp_PlotLineWidth = 1

# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---

if plot_General:
    targetILat = General_targetILat
    plot_Dispersive = False
elif plot_Dispersive:
    targetILat = Disp_targetILat


prgMsg('Loading Data')
rocketAttrs, b, c = ACES_mission_dicts()

# delta B
inputMagFiles_high = glob('C:\Data\ACESII\L3\deltaB\high\*Field_Aligned*')[0]
data_dict_mag_high = loadDictFromFile(inputFilePath=inputMagFiles_high, targetVar=[targetILat,'ILat'], wKeys=['B_e', 'B_r', 'B_p', 'ILat', 'Epoch', 'Alt'])
inputMagFiles_low = glob('C:\Data\ACESII\L3\deltaB\low\*Field_Aligned*')[0]
data_dict_mag_low = loadDictFromFile(inputFilePath=inputMagFiles_low, targetVar=[targetILat,'ILat'], wKeys=['B_e', 'B_r', 'B_p', 'ILat', 'Epoch', 'Alt'])

# delta E
inputEFIFiles_low = glob('C:\Data\ACESII\L3\deltaE\low\*Field_Aligned*')[0]
data_dict_Efield_low = loadDictFromFile(inputFilePath=inputEFIFiles_low, targetVar=[targetILat,'ILat'], wKeys=['E_e', 'E_r', 'E_p', 'ILat', 'Epoch', 'Alt'])

# EEPAA Particle Data
inputEEPAA_low = glob('C:\Data\ACESII\L2\low\*eepaa_fullCal*')[0]
data_dict_eepaa_low = loadDictFromFile(inputFilePath=inputEEPAA_low, targetVar=[targetILat,'ILat'], wKeys=['Differential_Energy_Flux', 'ILat', 'Epoch', 'Alt'])
inputEEPAA_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0]
data_dict_eepaa_high = loadDictFromFile(inputFilePath=inputEEPAA_high, targetVar=[targetILat,'ILat'], wKeys=['Differential_Energy_Flux', 'ILat', 'Epoch', 'Alt'])
Done(start_time)

############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

if plot_General:

    # --- Calculate Omni-Directional Flux ---
    prgMsg('Calculating OmniFlux')

    omniDirFlux_low = np.zeros(shape=(len(data_dict_eepaa_low['Differential_Energy_Flux'][0]), len(data_dict_eepaa_low['Energy'][0])))
    for tme in range(len(data_dict_eepaa_low['Epoch'][0])):
        for engy in range(len(data_dict_eepaa_low['Energy'][0])):
            sumValues = data_dict_eepaa_low['Differential_Energy_Flux'][0][tme, :, engy]
            sumValues = sum(sumValues[sumValues>0])
            omniDirFlux_low[tme][engy] = sumValues

    omniDirFlux_high = np.zeros(shape=(len(data_dict_eepaa_high['Differential_Energy_Flux'][0]), len(data_dict_eepaa_high['Energy'][0])))
    for tme in range(len(data_dict_eepaa_high['Epoch'][0])):
        for engy in range(len(data_dict_eepaa_high['Energy'][0])):
            sumValues = data_dict_eepaa_high['Differential_Energy_Flux'][0][tme, :, engy]
            sumValues = sum(sumValues[sumValues > 0])
            omniDirFlux_high[tme][engy] = sumValues

    Done(start_time)

    fig, ax = plt.subplots(4,sharex=True,height_ratios=[2,1,2,1])
    fig.set_figwidth(General_figure_width)
    fig.set_figheight(General_figure_height)
    fig.subplots_adjust(hspace=0) # remove the space between plots

    # ---HF EEPAA---
    cmap = ax[0].pcolormesh(data_dict_eepaa_high['ILat'][0],data_dict_eepaa_high['Energy'][0],omniDirFlux_high.T, cmap=GeneralCmap,vmin=cbarMin,vmax=cbarMax, norm='log')
    # [xmin, ymin, dx, dy]
    cax = fig.add_axes([0.91, 0.11, 0.02, 0.77])
    cbar = plt.colorbar(cmap,cax=cax)
    # cbar.set_label('Omni-Dir. diff E. Flux \n' + '[cm$^{-2}$str$^{-1}$eV/eV]', rotation=-90, labelpad=20, fontsize=General_LabelFontSize)
    cbar.ax.minorticks_on()
    cbar.ax.tick_params(labelsize=cbarTickLabelSize)
    ax[0].set_ylabel('Energy [eV]', fontsize=General_LabelFontSize-3)
    ax[0].set_yscale('log')

    # --- delta B HF---
    ax[1].plot(data_dict_mag_high['ILat'][0],data_dict_mag_high['B_e'][0],color='blue',linewidth=General_PlotLineWidth)
    ax[1].set_ylabel('$\delta B_{e}$', fontsize=General_LabelFontSize, color='blue')
    ax[1].tick_params(axis='y', colors='blue')
    ax[1].set_ylim(-8,8)

    # --- LF EEPAA---
    cmap = ax[2].pcolormesh(data_dict_eepaa_low['ILat'][0], data_dict_eepaa_low['Energy'][0], omniDirFlux_low.T, cmap=GeneralCmap, vmin=cbarMin, vmax=cbarMax, norm='log')
    ax[2].set_ylabel('Energy [eV]', fontsize=General_LabelFontSize-3)
    ax[2].set_yscale('log')

    # --- delta B LF---
    ax[3].plot(data_dict_mag_low['ILat'][0], data_dict_mag_low['B_e'][0], color='blue',linewidth=General_PlotLineWidth)
    ax[3].set_ylabel('$\delta B_{e}$', fontsize=General_LabelFontSize, color='blue')
    ax[3].set_xlabel('ILat [deg]', fontsize=General_LabelFontSize)
    ax[3].tick_params(axis='y', colors='blue')
    ax[3].set_ylim(-8, 8)

    for i in range(4):
        ax[i].margins(0)
        ax[i].set_xlim(targetILat[0],targetILat[1])

    # --- SHOW PLOT ---
    plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Paper_Photos\Plot2\Plot2_ConjugacyStack.png', dpi=dpi)


if plot_Dispersive:

    for wRocket in [4,5]:
        magDicts = [data_dict_mag_high,data_dict_mag_low]
        magDict = magDicts[wRocket-4]
        eepaaDicts = [data_dict_eepaa_high,data_dict_eepaa_low]
        eepaaDict = eepaaDicts[wRocket-4]



        # --- Calculate Spectrogram ---
        spectrogramData = magDict['B_e'][0]
        windowType, npersegN, scalingType = 'hann', 64, 'spectrum'  # spectrogram toggles
        overlap = int(npersegN * (7 / 8))  # hanning filter overlap
        f, t, Sxx = spectrogram(spectrogramData,
                                fs=128,
                                window=windowType,
                                nperseg=npersegN,  # note: if ==None default size is 256
                                noverlap=overlap,
                                scaling=scalingType)  # scaling = density or scaling = spectrum

        # - determine a new ILat variable for the spectrogram -
        # first determine the times of the spectrogram
        startTime = pycdf.lib.datetime_to_tt2000(magDict['Epoch'][0][0])
        specTempTimes = [pycdf.lib.tt2000_to_datetime(int(startTime + 1E9*tme)) for tme in t]
        specIndicies = [np.abs(magDict['Epoch'][0] - specTempTimes[k]).argmin() for k in range(len(specTempTimes))]
        specILats = magDict['ILat'][0][specIndicies]
        specAlts = magDict['Alt'][0][specIndicies]
        specTimes = magDict['Epoch'][0][specIndicies]

        # --- get the total Differential Energy FLux over all energies ---
        totalDirFlux = np.zeros(shape=(len(eepaaDict['Differential_Energy_Flux'][0]), len(eepaaDict['Pitch_Angle'][0])))
        for tme in range(len(eepaaDict['Epoch'][0])):
            for ptch in range(len(eepaaDict['Pitch_Angle'][0])):
                sumValues = eepaaDict['Differential_Energy_Flux'][0][tme, ptch, :]
                sumValues = sum(sumValues[sumValues > 0])
                totalDirFlux[tme][ptch] = sumValues

        # --- PLOT EVERYTHING ---
        fig, ax = plt.subplots(4, sharex=True,height_ratios=[2,2,1,1])
        fig.set_figwidth(Dispersive_figure_width)
        fig.set_figheight(Dispersive_figure_height)
        fig.subplots_adjust(top=0.97,hspace=0.1)  # remove the space between plots

        # --- EEPAA Data 10deg ---
        cmap = ax[0].pcolormesh(eepaaDict['ILat'][0], eepaaDict['Energy'][0], eepaaDict['Differential_Energy_Flux'][0][:,Dispersive_wPitch,:].T, cmap=GeneralCmap, vmin=cbarMin, vmax=cbarMax, norm='log')
        ax[0].set_ylabel(rf'P.A. = {eepaaDict["Pitch_Angle"][0][Dispersive_wPitch]}$^\circ$' + '\n Energy [eV]', fontsize=Dispersive_LabelFontSize,labelpad=Dispersive_LabelPadding)
        ax[0].set_yscale('log')
        ax[0].set_ylim(33,1150)
        ax[0].tick_params(axis='both', labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)

        # --- EEPAA Data All Pitch ---
        ax[1].pcolormesh(eepaaDict['ILat'][0], eepaaDict['Pitch_Angle'][0], totalDirFlux.T, cmap=GeneralCmap, vmin=cbarMin, vmax=cbarMax, norm='log')
        ax[1].set_ylabel('P. A.\n[deg]', fontsize=Dispersive_LabelFontSize,labelpad=Dispersive_LabelPadding)
        ax[1].set_ylim(0,181)
        ax[1].margins(0)
        ax[1].tick_params(axis='both', labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)
        yticks = [0, 60, 120, 180]
        ax[1].set_yticks(yticks)
        ax[1].set_yticklabels([str(tick) for tick in yticks])

        cax = fig.add_axes([0.91, 0.415, 0.02, 0.556])
        cbar = plt.colorbar(cmap, cax=cax)
        cbar.ax.minorticks_on()
        cbar.ax.tick_params(labelsize=cbarTickLabelSize)

        # --- Be/Er ---
        ax[2].plot(magDict['ILat'][0], magDict['B_e'][0], color='blue', linewidth=Disp_PlotLineWidth)
        ax[2].set_ylabel('$\delta B_{e}$\n[nT]', fontsize=Dispersive_LabelFontSize, color='blue',labelpad=Dispersive_LabelPadding+10)
        ax[2].tick_params(axis='y', colors='blue')
        ax[2].set_ylim(-8,8)
        ax[2].tick_params(axis='both', labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)

        # Er
        if wRocket == 4:
            axEr = ax[2].twinx()
            axEr.set_ylabel('E-Field\nN/A', fontsize=Dispersive_LabelFontSize, color='red', rotation=-90,labelpad=Dispersive_LabelPadding+15)
            axEr.tick_params(axis='y', colors='red',labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)
            axEr.set_ylim(-8, 8)
        elif wRocket == 5:
            axEr = ax[2].twinx()
            axEr.plot(data_dict_Efield_low['ILat'][0],data_dict_Efield_low['E_r'][0], color='red', linewidth=Disp_PlotLineWidth)
            axEr.set_ylabel('$\delta E_{r}$\n[mV/m]', fontsize=Dispersive_LabelFontSize, color='red', rotation=-90,labelpad=Dispersive_LabelPadding+15)
            axEr.tick_params(axis='y', colors='red',labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)
            axEr.set_ylim(-8, 8)


        # --- Be Spectrogram ---
        cmap = ax[3].pcolormesh(specILats, f, Sxx, shading='nearest', vmin=specCbarMin, vmax=specCbarMax, cmap=spectrogramCmap, norm='log')
        ax[3].set_ylabel('Freq.\n[Hz]', fontsize=Dispersive_LabelFontSize, labelpad=Dispersive_LabelPadding+5)
        ax[3].set_xlabel('ILat [deg]\nAlt [km]\ntime (UTC)', fontsize=Dispersive_TickFontSize-3, weight='bold')
        ax[3].xaxis.set_label_coords(-0.105, -0.17)
        ax[3].set_ylim(DispersiveFreqlimits[0],DispersiveFreqlimits[1])

        # spectrogram Ticks
        yticks = [0, 4, 8, 12]
        ax[3].set_yticks(yticks)
        ax[3].set_yticklabels([str(tick) for tick in yticks])

        xticks_iLat = ax[3].get_xticks()
        xtick_indicies = np.array([np.abs(magDict['ILat'][0] - tick ).argmin() for tick in xticks_iLat])
        ILat_ticks = [str(round(tick,2)) for tick in xticks_iLat]
        Alt_ticks = [str(round(tick,1)) for tick in  magDict['Alt'][0][xtick_indicies]]
        time_ticks =  [tick.strftime("%H:%M:%S") for tick in  magDict['Epoch'][0][xtick_indicies]]
        tickLabels = [ f'{ILat_ticks[k]}\n{Alt_ticks[k]}\n{time_ticks[k]}' for k in range(len(xtick_indicies)) ]
        ax[3].set_xticklabels(tickLabels)
        ax[3].tick_params(axis='y', labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)
        ax[3].tick_params(axis='x', labelsize=Dispersive_TickFontSize-2, length=Dispersive_TickLength, width=Dispersive_TickWidth)

        # spectrogram colorbar
        cax = fig.add_axes([0.91, 0.11, 0.02, 0.135])
        cbar = plt.colorbar(cmap, cax=cax)
        cbar.ax.minorticks_on()
        cbar.ax.tick_params(labelsize=cbarTickLabelSize)

        # output the figure
        outputMod = 'HF' if wRocket==4 else 'LF'
        plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot2\Plot2_{outputMod}dispersive.png', dpi=dpi)



