# --- Prove_Alfvenic.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Do an analysis on the waves observed to determine
# if they are alfven waves.


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
# --- TOGGLES ---
# --- --- --- ---
modifier = ''
inputPath_modifier_B = 'science\deltaB' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_E = 'science\deltaE'
inputPath_modifier_attitude = 'attitude'
inputPath_modifier_density = 'science\langmuir'
inputPath_modifier_AlfvenSpeed = 'science\AlfvenSpeed_rkt'
inputPath_modifier_Poynting = 'science\PoyntingFlux'

# --- Reduce data ---
reduceData = False
targetTimes = [dt.datetime(2022,11,20,17,24,10,00),dt.datetime(2022,11,20,17,25,50,00)]

# --- PLOTTING ---
SECTION_ShowWavesAreFieldAligned = True
SECTION_deltaBdeltaEPoyntingDensity = False
SECTION_deltaBdeltaEFFTspeed = False if not reduceData else False
SECTION_spectrogramData = False

spinFreq = [0.6442441031179716, 0.55]

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import InterpolateDataDict,dateTimetoTT2000,butter_filter
from numpy.fft import rfft, fftfreq
from scipy.signal import spectrogram

def Prove_Alfvenic():

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()

    inputFile_attitude_low = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[1]}{modifier}\*.cdf')[0]
    inputFile_attitude_high = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[0]}{modifier}\*.cdf')[0]

    inputFile_deltaE_low = glob(f'{rocketFolderPath}{inputPath_modifier_E}\{fliers[1]}{modifier}\*Field_Aligned*')[0]

    inputFile_deltaB_high = glob(f'{rocketFolderPath}{inputPath_modifier_B}\{fliers[0]}{modifier}\*Field_Aligned*')[0]
    inputFile_deltaB_low = glob(f'{rocketFolderPath}{inputPath_modifier_B}\{fliers[1]}{modifier}\*Field_Aligned*')[0]

    inputFile_density_low = glob(f'{rocketFolderPath}{inputPath_modifier_density}\{fliers[1]}{modifier}\*fixed*')[0]
    inputFile_density_high = glob(f'{rocketFolderPath}{inputPath_modifier_density}\{fliers[0]}{modifier}\*fixed*')[0]

    inputFile_AlfvenSpeed_low = glob(f'{rocketFolderPath}{inputPath_modifier_AlfvenSpeed}\{fliers[1]}{modifier}\*.cdf*')[0]
    inputFile_AlfvenSpeed_high = glob(f'{rocketFolderPath}{inputPath_modifier_AlfvenSpeed}\{fliers[0]}{modifier}\*.cdf*')[0]

    inputFile_Poynting_low = glob(f'{rocketFolderPath}{inputPath_modifier_Poynting}\{fliers[1]}{modifier}\*Field_Aligned*')[0]
    inputFile_Poynting_high = glob(f'{rocketFolderPath}{inputPath_modifier_Poynting}\{fliers[0]}{modifier}\*Field_Aligned*')[0]

    # --- get the data from the Attitude files ---
    prgMsg(f'Loading data from {inputPath_modifier_attitude} Files')
    data_dict_attitude_low = loadDictFromFile(inputFile_attitude_low, {},reduceData,targetTimes)
    data_dict_attitude_high = loadDictFromFile(inputFile_attitude_high, {},reduceData,targetTimes)
    Done(start_time)

    # --- get the data from the E-Field files ---
    prgMsg(f'Loading data from {inputPath_modifier_E} Files')
    data_dict_deltaE_low = loadDictFromFile(inputFile_deltaE_low, {},reduceData,targetTimes)
    Done(start_time)

    # --- get the data from the B-Field files ---
    prgMsg(f'Loading data from {inputPath_modifier_B} Files')
    data_dict_deltaB_low = loadDictFromFile(inputFile_deltaB_low, {},reduceData,targetTimes)
    data_dict_deltaB_high = loadDictFromFile(inputFile_deltaB_high, {},reduceData,targetTimes)
    Done(start_time)

    # --- get the data from the Alfven files ---
    prgMsg(f'Loading data from {inputPath_modifier_AlfvenSpeed} Files')
    data_dict_Alfven_low = loadDictFromFile(inputFile_AlfvenSpeed_low, {},reduceData,targetTimes)
    data_dict_Alfven_high = loadDictFromFile(inputFile_AlfvenSpeed_high, {},reduceData,targetTimes)
    Done(start_time)

    # --- get the data from the Density files ---
    prgMsg(f'Loading data from {inputPath_modifier_density} Files')
    data_dict_density_low = loadDictFromFile(inputFile_density_low, {}, reduceData, targetTimes)
    data_dict_density_high = loadDictFromFile(inputFile_density_high, {}, reduceData, targetTimes)
    Done(start_time)

    DataDicts = [deepcopy(data_dict_deltaE_low),deepcopy(data_dict_deltaB_high),deepcopy(data_dict_deltaB_low),deepcopy(data_dict_Alfven_high)]

    # --- get the data from the Poynting files ---
    prgMsg(f'Loading data from {inputPath_modifier_Poynting} Files')
    data_dict_Poynting_low = loadDictFromFile(inputFile_Poynting_low, {},reduceData,targetTimes)
    data_dict_Poynting_high = loadDictFromFile(inputFile_Poynting_high, {},reduceData,targetTimes)
    Done(start_time)


    #######################################
    # --- GET THE DATA ON THE MAG TIMES ---
    #######################################

    prgMsg('Reducing Data onto mag times')
    data_dict_densityInterp_low = InterpolateDataDict(InputDataDict=data_dict_density_low,
                                                  InputEpochArray=dateTimetoTT2000(data_dict_density_low['Epoch'][0],False),
                                                  wKeys=['ni', 'Epoch'],
                                                  targetEpochArray=dateTimetoTT2000(data_dict_deltaB_low['Epoch'][0], False))

    # NOTE: YOU need to have reducedData == True and use roughly the whole flight for E/B ratio plot. OTHERWISE
    # the last handful of values in data_Dict_Density['fixed_epoch'] throws and error since the epoch == fillvals

    data_dict_densityInterp_high = InterpolateDataDict(InputDataDict=data_dict_density_high,
                                                      InputEpochArray=dateTimetoTT2000(data_dict_density_high['Epoch'][0], False),
                                                      wKeys=['ni', 'Epoch'],
                                                      targetEpochArray=dateTimetoTT2000(data_dict_deltaB_high['Epoch'][0], False))

    data_dict_deltaEInterp_low = InterpolateDataDict(InputDataDict=data_dict_deltaE_low,
                                                      InputEpochArray=dateTimetoTT2000(data_dict_deltaE_low['Epoch'][0], False),
                                                      wKeys=[],
                                                      targetEpochArray=dateTimetoTT2000(data_dict_deltaB_low['Epoch'][0], False))
    Done(start_time)



    ###############################################
    # --- SHOW THE WAVES ARE VERY FIELD ALIGNED ---
    ###############################################
    if SECTION_ShowWavesAreFieldAligned:
        # make a plot of the e,p,r coordinates for High and Low Flyer showing:
        # [1] The wave components are mostly oscillating in the perpendicular direction
        # [2] The pointing flux is primarily towards the earth

        # Low Flyer Plot
        fig, ax = plt.subplots(nrows=3,ncols=2)
        fig.suptitle('ACESII 36364 Components')
        ax[0, 0].plot(data_dict_deltaB_low['Epoch'][0], data_dict_deltaB_low['B_e'][0])
        ax[0, 0].set_ylabel('$\delta B_{e}$')
        ax[1, 0].plot(data_dict_deltaB_low['Epoch'][0], data_dict_deltaB_low['B_p'][0])
        ax[1, 0].set_ylabel('$\delta B_{p}$')
        ax[2, 0].plot(data_dict_deltaB_low['Epoch'][0], data_dict_deltaB_low['B_r'][0])
        ax[2, 0].set_ylabel('$\delta B_{r}$')

        ax[0, 1].plot(data_dict_deltaE_low['Epoch'][0], data_dict_deltaE_low['E_e'][0])
        ax[0, 1].set_ylabel('$\delta E_{e}$')
        ax[1, 1].plot(data_dict_deltaE_low['Epoch'][0], data_dict_deltaE_low['E_p'][0])
        ax[1, 1].set_ylabel('$\delta E_{p}$')
        ax[2, 1].plot(data_dict_deltaE_low['Epoch'][0], data_dict_deltaE_low['E_r'][0])
        ax[2, 1].set_ylabel('$\delta E_{r}$')
        plt.show()


        fig, ax = plt.subplots(nrows=3,ncols=1)
        fig.suptitle('ACESII 36359 Components')
        ax[0].plot(data_dict_deltaB_high['Epoch'][0], data_dict_deltaB_high['B_e'][0])
        ax[0].set_ylabel('B_e')
        ax[1].plot(data_dict_deltaB_high['Epoch'][0], data_dict_deltaB_high['B_p'][0])
        ax[0].set_ylabel('B_p')
        ax[2].plot(data_dict_deltaB_high['Epoch'][0], data_dict_deltaB_high['B_r'][0])
        ax[0].set_ylabel('B_r')
        plt.show()

    ############################################
    # --- deltaB,deltaE,Poynting,DensityPlot ---
    ############################################

    if SECTION_deltaBdeltaEPoyntingDensity:

        ###### LOW FLYER ######
        fig, ax = plt.subplots(nrows=4,ncols=1)
        plt.subplots_adjust(wspace=0, hspace=0)
        fig.suptitle('ACESII 36364')

        # B_e, E_r
        ax[0].plot(data_dict_deltaEInterp_low['Epoch'][0], data_dict_deltaEInterp_low['E_r'][0], color='red', label='E (+r)')
        ax[0].set_ylabel('E (mV/m)')
        ax[0].set_ylim(-3, 3)
        ax[0].legend(loc='upper left')
        ax[0].grid(True)
        axBe = ax[0].twinx()
        axBe.plot(data_dict_deltaB_low['Epoch'][0], data_dict_deltaB_low['B_e'][0], color='tab:blue', label='B (+e)')
        axBe.set_ylabel('B (nT)')
        axBe.set_ylim(-6, 6)
        axBe.legend(loc='upper right')

        # B_r, E_e
        ax[1].plot(data_dict_deltaEInterp_low['Epoch'][0], data_dict_deltaEInterp_low['E_e'][0], color='red', label='E (+e)')
        ax[1].set_ylabel('E (mV/m)')
        ax[1].set_ylim(-3, 3)
        ax[1].grid(True)
        ax[1].legend(loc='upper left')
        axBr = ax[1].twinx()
        axBr.plot(data_dict_deltaB_low['Epoch'][0], data_dict_deltaB_low['B_r'][0], color='tab:blue', label='B (+r)')
        axBr.set_ylabel('B (nT)')
        axBr.set_ylim(-6, 6)
        axBr.legend(loc='upper right')

        # Poynting Flux
        Sflux = ax[2].plot(data_dict_Poynting_low['Epoch'][0],data_dict_Poynting_low['S_p'][0])
        ax[2].set_ylabel('S_p (W/m)')
        ax[2].grid(True)

        # Density

        # filter the density
        DensityFilter = butter_filter(data_dict_densityInterp_low['ni'][0],
                                      lowcutoff=1,
                                      highcutoff=10,
                                      filtertype='Bandstop',
                                      order=4,
                                      fs=1000)

        ax[3].plot(data_dict_densityInterp_low['Epoch'][0], DensityFilter)
        ax[3].set_ylabel('$n_{i} Filtered$ ' + '$(cm^{-3})$')
        ax[3].set_xlabel('Time UTC on November 20, 2022')
        ax[3].grid(True)
        plt.show()


        ###### HIGH FLYER ######

        fig, ax = plt.subplots(nrows=4, ncols=1)
        plt.subplots_adjust(wspace=0, hspace=0)
        fig.suptitle('ACESII 36359')

        # B_e, E_r
        ax[0].grid(True)
        ax[0].plot(data_dict_deltaB_high['Epoch'][0], data_dict_deltaB_high['B_e'][0], color='tab:blue', label='B (+e)')
        ax[0].set_ylabel('B (nT)')
        ax[0].set_ylim(-15, 15)
        ax[0].legend(loc='upper left')

        # B_r, E_e
        ax[1].grid(True)
        ax[1].plot(data_dict_deltaB_high['Epoch'][0], data_dict_deltaB_high['B_r'][0], color='tab:blue', label='B (+r)')
        ax[1].set_ylabel('B (nT)')
        ax[1].set_ylim(-15, 15)
        ax[1].legend(loc='upper left')

        # Poynting Flux
        Sflux = ax[2].plot(data_dict_Poynting_high['Epoch'][0], data_dict_Poynting_high['S_p'][0])
        ax[2].set_ylabel('S_p (W/m)')
        ax[2].grid(True)

        # Density

        # filter the density
        DensityFilter = butter_filter(data_dict_densityInterp_high['ni'][0],
                                      lowcutoff=1,
                                      highcutoff=10,
                                      filtertype='Bandstop',
                                      order=4,
                                      fs=1000)

        ax[3].plot(data_dict_densityInterp_high['Epoch'][0], DensityFilter)
        ax[3].set_ylabel('$n_{i} Filtered$ ' + '$(cm^{-3})$')
        ax[3].set_xlabel('Time UTC on November 20, 2022')
        ax[3].grid(True)
        plt.show()

    ###################################
    # --- FFT deltaBdeltaE E/B PLOT ---
    ###################################

    if SECTION_deltaBdeltaEFFTspeed:

        prgMsg('Separating Data into regions')
        # separate the data into auroral, quiet region and alfvenic regions

        # --- Quiet region ---
        targetTimes_quiet = [dt.datetime(2022, 11, 20, 17, 24, 28, 00), dt.datetime(2022, 11, 20, 17, 24, 40, 00)]
        DataDicts_quiet = []
        for i in range(len(DataDicts)):
            tempDict= deepcopy(DataDicts[i])

            lowCut,highCut = np.abs(tempDict['Epoch'][0] - targetTimes_quiet[0]).argmin(),np.abs(tempDict['Epoch'][0] - targetTimes_quiet[1]).argmin()

            for key,val in tempDict.items():
                tempDict[key][0] = tempDict[key][0][lowCut:highCut]

            DataDicts_quiet.append(tempDict)

        # --- Alfvenic region ---
        targetTimes_alfven = [dt.datetime(2022, 11, 20, 17, 24, 30, 00), dt.datetime(2022, 11, 20, 17, 25, 34, 00)]
        DataDicts_Alfvenic = []
        for i in range(len(DataDicts)):
            tempDict = deepcopy(DataDicts[i])
            lowCut, highCut = np.abs(tempDict['Epoch'][0] - targetTimes_alfven[0]).argmin(), np.abs(tempDict['Epoch'][0] - targetTimes_alfven[1]).argmin()

            for key, val in tempDict.items():
                tempDict[key][0] = tempDict[key][0][lowCut:highCut]

            DataDicts_Alfvenic.append(tempDict)

        # ---auroral region ---
        targetTimes_auroral = [dt.datetime(2022,11,20,17,25,34,00),dt.datetime(2022,11,20,17,26,30,00)]

        Done(start_time)

        #######################
        # --- MAKE THE PLOT ---
        #######################

        prgMsg('Making Plot')

        wRegion = ['Quiet','Alfvenic']
        wDataDicts = [DataDicts_quiet,DataDicts_Alfvenic]

        ##### LOW FLYER ####
        fig, ax = plt.subplots(nrows=2, ncols=2)
        fig.suptitle('ACESII 36364')

        for j in range(2):

            # plt.subplots_adjust(wspace=0, hspace=0)
            Data = wDataDicts[j]
            deltaB = Data[2]
            deltaE = Data[0]
            Alfspeed = Data[3]

            # Br
            N, T = len(deltaB['Epoch'][0]), 1 / 128
            xf = fftfreq(N, T)[:N // 2]

            yf_Br = rfft(deltaB['B_r'][0])
            FFT_Br = 2.0 / N * np.abs(yf_Br[0:N // 2])

            yf_Be = rfft(deltaB['B_e'][0])
            FFT_Be = 2.0 / N * np.abs(yf_Be[0:N // 2])

            yf_Bp = rfft( deltaB['B_p'][0])
            FFT_Bp = 2.0 / N * np.abs(yf_Bp[0:N // 2])

            yf_Er = rfft(deltaE['E_r'][0])
            FFT_Er = 2.0 / N * np.abs(yf_Er[0:N // 2])

            yf_Ee = rfft(deltaE['E_e'][0])
            FFT_Ee = 2.0 / N * np.abs(yf_Ee[0:N // 2])

            yf_Ep = rfft(deltaE['E_p'][0])
            FFT_Ep = 2.0 / N * np.abs(yf_Ep[0:N // 2])

            # calculate the magnitudes
            FFT_Bmag = np.array([np.linalg.norm([FFT_Be[i], FFT_Bp[i], FFT_Br[i]]) for i in range(len(FFT_Br))])
            FFT_Emag = np.array([np.linalg.norm([FFT_Ee[i], FFT_Ep[i], FFT_Er[i]]) for i in range(len(FFT_Br))])

            # model Alfven Speed taken from "modelAlfvenSpeedAndResonance.py"
            modelAlfvenSpeed = [0,1947540.1772956662]

            # take the ratio
            EBratio = (1E6)*np.array([FFT_Emag[i]/FFT_Bmag[i] for i in range(len(FFT_Bmag))])

            # FFT mag Axis
            colorAx='tab:blue'
            ax[0, j].set_title(f'{wRegion[j]} Region ({targetTimes[0].strftime("%H:%M:%S")} to {targetTimes[1].strftime("%H:%M:%S")})')
            ax[0, j].plot(xf, FFT_Bmag, label='B',color=colorAx)
            ax[0, j].set_xlim(0, 10)
            ax[0, j].set_ylim(10E-5, 1E2)
            ax[0, j].set_yscale('log')
            ax[0, j].legend(loc='upper left')
            ax[0, j].set_ylabel('B [nT]')
            ax[0, j].spines['left'].set_color(colorAx)
            ax[0, j].yaxis.label.set_color(colorAx)
            ax[0, j].grid(visible=True, which='both', axis='both')
            # ax[0].axes.xaxis.set_visible(False)

            colorAx = 'tab:red'
            axEmag = ax[0, j].twinx()
            axEmag.plot(xf, FFT_Emag, label='E',color='tab:red')
            axEmag.set_yscale('log')
            axEmag.set_ylim(10E-5, 1E2)
            axEmag.set_ylabel('E [mV/m]')
            axEmag.spines['right'].set_color(colorAx)
            axEmag.yaxis.label.set_color(colorAx)
            axEmag.legend(loc='upper right')

            # E/B ratio axis

            # Find the average value of the alfven speed within this time

            AlfvenSpeedAvg = 10*sum(Alfspeed['AlfvenSpeed'][0])/len(Alfspeed['AlfvenSpeed'][0])
            ax[1, j].plot(xf, EBratio,label='E/B')
            ax[1, j].set_xlim(0, 10)
            ax[1, j].set_xlabel('Frequency [Hz]')

            ax[1, j].axhline(AlfvenSpeedAvg,color='red',linestyle='--')
            ax[1, j].text(x=8,y=(1.1)*AlfvenSpeedAvg,s='Alfven Speed Calc (Avg)')

            ax[1, j].axhline(modelAlfvenSpeed[1], color='tab:gray', linestyle='--')
            ax[1, j].text(x=8, y=(1.1) * modelAlfvenSpeed[1], s='Alfven Speed Model')

            # ax[1].set_ylim(1E4,1E8)
            ax[1, j].set_yscale('log')
            ax[1, j].set_ylabel('E/B [m/s]')
            ax[1, j].legend()
        plt.show()

        Done(start_time)

    ###########################################
    # --- WHICH FREQUENCIES ARE RESPONSIBLE ---
    ###########################################

    if SECTION_spectrogramData:

        #################
        # High Flyer PLOT
        #################
        fig = plt.figure()
        fig.suptitle('ACESII 36359')
        gs0 = fig.add_gridspec(nrows=4,ncols=2,width_ratios=[1,0.01])

        ax0 = fig.add_subplot(gs0[0, 0])
        ax1 = fig.add_subplot(gs0[1, 0])
        ax2 = fig.add_subplot(gs0[2, 0])
        ax3 = fig.add_subplot(gs0[3, 0])
        axCbar = fig.add_subplot(gs0[1:, 1])

        # Be,Br plot
        ax0.plot(data_dict_deltaB_high['Epoch'][0],data_dict_deltaB_high['B_e'][0],color='tab:blue',label='B_e')
        ax0.set_ylabel('B_e')
        ax0.legend(loc='upper left')

        axBr = ax0.twinx()
        axBr.plot(data_dict_deltaB_high['Epoch'][0], data_dict_deltaB_high['B_r'][0], color='tab:red',label='B_r')
        axBr.set_ylabel('B_r')
        axBr.legend(loc='upper right')

        # component spectrograms
        compNames = ['B_e','B_p','B_r']
        axes = [ax1,ax2,ax3]
        for i in range(3):
            # Be spectrogram
            data = data_dict_deltaB_high[compNames[i]][0]
            windowType, npersegN, scalingType = 'hann', 64, 'spectrum'  # spectrogram toggles
            overlap = int(npersegN * (7 / 8))  # hanning filter overlap
            f, t, Sxx = spectrogram(data, fs=128,
                                    window=windowType,
                                    nperseg=npersegN,  # note: if ==None default size is 256
                                    noverlap=overlap,
                                    scaling=scalingType)  # scaling = density or scaling = spectrum

            spectrogramPlot = axes[i].pcolormesh(t, f, Sxx, shading='nearest', vmin=1E-2, vmax=100, cmap='turbo',norm='log')

            axes[i].set_ylim(-0.1, 10)
            axes[i].set_ylabel('Frequency [Hz]')
            axes[i].set_xlabel('Time [Sec]')
            axes[i].axhline(spinFreq[0],color='white',linestyle='--',alpha=0.8)
            axes[i].text(0.5,1.4*spinFreq[0],'Spin Frequency',color='black')

        cbar = plt.colorbar(spectrogramPlot, cax=axCbar)
        cbar.set_label('Power Density [nT^2/Hz]')
        plt.show()

        #################
        # Low Flyer PLOT
        #################
        fig = plt.figure()
        fig.suptitle('ACESII 36364')
        gs0 = fig.add_gridspec(nrows=4, ncols=2, width_ratios=[1, 0.01])

        ax0 = fig.add_subplot(gs0[0, 0])
        ax1 = fig.add_subplot(gs0[1, 0])
        ax2 = fig.add_subplot(gs0[2, 0])
        ax3 = fig.add_subplot(gs0[3, 0])
        axCbar = fig.add_subplot(gs0[1:4, 1])

        # Be,Br plot
        ax0.plot(data_dict_deltaB_low['Epoch'][0], data_dict_deltaB_low['B_e'][0], color='tab:blue', label='B_e')
        ax0.set_ylabel('B_e')
        ax0.legend(loc='upper left')

        axBr = ax0.twinx()
        axBr.plot(data_dict_deltaB_low['Epoch'][0], data_dict_deltaB_low['B_r'][0], color='tab:red', label='B_r')
        axBr.set_ylabel('B_r')
        axBr.legend(loc='upper right')

        # component spectrograms
        compNames = ['B_e', 'B_p', 'B_r']
        axes = [ax1, ax2, ax3]
        for i in range(3):
            # Be spectrogram
            data = data_dict_deltaB_low[compNames[i]][0]
            windowType, npersegN, scalingType = 'hann', 64, 'density'  # spectrogram toggles
            overlap = int(npersegN * (7 / 8))  # hanning filter overlap
            f, t, Sxx = spectrogram(data, fs=128,
                                    window=windowType,
                                    nperseg=npersegN,  # note: if ==None default size is 256
                                    noverlap=overlap,
                                    scaling=scalingType)  # scaling = density or scaling = spectrum

            spectrogramPlot = axes[i].pcolormesh(t, f, Sxx, shading='nearest', vmin=1E-2, vmax=100, cmap='turbo', norm='log')
            axes[i].set_ylim(-0.1, 10)
            axes[i].set_ylabel('Frequency [Hz]')
            axes[i].set_xlabel('Time [Sec]')

            axes[i].axhline(spinFreq[0], color='white', linestyle='--', alpha=0.8)
            axes[i].text(0.5, 1.4 * spinFreq[0], 'Spin Frequency', color='Black')

        cbar = plt.colorbar(spectrogramPlot, cax=axCbar)
        cbar.set_label('Power Density [nT^2/Hz]')

        plt.show()




    #######################
    # --- WAVE GEOMETRY ---
    #######################

    # make a plot of the polarization of the waves
    # note that deltaB fluctuations should occur when deltaE also occur (hence electromagnetic)
    # and the correlation of the angle between E and B is high













# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
Prove_Alfvenic()

