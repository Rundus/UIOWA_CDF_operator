# # --- MAG_InteractiveFiltering_Despin.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: Script to handle the despin procedure in the magnetometer data.
# Step [1] Define a Region of interest and reduce data to only that window.
# Step [2] Using a functional form to fit spinning and coning, fit the 3 axes of data to get nominal spin and coning rate
# Step [3]: Reduce the data to a larger time window than step [1]; Since we are going to filter and
# perform singular spectral analysis, we should pad our dataset in time (>5 seconds) in order to
# prevent the unavoidable error that occurs near the data edges from seeping into our region of interest
# Step [4] Apply a reverse spin rotation (based on the rocket's spin rate) using a modelled DCM to eliminate a lot of the
# spin frequency power.
# Step [5] Using an ENU coordinate IGRF model, introduce a phase offset in step [4] that aligns the rocket's axes
# approximately to an inertial frame.
# Step [6] Apply a simple HighPass filter (~ 1 Hz) to knock out the fundamental coning an spin harmonics
# Step [7] Perform a Singular Spectral Analysis that picks out the harmonics of the spin frequencies, leaving only the
# "noise" part, which should be our signal

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

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_attitude = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'l2' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder


# --- --- --- SSA --- --- ---
SECTION_SSA = True
SSA_window_Size = 501
calculateSSA = True # calculate the SSA components and store them. THIS TOGGLE REQUIRES unSpinData and filterData both == True
###################
subSECTION_groupSSAData = True
wAxesSSA = 1 # 0 -> X, 1 -> Y, 2 -> Z
justPrintSSAFiles = False # TELLS YOU WHICH SSA FILES TO LOAD for the plotting
wSSAFile = 4  # select a specific SSA file to plot
reduceTimePercent = 5 # kill this percent of data on either end AFTER the SSA has been calculated
plotGroupedComponents = True
plotENUSpectrogram = False
plotwCorMatrix = False


# --- --- --- GENERAL --- --- ---
targetTimes = [dt.datetime(2022, 11, 20, 17, 24, 30, 000000), dt.datetime(2022, 11, 20, 17, 25, 30, 000000)] # ACTUAL WINDOW for WL 501
targetTimes = [dt.datetime(2022, 11, 20, 17, 24, 30, 000000), dt.datetime(2022, 11, 20, 17, 25, 30, 000000)] # ACTUAL WINDOW for WL 501
# --- --- --- FITFUNC --- --- ---
SECTION_fitFuncPlot = False
wAxes = 'By'
useFitResultsforP0 = False
# --- --- --- unSPIN --- --- ---
SECTION_unSpinData = False if not calculateSSA else True
useAttitudeDCM = True
plotIGRFcompare = False
# --- --- --- FILTERING --- --- ---
SECTION_filterData = False if not calculateSSA else True
plotFilteredAxes = False
lowCut_toggle, highcut_toggle, filttype_toggle, order_toggle = 1, 1.5, 'Highpass', 4 # filter toggles
windowType, npersegN, scalingType  = 'hann', 128, 'density' # spectrogram toggles
overlap = int(npersegN*(7/8)) # hanning filter overlap
# --- --- --- OUTPUT --- --- ---
outputData = True
# --- --- --- --- --- ---


# --- FIT RESULTS ---
fitResults = {
    'Bx': {'Spin Amp': 25.42873940404161, 'Spin Freq': 0.6463295881639182, 'Spin Phase': 91.9759995936283, 'Cone Amp': 625.8772357084948, 'Cone Freq': 0.05294818121871208, 'Cone Phase': -138.77308595997619, 'Offset': -44919.748937299344},
    'By': {'Spin Amp': 7.378420193701481, 'Spin Freq': 0.6442248190622027, 'Spin Phase': 109.20255873087793, 'Cone Amp': 1380.5616077430786, 'Cone Freq': 0.02700105226961604, 'Cone Phase': 109.87799606103452, 'Offset': -139.74554466082876},
    'Bz': {'Spin Amp': 8.095746809541962, 'Spin Freq': 0.6442537451458561, 'Spin Phase': 19.11852573798773, 'Cone Amp': 1257.0313161879794, 'Cone Freq': 0.026874206798816504, 'Cone Phase': -69.78175516947503, 'Offset': 32.456720919269245}
}

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import pandas as pd
from pyIGRF import igrf_value
from ACESII_code.myImports import *
from matplotlib.widgets import Slider
from ACESII_code.class_var_func import butter_filter
from numpy.fft import rfft, fftfreq
from ACESII_code.class_var_func import DCM
from scipy.interpolate import CubicSpline
from ACESII_code.supportCode.Support_Libraries.pymssa import MSSA
from scipy.signal import spectrogram



def RingCore_L1_to_L2_Despin(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L0_ACES_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*RingCore_rktFrm*')
    inputFiles_attitude = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[wflyer]}{modifier}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]

    fileoutName = f'ACESII_{rocketID}_{outputPath_modifier.lower()}_RingCore_DeSpun'


    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'DeSpining RingCore Data' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the Magnetometer file ---
        prgMsg(f'Loading data from {inputPath_modifier} RingCore Files')
        data_dict_mag = loadDictFromFile(inputFiles[wFile],{})
        Done(start_time)

        # --- get the data from the Magnetometer file ---
        prgMsg(f'Loading data from {inputPath_modifier_attitude} Files')
        data_dict_attitude = loadDictFromFile(inputFiles_attitude[0], {})
        Done(start_time)

        ########################
        # --- Reduce dataset ---
        ########################
        prgMsg('Reducing Dataset')

        # --- apply reduction to mag data ---
        lowCutoff, highCutoff = np.abs(data_dict_mag['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_mag['Epoch'][0] - targetTimes[1]).argmin()
        for key, val in data_dict_mag.items():
            data_dict_mag[key][0] = np.array(data_dict_mag[key][0][lowCutoff:highCutoff])

        data_dict_mag['Epoch'][0] = np.array([ pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag['Epoch'][0]])
        Epoch_seconds = np.array([(tme - data_dict_mag['Epoch'][0][0]) / 1E9 for tme in data_dict_mag['Epoch'][0]])
        Epoch_dt = np.array([ pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_mag['Epoch'][0]])

        # --- apply reduction to attitude data ---
        lowCutoff, highCutoff = np.abs(data_dict_attitude['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_attitude['Epoch'][0] - targetTimes[1]).argmin()
        for key, val in data_dict_attitude.items():
            data_dict_attitude[key][0] = np.array(data_dict_attitude[key][0][lowCutoff:highCutoff])

        data_dict_attitude['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude['Epoch'][0]])
        Done(start_time)

        # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
        # --- interpolate attitude data up to magnetometer epoch ---
        # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

        # results from RinCore_L1_DetermineAdjustments_Before_Despin that suggest a linear offset of the attitude data in time
        # offsetResults_intercept = [104944444.44444445, 117000000]
        offsetResults_intercept = [0, 117000000]
        # offsetResults_VectorScalar = np.array([[0, 0, -234.47368421052633],[0,0,-355.47368421052633]])
        offsetResults_VectorScalar = np.array([[0, 0, 0], [0, 0, 0]])
        Epoch_attitude_loop = np.array([int(tme + offsetResults_intercept[wRocket - 4]) for tme in data_dict_attitude['Epoch'][0]])

        dataKeys = ['Epoch', 'Alt', 'Latgd', 'Long', 'Y_Az', 'a11', 'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33']
        dataKeysVal = [deepcopy(data_dict_mag['Epoch'][0]), [], [], [], [], [], [], [], [], [], [], [], [], []]
        attitudeData = [deepcopy(data_dict_attitude[key][0]) for key in dataKeys]  # a list to contain the attitude only the data that I care about
        dataInterp_dict_attitude = {key: value for key, value in zip(dataKeys, dataKeysVal)}

        counter = 0
        for key, newDataList in dataInterp_dict_attitude.items():
            if key != 'Epoch':
                # --- cubic interpolation ---
                splCub = CubicSpline(Epoch_attitude_loop, attitudeData[counter])

                # evaluate the interpolation at all the epoch_mag points
                dataInterp_dict_attitude[key] = np.array([splCub(timeVal) for timeVal in dataInterp_dict_attitude['Epoch']])
            counter += 1

        # FORMAT: [Bx, By, Bz]
        # coneFreq =  sum([0.05294816949860347/2,0.02700105226961604,0.026874206799143967 ])/3
        coneFreq = sum([fitResults['By']['Cone Freq'], fitResults['Bz']['Cone Freq']]) / 2
        spinFreq = sum([fitResults['Bz']['Spin Freq'], fitResults['By']['Spin Freq'], fitResults['Bz']['Spin Freq']]) / 3

        if SECTION_fitFuncPlot:
            # --- define the fit function ---
            if wAxes in ['By', 'Bz']:
                def fitFunc(x, spinAmp, spinFreq, spinPhase, coneAmp, coneFreq, conePhase, offset):
                    return (spinAmp * np.sin(2 * np.pi * spinFreq * x + np.radians(spinPhase)) *
                            coneAmp * np.cos(2 * np.pi * coneFreq * x + np.radians(conePhase)) +
                            offset)
            elif wAxes in ['Bx']:
                # def fitFunc(x, coneAmp, coneFreq, conePhase, offset):
                #     return coneAmp * np.cos(2*np.pi*coneFreq * x + np.radians(conePhase)) + offset

                def fitFunc(x, spinAmp, spinFreq, spinPhase, coneAmp, coneFreq, conePhase, offset):
                    return (spinAmp * np.sin(2 * np.pi * spinFreq * x + np.radians(spinPhase)) +
                            coneAmp * np.cos(2 * np.pi * coneFreq * x + np.radians(conePhase)) +
                            offset)

                # fitNames = ['Cone Amp', 'Cone Freq', 'Cone Phase', 'Offset']
                # p0 = [600, 0.05, -97, -46000]

            if useFitResultsforP0:
                fitNames, p0 = [], []
                for key, val in fitResults[wAxes].items():
                    fitNames.append(key)
                    p0.append(val)
            else:
                fitNames = ['Spin Amp', 'Spin Freq', 'Spin Phase', 'Cone Amp', 'Cone Freq', 'Cone Phase', 'Offset']
                p0 = [8.26770739, 0.6442537451403045, 90, 1300, 0.027, 0, 0]

            # --- organize the data ---
            xData = Epoch_seconds
            yData = np.array(data_dict_mag[wAxes][0])

            # --- Fit Curve to Raw Data ---

            #    spinAmp, spinFreq, spinPhase, coneAmp, coneFreq, conePhase, offset
            # bounds = ((1,      0.55,   -360,       1000,     0.01,      -360,  -1000),
            #          (12,      0.75,    360,       1300,     0.05,       360,  1000))
            fit_params, cov = scipy.optimize.curve_fit(fitFunc, xData, yData, p0=p0, maxfev=100000)

            yData_fit = np.array([fitFunc(val, *fit_params) for val in xData])
            sub = yData_fit - yData

            for key, val in zip(fitNames, fit_params):
                print(key, ':', val)

            print({key: val for key, val in zip(fitNames, fit_params)})

            # --- Initial filter of the data ---
            yData_filtered = butter_filter(sub, lowcutoff=0.5, highcutoff=1, order=1, filtertype='bandstop', fs=128)

            # --- FFT of the Data ---
            N, T = len(sub), 1 / 128
            yf_filt, xf = rfft(yData_filtered), fftfreq(N, T)[:N // 2]
            yf_sub = rfft(sub)

            ###############################
            # --- Plot the Initial Data ---
            ###############################
            fig, ax = plt.subplots(4)
            ax[0].plot(xData, yData)
            ax[0].plot(xData, yData_fit)
            ax[1].plot(xData, sub)
            filteredData, = ax[2].plot(xData, yData_filtered)
            ax[3].plot(xf, 2.0 / N * np.abs(yf_sub[0:N // 2]), color='red', alpha=0.4, linestyle='--',
                       label='Subtraction FFT unfiltered')
            FFT_filtered_plot, = ax[3].plot(xf, 2.0 / N * np.abs(yf_filt[0:N // 2]), label='Subtraction FFT Filtered')

            # --- Adjustments ---
            fig.subplots_adjust(left=0.25, bottom=0.25)  # adjust the main plot to make room for the sliders
            ax[0].set_ylabel(wAxes + ' [nT]')
            ax[1].set_ylabel(wAxes + '- $B_{fit}$')
            ax[2].set_ylabel('Filtered Subtraction')
            ax[3].set_ylabel('FFT Power')
            ax[3].set_xlim(-0.1, 5)
            fig.suptitle(f'{wAxes}\n Type: Band Stop\n')
            plt.legend()

            #################
            # --- SLIDERS ---
            #################
            axfilter_cutoff_low = fig.add_axes([0.05, 0.25, 0.0225, 0.63])
            slider_cutoff_low = Slider(ax=axfilter_cutoff_low, label='cutoff_low', valmin=0.0001, valmax=8, valinit=1, orientation="vertical")

            axfilter_cutoff_high = fig.add_axes([0.15, 0.25, 0.0225, 0.63])
            slider_cutoff_high = Slider(ax=axfilter_cutoff_high, label='cutoff_high', valmin=0.001, valmax=63, valinit=1.5, orientation="vertical")

            axfilter_order = fig.add_axes([0.25, 0.1, 0.65, 0.03])
            slider_filter_order = Slider(ax=axfilter_order, label='order', valmin=1, valmax=15, valinit=1)

            # --- SLIDER update function ---
            def f(cutoff_low, cutoff_high, order):
                sub = np.array([fitFunc(val, *fit_params) for val in xData]) - yData
                updated_filteredData = butter_filter(sub, lowcutoff=cutoff_low, highcutoff=cutoff_high,
                                                     order=int(order),
                                                     filtertype='bandpass', fs=128)
                yf = rfft(updated_filteredData)
                return updated_filteredData, yf

            def update(val):
                # calculate the newly filtered data and its FFT
                newData, yf_new = f(slider_cutoff_low.val, slider_cutoff_high.val, slider_filter_order.val)

                # update the newly filtered data
                filteredData.set_ydata(newData)
                FFT_filtered_plot.set_ydata(2.0 / N * np.abs(yf_new[0:N // 2]))

                # adjust the y-scale of filtered data
                ax[2].set_ylim(min(newData), max(newData))

                # update canvas
                fig.canvas.draw_idle()

            # --- register the update function with each slider ---
            slider_filter_order.on_changed(update)
            slider_cutoff_high.on_changed(update)
            slider_cutoff_low.on_changed(update)
            plt.show()

        if SECTION_unSpinData:


            # --- --- --- --- --- --- --- --- --- -
            # --- REVERSE ROTATE TO REMOVE SPIN ---
            # --- --- --- --- --- --- --- --- --- -

            # define the Yaw,Pitch,Roll angles to use over the timeseries. USE DEGREES since DCM does a radian conversion
            B_rkt = np.array([[data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]] for i in range(len(Epoch_seconds))])

            if useAttitudeDCM:
                # form the DCM matrix

                DCMmat = np.array([
                    [[dataInterp_dict_attitude['a11'][i], dataInterp_dict_attitude['a12'][i], dataInterp_dict_attitude['a13'][i]],
                     [dataInterp_dict_attitude['a21'][i], dataInterp_dict_attitude['a22'][i], dataInterp_dict_attitude['a23'][i]],
                     [dataInterp_dict_attitude['a31'][i], dataInterp_dict_attitude['a32'][i], dataInterp_dict_attitude['a33'][i]]]
                    for i in range(len(Epoch_seconds))
                ])


                B_noSpin = np.array([np.matmul(DCMmat[i], B_rkt[i]) for i in range(len(Epoch_seconds))]) + offsetResults_VectorScalar[wRocket-4]
                mapping = [0, 1, 2]
            else:

                # My choice of mapping:
                # B_rkt_Y --> East
                # B_rkt_Z --> North
                # B_rkt_X --> Up
                mapping = [1, 2, 0]
                initPhase = -1 * data_dict_attitude['Y_Az'][0][0] - 90  # 152.44deg to rotate Y axis aligned to North (from attitude Y_Az at t=0), then -90 to rotate to east

                def getB_noSpin(initialPhase, freq):
                    # get my DCM over time
                    Rolls = np.array([-1*np.degrees(2*np.pi*(freq) * tme) + initialPhase for tme in Epoch_seconds])
                    Yaws = np.array([0 for tme in Epoch_seconds])
                    Pitchs = np.array([0 for tme in Epoch_seconds])
                    DCMmat = np.array([DCM(Rolls[i], Pitchs[i], Yaws[i]) for i in range(len(Epoch_seconds))])

                    # Apply DCM to my data to un-spin it
                    return np.array([np.matmul(DCMmat[i], B_rkt[i]) for i in range(len(Epoch_seconds))])

                # Apply DCM to my data to un-spin it
                B_noSpin = getB_noSpin(initialPhase=initPhase, freq=spinFreq)

            # --- --- --- --- --- --- --- --- --- --- --
            # --- ADJUST INITIAL PHASE TO MATCH IGRF ---
            # --- --- --- --- --- --- --- --- --- --- --

            # --- get IGRF ENU ---
            date = 2022 + 323 / 365  # Corresponds to 11/20/2022
            ### IGRF info ###
            # [3] North Comp (+ N | - S)
            # [4] East Comp (+ E | - W)
            # [5] Vertical Comp (+ D | - U)
            # [6] Total Field

            IGRF = np.array([igrf_value(dataInterp_dict_attitude['Latgd'][i], dataInterp_dict_attitude['Long'][i],dataInterp_dict_attitude['Alt'][i]/1000, date) for i in range(len(dataInterp_dict_attitude['Epoch']))])
            IGRF_ENU = np.array([[vec[4], vec[3], -1 * vec[5]] for vec in IGRF])

            if plotIGRFcompare:
                # compare ENU IGRF to my unspun rocket data
                comps_IGRF = ['B_East', 'B_North', 'B_Up']
                fig, ax = plt.subplots(3)

                if useAttitudeDCM:
                    fig.suptitle('Attitude DCM')
                else:
                    fig.suptitle(f'Initial Phase: {initPhase}'+'$^{\circ}$')

                # East (B_rkt_Y)
                B_rkt_east_plot, = ax[0].plot(Epoch_seconds, B_noSpin[:, mapping[0]])
                ax[0].plot(Epoch_seconds, IGRF_ENU[:, 0],label='IGRF East')
                ax[0].set_ylabel(f'Rkt_X [nT]')
                # ax[0].set_ylim(-10000, 10000)

                # North (B_rkt_Z)
                i = 1
                B_rkt_north_plot, = ax[1].plot(Epoch_seconds, B_noSpin[:, mapping[1]])
                ax[1].plot(Epoch_seconds, IGRF_ENU[:, 1], label='IGRF North')
                ax[1].set_ylabel(f'Rkt_Y [nT]')
                # ax[1].set_ylim(-10000, 10000)

                # Up (B_rkt_X)
                B_rkt_up_plot, = ax[2].plot(Epoch_seconds, B_noSpin[:, mapping[2]])
                ax[2].plot(Epoch_seconds, IGRF_ENU[:, 2], label='IGRF Up')
                ax[2].set_ylabel(f'Rkt_Z [nT]')
                ax[2].set_xlabel('Seconds From 17:24:00')

                # --- --- --- --
                # --- SLIDER ---
                # --- --- --- --

                fig.subplots_adjust(left=0.25, bottom=0.25)

                # Make a horizontal slider to control the initial phae.
                axPhase = fig.add_axes([0.25, 0.1, 0.65, 0.03])
                phase_slider = Slider(ax=axPhase, label='Init Phase', valmin=-360, valmax=360, valinit=0)

                # Make a horizontal slider to control the frequency.
                axFreq = fig.add_axes([0.25, 0.15, 0.65, 0.03, ])
                freq_slider = Slider(ax=axFreq, label='Frequency', valmin=0.0001, valmax=1, valinit=0)

                def update(val):
                    newData = getB_noSpin(phase_slider.val, freq_slider.val)
                    B_rkt_east_plot.set_ydata(newData[:, 1])
                    B_rkt_north_plot.set_ydata(newData[:, 2])
                    B_rkt_up_plot.set_ydata(newData[:, 0])
                    fig.canvas.draw_idle()

                # register the update function with each slider
                phase_slider.on_changed(update)
                freq_slider.on_changed(update)

                plt.show()

            data_for_output = B_noSpin

        if SECTION_filterData:
            # --- --- --- --- -
            # --- FILTERING ---
            # --- --- --- --- -
            prgMsg('Filtering Data')

            # Apply Highpass filter to data
            B_noSpin = data_for_output # get the data
            B_rkt_filtered = []
            for i in range(3):
                B_rkt_filtered.append(butter_filter(B_noSpin[:, i], lowcutoff=lowCut_toggle, highcutoff=highcut_toggle, filtertype=filttype_toggle, order=order_toggle, fs=128))

            B_rkt_filtered = np.array(B_rkt_filtered)

            comps = ['B_east', 'B_north', 'B_up'] if useAttitudeDCM else ['Bx_noSpin', 'By_noSpin', 'Bz_noSpin']

            if plotFilteredAxes:

                for i in range(3):
                    ###############################
                    # --- Plot the Initial Data ---
                    ###############################
                    fig, ax = plt.subplots(nrows=4,ncols=1,constrained_layout=True)

                    ax[0].plot(Epoch_dt, B_noSpin[:, i], label='noSpin')
                    ax[0].set_ylabel(f'{comps[i]}')

                    # --- FFT noSpin ---
                    N, T = len(B_noSpin[:, i]), 1 / 128
                    yf_rawData = rfft(B_noSpin[:, i])

                    # --- Highpass Filter ---
                    filteredData = butter_filter(B_noSpin[:, i], lowcutoff=lowCut_toggle, highcutoff=highcut_toggle, filtertype=filttype_toggle, order=order_toggle, fs=128)
                    filteredDataPlot, = ax[1].plot(Epoch_dt, filteredData, color='red')
                    ax[1].set_ylabel(f'{comps[i]}_filtered')
                    ax[1].set_xlabel('Time [s]')

                    # --- FFT filtered ---
                    N, T = len(B_noSpin[:, i]), 1 / 128
                    yf_filtered = rfft(filteredData)
                    xf = fftfreq(N, T)[:N // 2]
                    FFT_filtered_plot, = ax[2].plot(xf, 2.0 / N * np.abs(yf_filtered[0:N // 2]))
                    ax[2].plot(xf, 2.0 / N * np.abs(yf_rawData[0:N // 2]), color='orange')
                    ax[2].set_ylabel('FFT Power')
                    ax[2].set_xlabel('Frequency [Hz]')
                    ax[2].set_xlim(-0.1, 5)
                    ax[2].set_ylim(-0.1, 5)

                    # --- PERIODOGRAM ---
                    f, t, Sxx = spectrogram(filteredData, fs=128,
                                window=windowType,
                                nperseg=npersegN, # note: if ==None default size is 256
                                noverlap=overlap,
                                scaling=scalingType) # scaling = density or scaling = spectrum
                    spectrogramPlot = ax[3].pcolormesh(t, f, Sxx, shading='nearest',vmin=0,vmax=1,cmap='turbo')
                    cbar = plt.colorbar(spectrogramPlot,ax=ax[3])
                    ax[3].set_ylim(-0.1, 15)
                    ax[3].set_ylabel('Frequency [Hz]')
                    ax[3].set_xlabel('Time [Sec]')

                    fig.subplots_adjust(left=0.1, bottom=0.1)
                    fig.suptitle(f'{comps[i]}\n'
                                 f'Type: {filttype_toggle}')

                    #################
                    # --- SLIDERS ---
                    #################
                    # fig.add_axes([x-pos,y-pos,width,height] (x,y) of bottom left corner
                    axfilter_cutoff_low = fig.add_axes([0.92, 0.3, 0.02, 0.63])
                    slider_cutoff_low = Slider(ax=axfilter_cutoff_low, label='lowFq', valmin=0.0001, valmax=8, valinit=1, orientation="vertical")

                    axfilter_cutoff_high = fig.add_axes([0.97, 0.3, 0.02, 0.63])
                    slider_cutoff_high = Slider(ax=axfilter_cutoff_high, label='highFq', valmin=0.001, valmax=63, valinit=1.5, orientation="vertical")

                    axfilter_order = fig.add_axes([0.945, 0.3, 0.02, 0.56])
                    slider_filter_order = Slider(ax=axfilter_order, label='order', valmin=1, valmax=15, valinit=1, orientation="vertical")

                    def f(cutoff_low, cutoff_high, order):
                        updated_filteredData = butter_filter(B_noSpin[:, i], lowcutoff=cutoff_low, highcutoff=cutoff_high, order=int(order), filtertype=filttype_toggle, fs=128)
                        yf = rfft(updated_filteredData)
                        f, t, Sxx = spectrogram(updated_filteredData, fs=128,
                                                window=windowType,
                                                nperseg=npersegN,  # note: if ==None default size is 256
                                                noverlap=overlap,
                                                scaling=scalingType)  # scaling = density or scaling = spectrum
                        return updated_filteredData, yf,f,t,Sxx

                    def update(val):
                        # calculate the newly filtered data and its FFT
                        newData, yf_new,f_new,t_new,Sxx_new = f(slider_cutoff_low.val, slider_cutoff_high.val, slider_filter_order.val)

                        # update the newly filtered data
                        filteredDataPlot.set_ydata(newData)
                        FFT_filtered_plot.set_ydata(2.0 / N * np.abs(yf_new[0:N // 2]))

                        # adjust the y-scale of filtered data
                        ax[1].set_ylim(min(newData), max(newData))

                        # update the spectrogram
                        spectrogramPlot.set_array(Sxx_new)

                        # update canvas
                        fig.canvas.draw_idle()

                    slider_filter_order.on_changed(update)
                    slider_cutoff_high.on_changed(update)
                    slider_cutoff_low.on_changed(update)
                    plt.show()

            # format data for output
            data_for_output = np.array([ [B_rkt_filtered[0][i],B_rkt_filtered[1][i],B_rkt_filtered[2][i]] for i in range(len(B_rkt_filtered[0]))])
            Done(start_time)

        if SECTION_SSA:

            # output file location for MSSA
            outputPathSSA = f'{rocketFolderPath}\\science\despinSSAcomponents\\{fliers[wflyer]}\\{fileoutName}_SSAcomponents_WL{SSA_window_Size}.cdf'

            # name of the components
            compNames = ['B_east_SSA', 'B_north_SSA', 'B_up_SSAcomps'] if useAttitudeDCM else ['Bx_SSAcomps', 'By_SSAcomps', 'Bz_SSAcomps']

            # create the MSSA object
            mssa = MSSA(n_components=None, window_size=SSA_window_Size, verbose=False)

            # --- perform the mSSA or SSA on specific components---
            if calculateSSA:

                prgMsg('Calculating SSA components')

                # convert data to pandas dataframe
                dataFormatted = {
                    'Bx': data_for_output[:, 0],
                    'By': data_for_output[:, 1],
                    'Bz': data_for_output[:, 2]}

                data = pd.DataFrame(dataFormatted)

                # calculate the mSSA
                mssa.fit(data)

                # get the mSSA components
                components = mssa.components_

                # --- output the data ---
                example_attrs = {'LABLAXIS': None, 'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None,
                         'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': None,
                         'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}

                data_dict_SSAcomps = {}

                for i in range(3):
                    dataToOutput = np.array(components[i, :, :])
                    attrs = deepcopy(example_attrs)
                    attrs['LABLAXIS'] = compNames[i]
                    attrs['VALIDMIN'] = dataToOutput.min()
                    attrs['VALIDMAX'] = dataToOutput.max()
                    data_dict_SSAcomps = {**data_dict_SSAcomps, **{compNames[i]:[dataToOutput, attrs]}}

                outputCDFdata(outputPathSSA, data_dict_SSAcomps, outputModelData, globalAttrsMod, 'RingCore')
                Done(start_time)
            elif subSECTION_groupSSAData:

                # load components data from file
                SSAFiles = glob(f'{rocketFolderPath}\\science\despinSSAcomponents\\{fliers[wflyer]}\\*.cdf*')

                if justPrintSSAFiles:
                    ssa_names = [ssafile.replace(f'{rocketFolderPath}\\science\despinSSAcomponents\\{fliers[wflyer]}\\', '') for ssafile in SSAFiles]

                    for i, file in enumerate(ssa_names):
                        print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, ssa_names[i], round(getsize(SSAFiles[i]) / (10 ** 6), 1)))

                else:
                    # find the windowsize in the file I've chosen
                    windowSize = int(''.join(x for x in SSAFiles[wSSAFile].replace('36359','').replace('36364','').replace('l2','') if x.isdigit()))
                    data_dict_SSA = loadDictFromFile(SSAFiles[wSSAFile],{})

                    prgMsg('Grouping mSSA elements')
                    from ACESII_code.Processing.SSAgrouping import groupings
                    groupings = groupings(wRocket=wRocket, SSA_window_Size=windowSize)

                    # get all the SSA components for the three axes
                    B_SSA = [data_dict_SSA[compNames[0]][0], data_dict_SSA[compNames[1]][0], data_dict_SSA[compNames[2]][0]]

                    # --- Plot the FFT and groupings for one wAxes axes ---
                    if plotGroupedComponents:

                        # --- Plot it ---
                        fig, ax = plt.subplots(nrows=len(groupings),ncols=2)
                        fig.suptitle(compNames[wAxesSSA] + f'\n Window Length: {windowSize}')

                        # loop over all the groups in grouping
                        for i in range(len(groupings)):
                            data = np.array(B_SSA[wAxesSSA])

                            # combine the groupings
                            plotThisData = np.array([0 for i in range(len(data[:,0]))],dtype='float64')
                            for j, componentIndex in enumerate(groupings[i]):
                                plotThisData += np.array(data[:, componentIndex])

                            # reduce the last "X" percent of the data on either end to eliminate the SSA effect
                            percent = reduceTimePercent/100
                            percentLow, percentHigh = int((percent)*len(plotThisData)), int((1-percent)*len(plotThisData))
                            plotThisData = plotThisData[percentLow:percentHigh]

                            Epoch_dt = np.array([pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_mag['Epoch'][0]])
                            Epoch_dt = Epoch_dt[percentLow:percentHigh]

                            # plot the component
                            ax[i, 0].plot(Epoch_dt, plotThisData)

                            if i == 0:
                                ax[i, 0].set_ylabel('Original')
                            elif i == len(groupings)-1:
                                ax[i, 0].set_ylabel('Physical Signal')
                                ax[i, 1].set_xlabel('Frequency [Hz]')
                                ax[i, 1].set_ylabel('FFT')
                            else:
                                ax[i, 0].set_ylabel('F{}'.format(i))

                            # calculate the FFT and plot it
                            N, T = len(plotThisData), 1 / 128
                            yf, xf = rfft(plotThisData), fftfreq(N, T)[:N // 2]
                            FFT = 2.0 / N * np.abs(yf[0:N // 2])
                            ax[i, 1].plot(xf, FFT)
                            ax[i, 1].vlines([spinFreq*(i+1) for i in range(50)],ymin=min(FFT), ymax=max(FFT),alpha=0.5, color='red')
                            ax[i, 1].set_xlim(-0.1, 10)
                            ax[i, 1].set_ylim(0, max(FFT))

                        plt.show()

                    # --- collect the group info for the LAST (noise) grouping ---
                    if plotENUSpectrogram:

                        # --- group the three components---
                        groupedData = []

                        for i in range(3):# loop over all three axes

                            data = np.array(B_SSA[i]) # get the SSA data for this axes

                            # --- collect the group info for the LAST (noise) grouping ---
                            plotThisData = np.array([0 for i in range(len(data[:, 0]))], dtype='float64')
                            for componentIndex in groupings[-1]:
                                plotThisData += np.array(data[:, componentIndex])

                            # reduce the last "X" percent of the data on either end to eliminate the SSA effect
                            percent = reduceTimePercent / 100
                            percentLow, percentHigh = int((percent) * len(plotThisData)), int((1 - percent) * len(plotThisData))
                            Epoch_dt = np.array([pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_mag['Epoch'][0]])
                            Epoch_dt = Epoch_dt[percentLow:percentHigh]

                            # append "noise" data to "groupedData"
                            groupedData.append(plotThisData[percentLow:percentHigh])

                        # --- Plot the Data ---
                        fig, ax = plt.subplots(nrows=3, ncols=2, constrained_layout=True)
                        fig.suptitle(f'ACES II {rocketID}\n'
                                     f'Window Length: {windowSize}')

                        for i in range(3):
                            # plot the SSA'd filtered data
                            ax[i, 0].plot(Epoch_dt, groupedData[i])
                            ax[i, 0].set_ylabel(compNames[i] + ' [nT]')

                            # plot the spectrogram of the data
                            f, t, Sxx = spectrogram(groupedData[i], fs=128,
                                                    window=windowType,
                                                    nperseg=npersegN,  # note: if ==None default size is 256
                                                    noverlap=overlap,
                                                    scaling=scalingType)  # scaling = density or scaling = spectrum

                            # convert spectrogram time to Epoch
                            spectroEpoch = np.array([pycdf.lib.tt2000_to_datetime(pycdf.lib.datetime_to_tt2000(Epoch_dt[0]) + int(1E9*tme)) for tme in t])

                            spectrogramPlot = ax[i, 1].pcolormesh(spectroEpoch, f, Sxx, shading='nearest', vmin=-0.1, vmax=10, cmap='turbo')
                            cbar = plt.colorbar(spectrogramPlot, ax=ax[i, 1])
                            cbar.set_label('$B^{2}$/Hz')
                            ax[i, 1].set_ylim(-0.1, 15)
                            ax[i, 1].set_ylabel('Frequency [Hz]')
                            ax[i, 1].set_xlabel('Time [Sec]')



                        plt.show()

                    if plotwCorMatrix:

                        # plot all three axes correlation matrix
                        mssa.fit(pd.DataFrame(
                            {'Data1': [i for i in range(len(Epoch_seconds))],
                             'Data2': [i for i in range(len(Epoch_seconds))],
                             'Data3': [i for i in range(len(Epoch_seconds))]}
                        ))

                        for i in range(3):
                            mssa.components_[i, :, :] = data_dict_SSA[compNames[i]][0]

                        for i in range(3):
                            # calculate correlation matrix
                            wcorr = np.abs(mssa.w_correlation(mssa.components_[i, :, :]))

                            # plot it
                            plt.title(compNames[i])
                            ax = plt.imshow(wcorr, cmap='turbo')
                            plt.xlabel(r"$\tilde{F}_i$")
                            plt.ylabel(r"$\tilde{F}_j$")
                            plt.colorbar(ax.colorbar, fraction=0.045)
                            ax.colorbar.set_label("$W_{i,j}$")
                            plt.clim(0, 1)
                            plt.show()


                    # --- format the data for output ---
                    data_for_output = []
                    newComps = ['dB_east', 'dB_north', 'dB_up']
                    for k in range(len(newComps)): # loop through all the components, but only take the last grouping

                        data = np.array(B_SSA[k])

                        # combine the groupings of the last group set only
                        formattedData = np.array([0 for i in range(len(data[:, 0]))], dtype='float64')
                        for j, componentIndex in enumerate(groupings[-1]):
                            formattedData += np.array(data[:, componentIndex])

                        # reduce the last "X" percent of the data on either end to eliminate the SSA effect
                        percent = reduceTimePercent / 100
                        percentLow, percentHigh = int((percent) * len(formattedData)), int((1 - percent) * len(formattedData))
                        data_for_output.append(formattedData[percentLow:percentHigh])

                    Epoch_SSA = data_dict_mag['Epoch'][0][percentLow:percentHigh]

                    # calculate dBmag
                    data_for_output.append(np.array([np.linalg.norm([data_for_output[0][i],data_for_output[1][i],data_for_output[2][i]]) for i in range(len(data_for_output[0]))]))

                    # reformat the data for output
                    data_for_output = np.array([  [data_for_output[0][i],data_for_output[1][i],data_for_output[2][i],data_for_output[3][i]] for i in range(len(data_for_output[0])) ])

                Done(start_time)



        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('Creating output file')

            # create the output data_dict
            data_dict = deepcopy(data_dict_mag)
            comps = ['Bx', 'By', 'Bz', 'Bmag']
            newComps = ['dB_east', 'dB_north', 'dB_up','dBmag']

            # --- Magnetic Components ---
            # get the attributes of the old components and replace them
            for i, key in enumerate(comps):
                newAttrs = deepcopy(data_dict[key][1])
                newAttrs['LABLAXIS'] = newComps[i]

                # remove the old key
                del data_dict[key]

                # append the new key
                data_dict = {**data_dict, **{newComps[i] : [data_for_output[:, i], newAttrs]}}

            if subSECTION_groupSSAData:
                data_dict['Epoch'][0] = Epoch_SSA

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}.cdf'

            outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, 'RingCore')

            Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        RingCore_L1_to_L2_Despin(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            RingCore_L1_to_L2_Despin(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            RingCore_L1_to_L2_Despin(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)