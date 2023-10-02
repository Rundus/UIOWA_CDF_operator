# # --- EFI_L2_to_L2_dE.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION:

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

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
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'l2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_attitude = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\deltaE' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder


# --- --- --- Time region --- --- ---
useAlfvenRegion = True

if useAlfvenRegion:
    scienceRegions = [
        [dt.datetime(2022, 11, 20, 17, 24, 25, 000000), dt.datetime(2022, 11, 20, 17, 25, 18, 000000)],
        [dt.datetime(2022, 11, 20, 17, 24, 25, 000000), dt.datetime(2022, 11, 20, 17, 25, 18, 000000)]]# ACTUAL WINDOW for WL 501 High Flyer
else:
    scienceRegions = [
        [dt.datetime(2022, 11, 20, 17, 24, 25, 000000), dt.datetime(2022, 11, 20, 17, 27, 50, 000000)],
        [dt.datetime(2022, 11, 20, 17, 24, 00, 000000), dt.datetime(2022, 11, 20, 17, 27, 50, 000000)]
    ]

EFIsampleFreq = 4000
convertTomVpm = True # converts from V/m to mV/m

# --- --- --- SSA --- --- ---
SECTION_SSA = True
SSA_window_Size = 501
calculateSSA = False # calculate the SSA components and store them. THIS TOGGLE REQUIRES unSpinData and filterData both == True
###################
subSECTION_groupSSAData = True if not calculateSSA else False
wAxesSSA = 0 # 0 -> X, 1 -> Y, 2 -> Z
justPrintSSAFiles = False # TELLS YOU WHICH SSA FILES TO LOAD for the plotting
wSSAFile = 0  # select a specific SSA file to plot
reduceTimePercent = 2 # kill this percent of data on either end AFTER the SSA has been calculated
plotGroupedComponents = False
plotENUSpectrogram = False
plotwCorMatrix = False
# --- --- --- FILTERING --- --- ---
SECTION_filterData = False if not calculateSSA else True
plotFilteredAxes = False
lowCut_toggle, highcut_toggle, filttype_toggle, order_toggle = 1, 63, 'Bandpass', 3 # filter toggles LOW FLYER
windowType, npersegN, scalingType  = 'hann', EFIsampleFreq, 'density' # spectrogram toggles
overlap = int(npersegN*(7/8)) # hanning filter overlap
# --- --- --- OUTPUT --- --- ---
outputData_dE = True if not calculateSSA else False
# --- --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import pandas as pd
from ACESII_code.myImports import *
from matplotlib.widgets import Slider
from ACESII_code.class_var_func import butter_filter
from numpy.fft import rfft, fftfreq
from scipy.interpolate import CubicSpline
from ACESII_code.supportCode.Support_Libraries.pymssa import MSSA
from scipy.signal import spectrogram

def EFI_L2_to_L2_dE(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L0_ACES_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*E_Field*')
    inputFiles_attitude = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[wflyer]}{modifier}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    fileoutName_dE = f'ACESII_{rocketID}_l2_dE'


    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Isolating Electric Field Instruments' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the Magnetometer file ---
        prgMsg(f'Loading data from {inputPath_modifier} Electric Field Files')
        data_dict_elec = loadDictFromFile(inputFiles[wFile], {})
        Done(start_time)

        # --- get the data from the Magnetometer file ---
        prgMsg(f'Loading data from {inputPath_modifier_attitude} Files')
        data_dict_attitude = loadDictFromFile(inputFiles_attitude[0], {})
        Done(start_time)

        ########################
        # --- Reduce dataset ---
        ########################
        prgMsg('Reducing Dataset')

        targetTimes = scienceRegions[wRocket - 4]

        # --- apply reduction to E_Field data ---
        lowCutoff, highCutoff = np.abs(data_dict_elec['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_elec['Epoch'][0] - targetTimes[1]).argmin()
        for key, val in data_dict_elec.items():
            data_dict_elec[key][0] = np.array(data_dict_elec[key][0][lowCutoff:highCutoff])

        Epoch_dt = np.array([tme for tme in deepcopy(data_dict_elec['Epoch'][0])])
        data_dict_elec['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_elec['Epoch'][0]])
        Epoch_seconds = np.array([(tme - data_dict_elec['Epoch'][0][0]) / 1E9 for tme in data_dict_elec['Epoch'][0]])


        # --- apply reduction to attitude data ---
        lowCutoff, highCutoff = np.abs(data_dict_attitude['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_attitude['Epoch'][0] - targetTimes[1]).argmin()
        for key, val in data_dict_attitude.items():
            data_dict_attitude[key][0] = np.array(data_dict_attitude[key][0][lowCutoff:highCutoff])

        data_dict_attitude['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude['Epoch'][0]])
        Done(start_time)

        # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
        # --- interpolate attitude data up to magnetometer epoch ---
        # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
        prgMsg('Interpolating Attitude Data')
        Epoch_attitude_loop = data_dict_attitude['Epoch'][0]
        dataKeys = ['Epoch', 'Alt', 'Latgd', 'Long', 'Y_Az', 'a11', 'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33']
        dataKeysVal = [deepcopy(data_dict_elec['Epoch'][0]), [], [], [], [], [], [], [], [], [], [], [], [], []]
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

        # --- form a single variable containing the EFI data ---
        data_for_output = np.array([ [data_dict_elec['E_East'][0][i],data_dict_elec['E_North'][0][i], data_dict_elec['E_Up'][0][i]] for i in range(len(data_dict_elec['Epoch'][0]))])

        if convertTomVpm:
            data_for_output = 1000*data_for_output

        # --- get the rocket's spin frequency for SSA work ---
        fitResults = {
            'Bx': {'Spin Amp': 25.42873940404161, 'Spin Freq': 0.6463295881639182, 'Spin Phase': 91.9759995936283, 'Cone Amp': 625.8772357084948, 'Cone Freq': 0.05294818121871208, 'Cone Phase': -138.77308595997619, 'Offset': -44919.748937299344},
            'By': {'Spin Amp': 7.378420193701481, 'Spin Freq': 0.6442248190622027, 'Spin Phase': 109.20255873087793, 'Cone Amp': 1380.5616077430786, 'Cone Freq': 0.02700105226961604, 'Cone Phase': 109.87799606103452, 'Offset': -139.74554466082876},
            'Bz': {'Spin Amp': 8.095746809541962, 'Spin Freq': 0.6442537451458561, 'Spin Phase': 19.11852573798773, 'Cone Amp': 1257.0313161879794, 'Cone Freq': 0.026874206798816504, 'Cone Phase': -69.78175516947503, 'Offset': 32.456720919269245}
        }
        spinFreq = sum([fitResults['Bz']['Spin Freq'], fitResults['By']['Spin Freq'], fitResults['Bz']['Spin Freq']]) / 3 if wRocket == 4 else 0.55
        Done(start_time)

        if SECTION_filterData:
            # --- --- --- --- -
            # --- FILTERING ---
            # --- --- --- --- -
            prgMsg('Filtering Data')

            # Apply Highpass filter to data
            E_ENU = data_for_output # get the data
            E_filtered = []
            for i in range(3):
                E_filtered.append(butter_filter(E_ENU[:, i], lowcutoff=lowCut_toggle, highcutoff=highcut_toggle, filtertype=filttype_toggle, order=order_toggle, fs=EFIsampleFreq))

            E_filtered = np.array(E_filtered)
            unit = 'mV/m' if convertTomVpm else 'V/m'
            comps = [f'E_East [{unit}]', f'E_North [{unit}]', f'E_Up [{unit}]']

            if plotFilteredAxes:

                for i in range(3):
                    ###############################
                    # --- Plot the Initial Data ---
                    ###############################
                    fig, ax = plt.subplots(nrows=4,ncols=1,constrained_layout=True)

                    ax[0].plot(Epoch_dt, E_ENU[:, i], label='noSpin')
                    ax[0].set_ylabel(f'{comps[i]}')

                    # --- FFT noSpin ---
                    N, T = len(E_ENU[:, i]), 1 / EFIsampleFreq
                    yf_rawData = rfft(E_ENU[:, i])

                    # --- Highpass Filter ---
                    filteredData = butter_filter(E_ENU[:, i], lowcutoff=lowCut_toggle, highcutoff=highcut_toggle, filtertype=filttype_toggle, order=order_toggle, fs=EFIsampleFreq)
                    filteredDataPlot, = ax[1].plot(Epoch_dt, filteredData, color='red')
                    ax[1].set_ylabel(f'{comps[i]}_filtered')
                    ax[1].set_xlabel('Time [s]')

                    # --- FFT filtered ---
                    N, T = len(E_ENU[:, i]), 1 / EFIsampleFreq
                    yf_filtered = rfft(filteredData)
                    xf = fftfreq(N, T)[:N // 2]
                    FFT_filtered_plot, = ax[2].plot(xf, 2.0 / N * np.abs(yf_filtered[0:N // 2]))
                    ax[2].plot(xf, 2.0 / N * np.abs(yf_rawData[0:N // 2]), color='orange')
                    ax[2].set_ylabel('FFT Power')
                    ax[2].set_xlabel('Frequency [Hz]')
                    ax[2].set_xlim(-0.1, 50)
                    ax[2].set_ylim(-0.01, 0.5)

                    # --- PERIODOGRAM ---
                    f, t, Sxx = spectrogram(filteredData, fs=EFIsampleFreq,
                                window=windowType,
                                nperseg=npersegN, # note: if ==None default size is 256
                                noverlap=overlap,
                                scaling=scalingType) # scaling = density or scaling = spectrum
                    spectrogramPlot = ax[3].pcolormesh(t, f, Sxx, shading='nearest',vmin=0,vmax=0.07,cmap='turbo')
                    cbar = plt.colorbar(spectrogramPlot,ax=ax[3])
                    ax[3].set_ylim(-0.1, 100)
                    ax[3].set_ylabel('Frequency [Hz]')
                    ax[3].set_xlabel('Epoch')

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
                    slider_cutoff_high = Slider(ax=axfilter_cutoff_high, label='highFq', valmin=0.001, valmax=highcut_toggle, valinit=1.5, orientation="vertical")

                    axfilter_order = fig.add_axes([0.945, 0.3, 0.02, 0.56])
                    slider_filter_order = Slider(ax=axfilter_order, label='order', valmin=1, valmax=15, valinit=1, orientation="vertical")

                    def f(cutoff_low, cutoff_high, order):
                        updated_filteredData = butter_filter(E_ENU[:, i], lowcutoff=cutoff_low, highcutoff=cutoff_high, order=int(order), filtertype=filttype_toggle, fs=EFIsampleFreq)
                        yf = rfft(updated_filteredData)
                        f, t, Sxx = spectrogram(updated_filteredData, fs=EFIsampleFreq,
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
            data_for_output = np.array([ [E_filtered[0][i],E_filtered[1][i],E_filtered[2][i]] for i in range(len(E_filtered[0]))])
            Done(start_time)

        if SECTION_SSA:
            # output file location for MSSA
            outputPathSSA = f'{rocketFolderPath}\\science\SSAcomponents_E\\{fliers[wflyer]}\\{fileoutName_dE}_SSAcomponents_WL{SSA_window_Size}'

            if useAlfvenRegion:
                outputPathSSA = outputPathSSA + '_Alfven.cdf'
            else:
                outputPathSSA = outputPathSSA + '.cdf'

            # name of the components
            compNames = ['E_East_SSA', 'E_North_SSA', 'E_Up_SSA']

            # create the MSSA object
            mssa = MSSA(n_components=None, window_size=SSA_window_Size, verbose=False)

            # --- perform the mSSA or SSA on specific components---
            if calculateSSA:

                prgMsg('Calculating SSA components')

                # convert data to pandas dataframe
                dataFormatted = {
                    'E_East': data_for_output[:, 0],
                    'E_North': data_for_output[:, 1],
                    'E_Up': data_for_output[:, 2]}

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
                SSAFiles = glob(f'{rocketFolderPath}\\science\SSAcomponents_E\\{fliers[wflyer]}\\*.cdf*')

                if justPrintSSAFiles:
                    ssa_names = [ssafile.replace(f'{rocketFolderPath}\\science\SSAcomponents_E\\{fliers[wflyer]}\\', '') for ssafile in SSAFiles]

                    for i, file in enumerate(ssa_names):
                        print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, ssa_names[i], round(getsize(SSAFiles[i]) / (10 ** 6), 1)))

                else:
                    # find the windowsize in the file I've chosen
                    windowSize = int(''.join(x for x in SSAFiles[wSSAFile].replace('36359','').replace('36364','').replace('l2','') if x.isdigit()))
                    data_dict_SSA = loadDictFromFile(SSAFiles[wSSAFile],{})

                    prgMsg('Grouping mSSA elements')
                    from ACESII_code.Processing.SSAgrouping_E import groupings
                    groupings = groupings(wRocket=wRocket, SSA_window_Size=windowSize, useAlfvenRegion=useAlfvenRegion)

                    # get all the SSA components for the three axes
                    E_SSA = [data_dict_SSA[compNames[0]][0], data_dict_SSA[compNames[1]][0], data_dict_SSA[compNames[2]][0]]

                    # --- Plot the FFT and groupings for one wAxes axes ---
                    if plotGroupedComponents:

                        # --- Plot it ---
                        fig, ax = plt.subplots(nrows=len(groupings),ncols=2)
                        fig.suptitle(compNames[wAxesSSA] + f'\n Window Length: {windowSize}')

                        # loop over all the groups in grouping
                        for i in range(len(groupings)):
                            data = np.array(E_SSA[wAxesSSA])

                            # combine the groupings
                            plotThisData = np.array([0 for i in range(len(data[:,0]))],dtype='float64')
                            for j, componentIndex in enumerate(groupings[i]):
                                plotThisData += np.array(data[:, componentIndex])

                            # reduce the last "X" percent of the data on either end to eliminate the SSA effect
                            percent = reduceTimePercent/100
                            percentLow, percentHigh = int((percent)*len(plotThisData)), int((1-percent)*len(plotThisData))
                            plotThisData = plotThisData[percentLow:percentHigh]

                            Epoch_dt = np.array([pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_elec['Epoch'][0]])
                            Epoch_dt = Epoch_dt[percentLow:percentHigh]

                            # plot the component
                            ax[i, 0].plot(Epoch_dt, plotThisData)

                            if i == 0:
                                ax[i, 0].set_ylabel('Original')
                            elif i == 1:
                                ax[i, 0].set_ylabel('Harmonics')
                            elif i == len(groupings) - 2:
                                ax[i, 0].set_ylabel('Noise')
                            elif i == len(groupings)-1:
                                ax[i, 0].set_ylabel('Physical Signal')
                                ax[i, 1].set_xlabel('Frequency [Hz]')
                                ax[i, 1].set_ylabel('FFT')
                            else:
                                ax[i, 0].set_ylabel('F{}'.format(i))

                            # calculate the FFT and plot it
                            N, T = len(plotThisData), 1 / EFIsampleFreq
                            yf, xf = rfft(plotThisData), fftfreq(N, T)[:N // 2]
                            FFT = 2.0 / N * np.abs(yf[0:N // 2])
                            ax[i, 1].plot(xf, FFT)
                            ax[i, 1].vlines([spinFreq*(i+1) for i in range(50)],ymin=min(FFT), ymax=max(FFT),alpha=0.5, color='red')
                            ax[i, 1].set_xlim(-0.1, 30)
                            ax[i, 1].set_ylim(0, max(FFT))

                        plt.show()

                    # --- collect the group info for the LAST (noise) grouping ---
                    if plotENUSpectrogram:

                        # --- group the three components---
                        groupedData = []

                        for i in range(3):# loop over all three axes

                            data = np.array(E_SSA[i]) # get the SSA data for this axes

                            # --- collect the group info for the LAST (noise) grouping ---
                            plotThisData = np.array([0 for i in range(len(data[:, 0]))], dtype='float64')
                            for componentIndex in groupings[-1]:
                                plotThisData += np.array(data[:, componentIndex])

                            # reduce the last "X" percent of the data on either end to eliminate the SSA effect
                            percent = reduceTimePercent / 100
                            percentLow, percentHigh = int((percent) * len(plotThisData)), int((1 - percent) * len(plotThisData))
                            Epoch_dt = [ pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_elec['Epoch'][0][percentLow:percentHigh]]

                            # append "noise" data to "groupedData"
                            groupedData.append(plotThisData[percentLow:percentHigh])

                        # --- Plot the Data ---
                        fig, ax = plt.subplots(nrows=2, ncols=3, sharex=True,layout='constrained')
                        plt.subplots_adjust(wspace=0.1, hspace=0)
                        fig.suptitle(f'ACES II {rocketID}\n'
                                     f'Window Length: {windowSize}')

                        for i in range(3):
                            # plot the SSA'd filtered data
                            ax[0, i].plot(Epoch_dt, groupedData[i])
                            ax[0, i].set_title(compNames[i])
                            ax[0, i].set_ylim(-12,12)

                            # plot the spectrogram of the data
                            f, t, Sxx = spectrogram(groupedData[i], fs=EFIsampleFreq,
                                                    window=windowType,
                                                    nperseg=npersegN,  # note: if ==None default size is 256
                                                    noverlap=overlap,
                                                    scaling=scalingType)  # scaling = density or scaling = spectrum

                            # convert spectrogram time to Epoch
                            spectroEpoch = np.array([pycdf.lib.tt2000_to_datetime(pycdf.lib.datetime_to_tt2000(Epoch_dt[0]) + int(1E9*tme)) for tme in t])
                            spectrogramPlot = ax[1, i].pcolormesh(spectroEpoch, f, Sxx, shading='nearest', vmin=-0.1, vmax=14, cmap='turbo')
                            if i == 2:
                                cbar = plt.colorbar(spectrogramPlot, ax=ax[:, i])
                                cbar.set_label('$B^{2}$/Hz')

                            if i == 0:
                                ax[0, i].set_ylabel('Magnetic Field [nT]')
                                ax[1, i].set_ylabel('Frequency [Hz]')

                            ax[1, i].set_ylim(-0.1, 15)
                            ax[1, i].set_xticklabels(ax[1, i].get_xticklabels(), rotation=45, ha='right')

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
                    newComps = ['dE_east', 'dE_north', 'dE_up']
                    for k in range(len(newComps)): # loop through all the components, but only take the last grouping

                        data = np.array(E_SSA[k])

                        # combine the groupings of the last group set only
                        formattedData = np.array([0 for i in range(len(data[:, 0]))], dtype='float64')
                        for j, componentIndex in enumerate(groupings[-1]):
                            formattedData += np.array(data[:, componentIndex])

                        # reduce the last "X" percent of the data on either end to eliminate the SSA effect
                        percent = reduceTimePercent / 100
                        percentLow, percentHigh = int((percent) * len(formattedData)), int((1 - percent) * len(formattedData))
                        data_for_output.append(formattedData[percentLow:percentHigh])

                    Epoch_SSA = data_dict_elec['Epoch'][0][percentLow:percentHigh]


                    # calculate dEmag
                    Emag = np.array([np.linalg.norm([data_for_output[0][i],data_for_output[1][i],data_for_output[2][i]]) for i in range(len(Epoch_SSA))])

                    # reformat the data for output
                    data_for_output = np.array([  [data_for_output[0][i],data_for_output[1][i],data_for_output[2][i], Emag[i]] for i in range(len(Epoch_SSA)) ])

                Done(start_time)

            # --- --- --- --- --- --- ---
            # --- WRITE OUT THE DATA ---
            # --- --- --- --- --- --- ---
            if outputData_dE:
                prgMsg('Creating dB output file')

                # create the output data_dict
                data_dict = deepcopy(data_dict_elec)
                comps = ['E_East', 'E_North', 'E_Up', 'E_mag']
                newCompsLabels = ['dE_East', 'dE_North', 'dE_Up', 'dE_mag']
                newLabels = ['&delta;E!BEast!n', '&delta;E!BNorth!n', '&delta;E!BUp!n', '&delta;E!Bmag!n']

                # --- Magnetic Components ---
                # get the attributes of the old components and replace them
                for i, key in enumerate(comps):
                    newAttrs = deepcopy(data_dict[key][1])
                    newAttrs['LABLAXIS'] = newLabels[i]

                    # remove the old key
                    del data_dict[key]

                    # append the new key
                    data_dict = {**data_dict, **{newCompsLabels[i]: [data_for_output[:, i], newAttrs]}}

                if subSECTION_groupSSAData:
                    data_dict['Epoch'][0] = Epoch_SSA

                outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName_dE}.cdf'

                outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, 'EFI')

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
        EFI_L2_to_L2_dE(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            EFI_L2_to_L2_dE(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            EFI_L2_to_L2_dE(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)