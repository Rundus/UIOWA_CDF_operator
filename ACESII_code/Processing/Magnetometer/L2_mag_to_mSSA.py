# --- L2_mag_to_mSSA.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Takes as input the despun B-Field data and filters, then mSSAs the data. The output are
# SSA component files found in \science\SSAcomponents_B. Grouping to get deltaB is done in mSSA_to_deltaB.py



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
wRocket = 4
modifier = ''
inputPath_modifier = 'l2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder


# --- --- --- WHICH DATA --- --- ---
# 0 --> ENU
# 1 --> model Subtracted
# 2 --> Field Aligned
wUseData = 0
# --- --- --- REDUCE DATA --- --- ---
reduceDataSet = True
# 0 --> Alfven Region
# 1 --> Kenton's Zoomed Region
# 2 --> The Whole flight, excluding outlier regions
wTargetTimes = 2
# --- --- --- FILTERING --- --- ---

# NOTE!!! I have outsourced the filtering of the Despun Data using MATLAB!. Keep SECTION_filterData = False
SECTION_filterData = True
plotFilteredAxes = False
lowCut_toggle, highcut_toggle, filttype_toggle, order_toggle = 0.5, 20, 'Bandpass', 4 # filter toggles LOW FLYER
# --- --- --- SSA --- --- ---
SECTION_SSA = True
SSA_window_Size = 1201




# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import pandas as pd
from matplotlib.widgets import Slider
from ACESII_code.class_var_func import butter_filter
from numpy.fft import rfft, fftfreq
from ACESII_code.supportCode.Support_Libraries.pymssa import MSSA
from scipy.signal import spectrogram



def RingCore_L1_to_L2_Despin(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L0_ACES_Quick(wflyer)

    if wUseData == 0: # ENU
        searchMod = 'ENU'
    elif wUseData == 1: # Model Subtracted
        searchMod = 'subModel'
    elif wUseData == 2: # Field Aligned
        searchMod = 'FieldAligned'

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*RingCore_DeSpun_{searchMod}.cdf*')
    # inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*Geo*')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]

    input_names_searchable = [ifile.replace(inputPath_modifier.lower() +'_', '') for ifile in input_names]


    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'DeSpining RingCore Data' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the Magnetometer file ---
        prgMsg(f'Loading data from {inputPath_modifier} RingCore Files')
        from ACESII_code.Processing.Magnetometer.SSAgrouping_and_target_times_B import timeWindow
        targetTimes = timeWindow(wTargetTimes, wRocket)
        data_dict_mag = loadDictFromFile(inputFiles[wFile],{},reduceData=reduceDataSet,targetTimes=targetTimes)
        Epoch_seconds = np.array([(tme - data_dict_mag['Epoch'][0][0]) / 1E9 for tme in data_dict_mag['Epoch'][0]])
        Epoch_dt = np.array([tme for tme in data_dict_mag['Epoch'][0]])
        Done(start_time)

        # prepare data for further processing
        comps = ['B_East', 'B_North', 'B_Up'] if wUseData in [0, 1] else ['B_e', 'B_p', 'B_r']
        # comps = ['Bx', 'By', 'Bz']
        data_for_output = np.array([[data_dict_mag[comps[0]][0][i],
                                     data_dict_mag[comps[1]][0][i],
                                     data_dict_mag[comps[2]][0][i]] for i in range(len(Epoch_seconds))])


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

            if plotFilteredAxes:

                windowType, npersegN, scalingType = 'hann', 128, 'density'  # spectrogram toggles
                overlap = int(npersegN * (7 / 8))  # hanning filter overlap

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
                    # ax[2].set_ylim(-0.1, )

                    # --- PERIODOGRAM ---
                    f, t, Sxx = spectrogram(filteredData, fs=128,
                                window=windowType,
                                nperseg=npersegN, # note: if ==None default size is 256
                                noverlap=overlap,
                                scaling=scalingType) # scaling = density or scaling = spectrum
                    spectrogramPlot = ax[3].pcolormesh(t, f, Sxx, shading='nearest',vmin=0,vmax=10,cmap='turbo')
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
                    slider_cutoff_low = Slider(ax=axfilter_cutoff_low, label='lowFq', valmin=0.0001, valmax=4, valinit=1, orientation="vertical")

                    axfilter_cutoff_high = fig.add_axes([0.97, 0.3, 0.02, 0.63])
                    slider_cutoff_high = Slider(ax=axfilter_cutoff_high, label='highFq', valmin=0.00001, valmax=63, valinit=2, orientation="vertical")

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
                        FFT = 2.0 / N * np.abs(yf_new[0:N // 2])
                        FFT_filtered_plot.set_ydata(FFT)

                        # adjust the y-scale of filtered data
                        ax[1].set_ylim(min(newData), max(newData))
                        ax[2].set_ylim(min(FFT), max(FFT))

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
            outputPathSSA = f'{rocketFolderPath}\\science\SSAcomponents_B\\{fliers[wflyer]}\\ACESII_{rocketID}_RingCore_SSAcomponents_{searchMod}_WL{SSA_window_Size}'

            if wTargetTimes == 0:
                outputPathSSA = outputPathSSA + '_Alfven.cdf'
            elif wTargetTimes == 1:
                outputPathSSA = outputPathSSA + '_kenton.cdf'
            else:
                outputPathSSA = outputPathSSA + '_flight.cdf'


            # name of the components
            compNames = [f'{comps[0]}_SSA', f'{comps[1]}_SSA', f'{comps[2]}_SSA']

            # create the MSSA object
            mssa = MSSA(n_components=None, window_size=SSA_window_Size, verbose=False)

            # --- perform the mSSA or SSA on specific components---
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
                dataToOutput = np.array(components[i, :, :]) # update the data for output to be the components
                attrs = deepcopy(example_attrs)
                attrs['LABLAXIS'] = compNames[i]
                attrs['VALIDMIN'] = dataToOutput.min()
                attrs['VALIDMAX'] = dataToOutput.max()
                data_dict_SSAcomps = {**data_dict_SSAcomps, **{compNames[i]:[dataToOutput, attrs]}}

            # add in the Epoch variable
            data_dict_SSAcomps = {**data_dict_SSAcomps, **{'Epoch':[data_dict_mag['Epoch'][0],data_dict_mag['Epoch'][1]]}}

            outputCDFdata(outputPathSSA, data_dict_SSAcomps, outputModelData, globalAttrsMod, 'RingCore')
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
    RingCore_L1_to_L2_Despin(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)

