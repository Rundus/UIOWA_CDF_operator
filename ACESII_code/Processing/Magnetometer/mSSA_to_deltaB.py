# # --- mSSA_to_deltaB.py ---
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
modifier = ''
inputPath_modifier = '\science\SSAcomponents_B' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier_db = 'science\deltaB'


# --- --- --- WHICH DATA --- --- ---
wSSAFile = 0
# --- --- --- SSA --- --- ---
wAxesSSA = 0 # 0 -> X, 1 -> Y, 2 -> Z
reduceTimePercent = 4 # kill this percent of data on either end AFTER the SSA has been calculated
plotGroupedComponents = False
plotSpectrogram = False
plotwCorMatrix = False
# --- --- --- OUTPUT --- --- ---
outputData_dB = True
# --- --- --- --- --- ---

spinFreq = 0.6442441031179716 if wRocket == 4 else 0.55


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import pandas as pd
from ACESII_code.myImports import *
from numpy.fft import rfft, fftfreq
from ACESII_code.supportCode.Support_Libraries.pymssa import MSSA
from scipy.signal import spectrogram



def mSSA_to_deltaB(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L0_ACES_Quick(wflyer)

    # get the specific input file
    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf*')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() + '_','').replace('_v00', '') for ifile in input_names]


    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'DeSpining RingCore Data' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')


        componentFile = inputFiles[wSSAFile].split('_SSAcomponents_')[1].replace('.cdf', '')
        temp1 = componentFile.split('_WL')
        temp2 = temp1[1].split('_')
        coordinates = temp1[0] # determine which coordinate system was used
        SSA_window_Size = int(temp2[0]) # determine which window length
        TimeWindow = temp2[1] # determine which time region was used

        if TimeWindow in ['Alfven']:
            wUseData = 0
        elif TimeWindow in ['kenton']:
            wUseData = 1
        elif TimeWindow in ['flight']:
            wUseData = 2

        # --- output file name ---
        fileoutName_dB = f'ACESII_{rocketID}_l2_RingCore_dB_{coordinates}_WL{SSA_window_Size}_{TimeWindow}'

        # --- get the data from the SSAcomponet file ---
        prgMsg(f'Loading data from {inputPath_modifier} RingCore Files')
        data_dict_SSAcomps = loadDictFromFile(inputFiles[wSSAFile],{},reduceData=False,targetTimes=[])
        Done(start_time)

        # prepare data for further processing
        comps = ['B_East', 'B_North', 'B_Up'] if coordinates in ['ENU','kenton'] else ['B_e', 'B_p', 'B_r']

        # name of the components in the SSA file
        compNames = [f'{comps[0]}_SSA', f'{comps[1]}_SSA', f'{comps[2]}_SSA']

        # create the MSSA object
        mssa = MSSA(n_components=None, window_size=SSA_window_Size, verbose=False)

        prgMsg('Grouping mSSA elements')
        from ACESII_code.Processing.Magnetometer.SSAgrouping_and_target_times_B import groupings
        groupings = groupings(wRocket=wRocket, SSA_window_Size=SSA_window_Size, wUseData=wUseData)

        # get all the SSA components for the three axes
        B_SSA = [data_dict_SSAcomps[compNames[0]][0], data_dict_SSAcomps[compNames[1]][0], data_dict_SSAcomps[compNames[2]][0]]

        # --- Plot the FFT and groupings for one wAxes axes ---
        if plotGroupedComponents:

            # --- Plot it ---
            fig, ax = plt.subplots(nrows=len(groupings),ncols=2)
            fig.suptitle(compNames[wAxesSSA] + f'\n Window Length: {SSA_window_Size}')

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

                Epoch_dt = data_dict_SSAcomps['Epoch'][0][percentLow:percentHigh]

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
                    ax[i, 1].set_ylabel(groupings[i])

                # calculate the FFT and plot it
                N, T = len(plotThisData), 1 / 128
                yf, xf = rfft(plotThisData), fftfreq(N, T)[:N // 2]
                FFT = 2.0 / N * np.abs(yf[0:N // 2])
                ax[i, 1].plot(xf, FFT)
                ax[i, 1].vlines([spinFreq*(i+1) for i in range(50)],ymin=min(FFT), ymax=max(FFT),alpha=0.5, color='red')
                ax[i, 1].set_xlim(-0.1, 10)
                ax[i, 1].set_ylim(0, max(FFT))

            plt.show()

        # --- plot spectrogram of the Data after the grouping---
        if plotSpectrogram:

            windowType, npersegN, scalingType = 'hann', 128, 'density'  # spectrogram toggles
            overlap = int(npersegN * (7 / 8))  # hanning filter overlap

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
                Epoch_dt = [ pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_SSAcomps['Epoch'][0][percentLow:percentHigh]]


                # append "noise" data to "groupedData"
                groupedData.append(plotThisData[percentLow:percentHigh])


            # --- Plot the Data ---
            fig, ax = plt.subplots(nrows=2, ncols=3, sharex=True,layout='constrained')
            plt.subplots_adjust(wspace=0.1, hspace=0)
            fig.suptitle(f'ACES II {rocketID}\n'
                         f'Window Length: {SSA_window_Size}')

            for i in range(3):
                # plot the SSA'd filtered data
                ax[0, i].plot(Epoch_dt, groupedData[i])
                ax[0, i].set_title(compNames[i])
                ax[0, i].set_ylim(-12,12)

                # plot the spectrogram of the data
                f, t, Sxx = spectrogram(groupedData[i], fs=128,
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

        # --- collect the group info for the LAST (noise) grouping ---
        if plotwCorMatrix:

            # plot all three axes correlation matrix
            mssa.fit(pd.DataFrame(
                {'Data1': [i for i in range(len(data_dict_SSAcomps['Epoch'][0]))],
                 'Data2': [i for i in range(len(data_dict_SSAcomps['Epoch'][0]))],
                 'Data3': [i for i in range(len(data_dict_SSAcomps['Epoch'][0]))]}
            ))

            for i in range(3):
                mssa.components_[i, :, :] = data_dict_SSAcomps[compNames[i]][0]

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
        newComps = ['dB_East', 'dB_North', 'dB_Up'] if wUseData in [0,1] else ['dB_e', 'dB_p', 'dB_r']
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

        Epoch_SSA = data_dict_SSAcomps['Epoch'][0][percentLow:percentHigh]

        # calculate dBmag
        data_for_output.append(np.array([np.linalg.norm([data_for_output[0][i],data_for_output[1][i],data_for_output[2][i]]) for i in range(len(data_for_output[0]))]))

        # reformat the data for output
        data_for_output = np.array([[data_for_output[0][i],data_for_output[1][i],data_for_output[2][i],data_for_output[3][i]] for i in range(len(data_for_output[0])) ])

        Done(start_time)

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData_dB:
            prgMsg('Creating dB output file')

            # create the output data_dict
            data_dict = deepcopy(data_dict_SSAcomps)
            newComps = ['dB_East', 'dB_North', 'dB_Up', 'dB_Mag']
            newLabels = ['&delta;B!BEast!n', '&delta;B!BNorth!n', '&delta;B!BUp!n', '&delta;B!BMag!n']

            # --- Magnetic Components ---
            # get the attributes of the old components and replace them
            for i, key in enumerate(compNames):
                newAttrs = deepcopy(data_dict[key][1])
                newAttrs['LABLAXIS'] = newLabels[i]
                newAttrs['DEPEND_0'] = 'Epoch'

                # remove the old key
                del data_dict[key]

                # append the new key
                data_dict = {**data_dict, **{newComps[i] : [data_for_output[:, i], newAttrs]}}

            # --- update the reduced epoch ---
            data_dict['Epoch'][0] = Epoch_SSA

            outputPath = f'{rocketFolderPath}{outputPath_modifier_db}\{fliers[wflyer]}\\{fileoutName_dB}.cdf'

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
    mSSA_to_deltaB(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)

