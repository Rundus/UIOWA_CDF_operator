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
wRocket = 4
modifier = '\\Stitching'
inputPath_modifier = '\L3\SSAcomponents_B' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier_db = 'L3\deltaB'

# --- --- --- WHICH DATA --- --- ---
wSSAFile_subset = 0
# --- --- --- SSA --- --- ---
wAxesSSA = 0 # 0 -> X, 1 -> Y, 2 -> Z
reduceTimePercent = 4 # kill this percent of data on either end AFTER the SSA has been calculated
plotGroupedComponents = True
# --- --- --- OUTPUT --- --- ---
Stitch_All_Files_Together_dB = False
# --- --- --- --- --- ---

spinFreq = 0.6442441031179716 if wRocket == 4 else 0.55


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.myImports import *
from numpy.fft import rfft, fftfreq
from ACESII_code.supportCode.Support_Libraries.pymssa import MSSA

def mSSA_to_deltaB(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L3'
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



        if Stitch_All_Files_Together_dB == False:

            componentFile = inputFiles[wSSAFile_subset].split('_SSAcomponents_')[1].replace('.cdf', '')
            temp1 = componentFile.split('_WL')
            temp2 = temp1[1].split('_')
            coordinates = temp1[0] # determine which coordinate system was used
            SSA_window_Size = int(temp2[0]) # determine which window length


            # --- get the data from the SSAcomponet file ---
            prgMsg(f'Loading data from {inputPath_modifier} RingCore Files')
            data_dict_SSAcomps = loadDictFromFile(inputFiles[wSSAFile_subset])
            Done(start_time)

            # prepare data for further processing
            comps = ['B_East', 'B_North', 'B_Up'] if coordinates in ['ENU','kenton'] else ['B_e', 'B_p', 'B_r']

            # name of the components in the SSA file
            compNames = [f'{comps[0]}_SSA', f'{comps[1]}_SSA', f'{comps[2]}_SSA']

            # create the MSSA object
            mssa = MSSA(n_components=None, window_size=SSA_window_Size, verbose=False)

            prgMsg('Grouping mSSA elements')
            from ACESII_code.Processing.Filtering.Stitching_mSSA.SSAgrouping_Stitch import groupings
            groupings = groupings(wRocket=wRocket, SSA_window_Size=SSA_window_Size, subsetNo=wSSAFile_subset)

            # get all the SSA components for the three axes
            B_SSA = [data_dict_SSAcomps[compNames[0]][0],
                     data_dict_SSAcomps[compNames[1]][0],
                     data_dict_SSAcomps[compNames[2]][0]]

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


        elif Stitch_All_Files_Together_dB == True:


            # --- format the data for output ---
            from ACESII_code.Processing.Filtering.Stitching_mSSA.SSAgrouping_Stitch import groupings
            Epoch_for_output = []
            data_for_output_X = []
            data_for_output_Y = []
            data_for_output_Z = []
            data_for_output_Bmag = []
            newComps = ['dB_East', 'dB_North', 'dB_Up'] if wSSAFile_subset in [0, 1] else ['dB_e', 'dB_p', 'dB_r']

            # Loop through all the SSAcomponent Files
            for j in range(len(inputFiles)):

                # load the data dict and get the epoch
                data_dict_SSAcomps = loadDictFromFile(inputFiles[j])

                # get all the info about the file
                componentFile = inputFiles[j].split('_SSAcomponents_')[1].replace('.cdf', '')
                temp1 = componentFile.split('_WL')
                temp2 = temp1[1].split('_')
                coordinates = temp1[0]  # determine which coordinate system was used
                SSA_window_Size = int(temp2[0])  # determine which window length
                comps = ['B_East', 'B_North', 'B_Up'] if coordinates in ['ENU', 'kenton'] else ['B_e', 'B_p', 'B_r']
                compNames = [f'{comps[0]}_SSA', f'{comps[1]}_SSA', f'{comps[2]}_SSA']

                # create the MSSA object
                mssa = MSSA(n_components=None, window_size=SSA_window_Size, verbose=False)

                # get the groupigs
                grouping = groupings(wRocket=wRocket, SSA_window_Size=SSA_window_Size, subsetNo=j)

                # combine the SSA components based on the groupings
                B_SSA = [data_dict_SSAcomps[compNames[0]][0],
                         data_dict_SSAcomps[compNames[1]][0],
                         data_dict_SSAcomps[compNames[2]][0]]

                # get all the data into
                for k in range(len(newComps)): # loop through all the components, but only take the last grouping
                    data = np.array(B_SSA[k])

                    # combine the groupings of the last group set only
                    formattedData = np.array([0 for i in range(len(data[:, 0]))], dtype='float64')
                    for j, componentIndex in enumerate(grouping[-1]):
                        formattedData += np.array(data[:, componentIndex])

                    if k == 0:
                        data_for_output_X.append(formattedData)
                    elif k == 1:
                        data_for_output_Y.append(formattedData)
                    elif k == 2:
                        data_for_output_Z.append(formattedData)

                # add in the Epoch
                Epoch_for_output.append(data_dict_SSAcomps['Epoch'][0])

                # determine Bmag
                data_for_output_Bmag.append([np.linalg.norm([data_for_output_X[j][i],
                                                             data_for_output_Y[j][i],
                                                             data_for_output_Z[j][i]]) for i in range(len(data_dict_SSAcomps['Epoch'][0]))])

            # collapse the dimensionality of the dataset and create a [4x1]xlen(data) object
            data_for_output = np.array([
                [item for sublist in data_for_output_X for item in sublist],
                [item for sublist in data_for_output_Y for item in sublist],
                [item for sublist in data_for_output_Z for item in sublist],
                [item for sublist in data_for_output_Bmag for item in sublist]]).T

            Epoch_SSA = np.array([item for sublist in Epoch_for_output for item in sublist])

            Done(start_time)

            # --- --- --- --- --- --- ---
            # --- WRITE OUT THE DATA ---
            # --- --- --- --- --- --- ---
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

            fileoutName_dB = f'ACESII_{rocketID}_l3_RingCore_dB_{coordinates}_WL{SSA_window_Size}_stitched_Flight'
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



if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    mSSA_to_deltaB(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)

