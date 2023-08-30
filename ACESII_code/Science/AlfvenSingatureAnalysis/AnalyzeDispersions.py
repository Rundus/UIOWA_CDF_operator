# --- AnalyzeDispersions.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Pull the L1 eepaa data and look at the specific locations of the alfven signature
# and perform an analysis on this dataset to study them further



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False

wRocket = 4

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_trajectory = 'trajectories'
outputPath_modifier = 'science\AlfvenSignatureAnalysis' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

######################################
# --- PLOT THE DISPERSION FEATURES ---
######################################
# plot all of the dispersion functions over a range of pitch angles (user input)
plotKeyDispersions = True
wDispersions = [2] # [] -> plot all dispersion traces, [#,#,#,...] plot specific ones. USE THE DISPERSION NUMBER NOT PYTHON -1 INDEX
wPitches = [2] # plots specific pitch angles by their index
isolateAlfvenSignature = True # removes unwanted data from the alfven signature
fitCurveToData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.optimize import curve_fit
from ACESII_code.class_var_func import Re


def AlfvenSignatureAnalysis(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer,wDispersions):


    wDispersions = [dispersionNo-1 for dispersionNo in wDispersions]

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')
    inputFiles_trajectory = glob(f'{rocketFolderPath}{inputPath_modifier_trajectory}\{fliers[wflyer]}{modifier}\*_ILat_ILong*')
    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    input_names_trajectory = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles_trajectory]
    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')
    dataFile_name_trajectory = inputFiles_trajectory[0].replace(f'{rocketFolderPath}{inputPath_modifier_trajectory}\{fliers[wflyer]}{modifier}\\', '')

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            descriptiorNam = ['EEPAA', 'LEESA', 'IEPAA', 'Langmuir_Probe']
            wInstr = [index, instr, descriptiorNam[index]]

    fileoutName = f'ACESII_{rocketID}_{wInstr[1]}_InvertedV_removed.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the ESA file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict = loadDictFromFile(inputFiles[wFile],{})
        data_dict['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_esa'][0][i]) for i in (range(len(data_dict['Epoch_esa'][0])))])
        Done(start_time)

        # --- get data from Trajectory Files ---
        prgMsg(f'Loading data from {inputPath_modifier_trajectory} Files')
        data_dict_traj = loadDictFromFile(inputFiles_trajectory[0], {})
        data_dict_traj['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch_esa'][0][i]) for i in (range(len(data_dict_traj['Epoch_esa'][0])))])
        Done(start_time)

        # Plot the key dispersion colormaps
        if plotKeyDispersions:
            from ACESII_code.Science.AlfvenSingatureAnalysis.dispersionAttributes import dispersionAttributes

            if len(wDispersions) == 0:
                nColumns = 3
                nRows = len(wPitches)
            else:
                nColumns = len(wDispersions)
                nRows = len(wPitches)

            fig = plt.figure(1, figsize=(20,20))
            axes = []

            # plot each dispersion
            counter = 1

            for rowI in range(nRows):
                for colI in range(nColumns):

                    wPitch = wPitches[rowI] # the pitch index you're plotting
                    wDispersion = wDispersions[colI] # the dispersion index you're plotting

                    # isolate a single dispersion trace
                    targetTimes = [pycdf.lib.datetime_to_tt2000(dispersionAttributes.keyDispersionTimes[wDispersion][0]),
                                   pycdf.lib.datetime_to_tt2000(dispersionAttributes.keyDispersionTimes[wDispersion][1])]
                    limitIndexes = [np.abs(data_dict['Epoch_esa'][0] - targetTimes[0]).argmin(),
                                    np.abs(data_dict['Epoch_esa'][0] - targetTimes[1]).argmin()]
                    EpochOneDis = np.array([data_dict['Epoch_esa'][0][i] for i in range(limitIndexes[0], limitIndexes[1])])
                    esaDataOneDis = np.array(data_dict[wInstr[1]][0][limitIndexes[0]:limitIndexes[1]])

                    # select only the energies from dispersionAttributes.energyLimits
                    esaDataOneDis = np.array([esaDataOneDis[tme][wPitch][dispersionAttributes.keyDispersionEnergyLimits[wDispersion]:] for tme in range(len(esaDataOneDis))],dtype='int')
                    Energy = data_dict['Energy'][0][dispersionAttributes.keyDispersionEnergyLimits[wDispersion]:]
                    Pitch = data_dict['Pitch_Angle'][0]

                    # convert time subset to start t=0 at the beginning of the dispersion
                    timeLaunch = EpochOneDis.min()
                    EpochOneDis = np.array((EpochOneDis - timeLaunch) / 1E9)


                    # mask the data by maskval
                    maskval = 1
                    esaDataOneDis = np.array([[(esaDataOneDis[k][j] - maskval) if (esaDataOneDis[k][j] - maskval) > 0 and esaDataOneDis[k][j] >= 0 else 0 for j in range(len(esaDataOneDis[k]))]for k in range(len(esaDataOneDis))])

                    # apply the isolation functions found in dispersionAttributes.py
                    if isolateAlfvenSignature:
                        esaDataOneDis = dispersionAttributes.isolationFunctions[f's{wDispersion+1}'](esaDataOneDis, Energy, EpochOneDis)

                    ########################
                    # create the x/y dataset
                    ########################
                    xData, yData, zData, xData_for_fitting, yData_for_fitting = [], [], [], [], []  # time, energy, counts
                    for tme in range(len(esaDataOneDis)):
                        for engy in range(len(esaDataOneDis[0])):
                            for j in range(int(esaDataOneDis[tme][engy])): # collect the data for the fitting
                                xData_for_fitting.append(EpochOneDis[tme])
                                yData_for_fitting.append(Energy[engy])

                            if int(esaDataOneDis[tme][engy]) != 0: # collect the data for the colorplot
                                xData.append(EpochOneDis[tme])
                                yData.append(Energy[engy])
                                zData.append(int(esaDataOneDis[tme][engy]))

                    ###############
                    # sort the data
                    ###############
                    xData, yData, zData = map(np.array, zip(*sorted(zip(xData, yData, zData))))
                    xData_for_fitting, yData_for_fitting = map(np.array, zip(*sorted(zip(xData_for_fitting,yData_for_fitting))))

                    ########################
                    # --- CREATE SUBPLOT ---
                    ########################
                    ax = fig.add_subplot(nRows, nColumns, counter)
                    axes.append(ax)

                    xx, yy = np.meshgrid(EpochOneDis,Energy)
                    cmap = ax.pcolormesh(xx, yy, esaDataOneDis.T, vmin=1, vmax=40, cmap='turbo')
                    # ax.set_yscale('log')
                    # cmap = ax.scatter(xData, yData, c=zData, vmin=0, vmax=20, cmap='plasma')

                    if rowI == 0:
                        ax.set_title(f's{wDispersion + 1}', fontsize=15)
                    if colI == 0:
                        ax.set_ylabel(fr'$\alpha$ = {Pitch[wPitch]}$^\circ$'+' \n Energy [eV]',fontsize=11)
                    if rowI == nRows-1:
                        ax.set_xlabel('Time [seconds]', fontsize=11)

                    if fitCurveToData:
                        # --- --- --- --- --- -
                        # --- CURVE FITTING ---
                        # --- --- --- --- --- -
                        def fitFunc(x, A, B, C):
                            y = np.cos(np.radians(Pitch[wPitch]))*(0.5*(m_e/q0)*(A/(x-B))**2)
                            return y

                        # params, cov = curve_fit(fitFunc, xData_for_fitting, yData_for_fitting, p0=[400,-1,-80], maxfev=10000)
                        params, cov = curve_fit(fitFunc, xData_for_fitting, yData_for_fitting, maxfev=10000)
                        xDataFit = np.linspace(xData_for_fitting.min(), xData_for_fitting.max(), 200)
                        yDataFit = np.array([fitFunc(x, *params) for x in xDataFit])
                        ax.plot(xDataFit, yDataFit, color='red', linewidth=5,label=f'$\Delta z$: {round((params[0]/(1000*Re)),2)  } $R_E$')
                        ax.legend(loc='upper right')

                    # keep this line at the end of the row/col loop!
                    counter += 1


            fig.suptitle(f'ACESII Alfven Dispersion',fontsize=25)
            plt.subplots_adjust(wspace=0)
            fig.colorbar(cmap, ax=axes,orientation="horizontal", pad=0.1,fraction=0.025,aspect=100,label='counts',norm='log')
            plt.show()




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
        AlfvenSignatureAnalysis(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer,wDispersions)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            AlfvenSignatureAnalysis(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer,wDispersions)
    else:
        for filesNo in wFiles:
            AlfvenSignatureAnalysis(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer,wDispersions)