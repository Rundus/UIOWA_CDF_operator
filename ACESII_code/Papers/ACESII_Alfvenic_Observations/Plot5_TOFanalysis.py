# --- AnalyzeDispersions.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Isolate Alfvenic Dispersions and perform cross-correlation analysis
# to determine values that are fitted linearly to get the height of the source region

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

import ACESII_code.class_var_func
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---
justPrintFileNames = False

wRocket = 4
modifier = ''
inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\AlfvenSignatureAnalysis' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
# plot all of the dispersion functions over a range of pitch angles (user input)
wDispersions = [10] # [] -> plot all dispersion traces, [#,#,#,...] plot specific ones. USE THE DISPERSION NUMBER NOT PYTHON -1 INDEX
wPitch = 2 # plots specific pitch angles by their index
# ---------------------------
justPlotKeyDispersions = False #IF ==TRUE no cross-correlation will occur
# ---------------------------
applyMaskVal = True
maskVal = 2 # apply a single mask to the dispersion feature
# ---------------------------
isolateAlfvenSignature = False # removes unwanted data from the alfven signature
# ---------------------------
plotCorrelationProcess = False
DetectorTimeResolution = 0.05 # in seconds
DetectorEnergyResolution = 0.18
# ---------------------------
correlationAnalysis = True
showErrorBars = False
weightLinearFitByCounts = True
outputCorrelationPlot = True
# ---------------------------
plotZsourceVsTime = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.signal import correlate, correlation_lags
from itertools import combinations
from ACESII_code.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes
from ACESII_code.class_var_func import Re

def AlfvenSignatureCrossCorrelation(wRocket, rocketFolderPath, justPrintFileNames,wDis):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    globalAttrsMod = rocketAttrs.globalAttributes[wRocket-4]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L1_ACES_Quick(wRocket-4)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}{modifier}\*eepaa_fullcal*')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}{modifier}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    # --- get the data from the ESA file ---
    prgMsg(f'Loading data from {inputPath_modifier} Files')
    data_dict = loadDictFromFile(inputFiles[0])
    Done(start_time)

    # --- --- --- --- --- --- ----
    # --- Isolate a Dispersion ---
    # --- --- --- --- --- --- ----
    prgMsg('Isolating Dispersion')

    wDispersion_key = f's{wDis}'
    Energy = data_dict['Energy'][0]
    Pitch = data_dict['Pitch_Angle'][0]
    lowCut, highCut = np.abs(data_dict['Epoch'][0] - dispersionAttributes.keyDispersionTimes[wDis-1][0]).argmin(), np.abs(data_dict['Epoch'][0] - dispersionAttributes.keyDispersionTimes[wDis-1][1]).argmin()
    Epoch_dis = dateTimetoTT2000(data_dict['Epoch'][0][lowCut:highCut+1],inverse=False)
    Epoch_dis = (np.array(Epoch_dis) - Epoch_dis[0]) / 1E9 # converted data to TIME SINCE START OF DISPERSION (needed for proper plot fitting)

    # calculate the center point (in time) of the dispersion
    whenSTEBoccured_time = data_dict['Epoch'][0][int((highCut+lowCut)/2)]
    whenSTEBoccured_Alt = data_dict['Alt'][0][int((highCut + lowCut) / 2)]

    eepaa_dis = data_dict['eepaa'][0][lowCut:highCut+1]

    if isolateAlfvenSignature:
        eepaa_dis = np.array(dispersionAttributes.isolationFunctions[wDispersion_key](eepaa_dis, Energy, Epoch_dis)) # pply the isolation functions found in dispersionAttributes.py

    if applyMaskVal:
        # --- Apply the Masking Value to the dataset ---
        # Note: This removes all Fillvals
        eepaa_dis[eepaa_dis<=maskVal] = 0 # applies mask and zeros-out anything below 0
        eepaa_dis = eepaa_dis.clip(min=0)
        eepaa_dis[eepaa_dis>1E15] = 0 # zeros-out anything above 1E15

    # --- Reduce data to Pitch Slice of Interest ---
    # gets one slice in pitch of the data and flips
    # data from dimensions 0 == Time to --> dimension 0 == Energy so that the first index gets energy for all time
    eepaa_dis_onePitch = np.array(eepaa_dis[:,wPitch,:]).T
    Done(start_time)

    #################################################
    # --- PLOT THE DISPERSION IN ORDER TO ISOLATE ---
    #################################################

    if justPlotKeyDispersions:

        fig, ax = plt.subplots(2)
        figure_width = 25
        figure_height = 15
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)
        fig.suptitle(f'STEB {wDis}\n Pitch = {Pitch[wPitch]}deg\n {data_dict["Epoch"][0][lowCut].strftime("%H:%M:%S.%f")} to {data_dict["Epoch"][0][highCut].strftime("%H:%M:%S.%f")}')
        cmap = ax[0].pcolormesh(Epoch_dis, Energy, np.array(eepaa_dis[:, wPitch, :]).T, cmap='turbo', vmin=0, vmax=np.array(eepaa_dis[:, wPitch, :]).max())
        ax[0].set_yscale('log')
        ax[0].set_ylim(Energy[-1], Energy[np.abs(Energy - dispersionAttributes.keyDispersionEnergyLimits[wDis-1][1]).argmin()-1])
        ax[1].pcolormesh(Epoch_dis, Energy, np.array(data_dict['eepaa'][0][lowCut:highCut+1,wPitch,:]).T, cmap='turbo', vmin=0, vmax=np.array(eepaa_dis[:, wPitch, :]).max())
        ax[1].set_yscale('log')
        ax[1].set_ylim(Energy[-1], Energy[np.abs(Energy - dispersionAttributes.keyDispersionEnergyLimits[wDis - 1][1]).argmin() - 1])
        cbar = plt.colorbar(cmap,ax=ax.ravel().tolist())
        plt.show()
    elif correlationAnalysis:

        ############################################
        # --- PERFORM CROSS CORRELATION ANALYSIS ---
        ############################################

        prgMsg('Performing Cross-correlation Analysis')

        # determine the energy pairs to be used in the correlation analysis via combination statistics
        engyIndexLow, engyIndexHigh = np.abs(Energy - dispersionAttributes.keyDispersionEnergyLimits[wDis-1][0]).argmin(),np.abs(Energy - dispersionAttributes.keyDispersionEnergyLimits[wDis-1][1]).argmin()
        engysToAnalyze = Energy[engyIndexHigh:engyIndexLow+1]
        engysIndicesToAnalyze = np.arange(engyIndexHigh,engyIndexLow+1,1)
        energyPairs = [comb for comb in combinations(engysIndicesToAnalyze,2)]


        # for each pair of Energies, perform the cross-correlation analysis:
        # Note: the peak value in the lag-time determines the deltaT value between the velocities
        deltaTs = []
        deltaVs = []

        errorT = []
        errorV = []
        errorZ = []

        for engyPair in energyPairs:

            higherEnergyIndex = engyPair[0]
            lowerEnergyIndex = engyPair[1]
            EnergySlow = q0*Energy[lowerEnergyIndex]
            EnergyFast = q0*Energy[higherEnergyIndex]
            velSlow = np.cos(np.radians(Pitch[wPitch]))*np.sqrt(2*EnergySlow/m_e)
            errorVelSlow = np.sqrt(DetectorEnergyResolution) * velSlow
            velFast = np.cos(np.radians(Pitch[wPitch]))*np.sqrt(2*EnergyFast/m_e)
            errorVelFast = np.sqrt(DetectorEnergyResolution) * velFast

            # check to make sure both datasets have at least one non-zero value
            if not np.any(eepaa_dis_onePitch[higherEnergyIndex]) or not np.any(eepaa_dis_onePitch[lowerEnergyIndex]):
                continue
            else:
                corr = np.array(correlate(eepaa_dis_onePitch[higherEnergyIndex], eepaa_dis_onePitch[lowerEnergyIndex]))
                corr_norm = corr/np.max(corr)
                lags = correlation_lags(len(eepaa_dis_onePitch[higherEnergyIndex]), len(eepaa_dis_onePitch[lowerEnergyIndex]))

            if plotCorrelationProcess:
                fig, ax = plt.subplots(4)
                figure_width = 10
                figure_height = 20
                fig.set_figwidth(figure_width)
                fig.set_figheight(figure_height)
                fig.suptitle(f'High Energy {Energy[higherEnergyIndex]} \n Low Energy {Energy[lowerEnergyIndex]}')
                ax[0].pcolormesh(Epoch_dis,Energy,np.array(eepaa_dis[:,wPitch,:]).T,cmap='turbo',vmin=0,vmax=40)
                ax[0].set_yscale('log')
                ax[0].set_ylim(Energy[-1],Energy[higherEnergyIndex-1])
                ax[0].axhline(Energy[lowerEnergyIndex], color='green')
                ax[0].axhline(Energy[higherEnergyIndex], color='red')
                ax[1].plot(eepaa_dis_onePitch[higherEnergyIndex])
                ax[1].set_ylabel('Higher Energy Counts')
                ax[2].plot(eepaa_dis_onePitch[lowerEnergyIndex])
                ax[2].set_ylabel('Lower Energy Counts')
                ax[3].plot(lags, corr)
                ax[3].set_ylabel('Correlation (normalized)')
                ax[3].set_xlabel('Lag')
                plt.tight_layout()
                plt.show()

            try:
                # Find the x,y value of the peak in the correlation output
                indexMax = np.where(corr == np.max(corr))[0][0]
                delayTime = -1*DetectorTimeResolution*lags[indexMax]
                velDiff = 1000*(np.sqrt(m_e) / (np.cos(np.radians(Pitch[wPitch]))*(np.sqrt(2)))) * (1/(np.sqrt(EnergySlow)) - 1/(np.sqrt(EnergyFast)))

                if weightLinearFitByCounts:
                    # find if counts_slow or counts_fast is more at the correlation peak
                    N_fast = np.array(eepaa_dis_onePitch[higherEnergyIndex]).max()
                    N_Slow = np.array(eepaa_dis_onePitch[lowerEnergyIndex]).max()

                    # whichever is less (N_Fast,N_slow) append that many values to the results in order to bias the linear fit
                    # towards values that are both higher in counts, instead of just one
                    for h in range(min(N_fast, N_Slow)):
                        deltaTs.append(delayTime)
                        deltaVs.append(velDiff)
                else:
                    deltaTs.append(delayTime)
                    deltaVs.append(velDiff)

                # calculate the error in the measurements
                errorT.append(DetectorTimeResolution)
                errorV.append(((DetectorEnergyResolution*np.sqrt(m_e)) / ((2**(3/2))*np.cos(np.radians(Pitch[wPitch])))) * np.sqrt(1/EnergySlow + 1/EnergyFast))

                errZ1 = (DetectorTimeResolution**2) * (1 / (1/velSlow - 1/velFast))
                errZ2 = (np.square(errorVelSlow))*(delayTime*np.square(velFast)/np.square(velSlow - velFast))**2
                errZ3 = (np.square(errorVelFast))*(delayTime*np.square(velSlow)/np.square(velSlow - velFast))**2
                errorZ.append(errZ1 + errZ2 + errZ3)

            except:
                print('indexMax = np.where(corr == np.max(corr))[0][0] Found no maximum. Correlation data is likely all zeros. Did you set your maskval too high?')

        deltaTs = np.array(deltaTs)
        deltaVs = np.array(deltaVs)

        # convert error in velocity to kilometers
        errorV = np.array(errorV)*1000

        Done(start_time)

        ##########################################
        # --- FIT THE DATA TO DERIVE TOF Z_ACC ---
        ##########################################
        prgMsg('Fitting Linear Data')

        # Fitted Model
        def fitFunc(x, a, b):
            return a * x + b

        params, cov = curve_fit(fitFunc, deltaVs, deltaTs)
        x_s = np.linspace(deltaVs.min(), deltaVs.max(), 20)
        fitData = [params[0] * x + params[1] for x in x_s]
        r_corr = round(np.corrcoef(deltaVs, deltaTs)[0, 1], 3)
        Done(start_time)

        #####################################
        # --- PLOT THE RESULTS OF THE FIT ---
        #####################################
        if outputCorrelationPlot:
            prgMsg('Outputting results of correlation analysis')

            fig, ax = plt.subplots()
            fig.suptitle(f'STEB {wDis}\n' + r'$\alpha$ =' + f'{Pitch[wPitch]} deg, maslVal = {maskVal}, N = {len(deltaVs)}')

            # Raw Data
            if showErrorBars:
                ax.errorbar(deltaVs, deltaTs, xerr=errorV, yerr=errorT, fmt='o', linewidth=2, capsize=6)
            else:
                ax.scatter(deltaVs, deltaTs)

            ax.set_ylabel('Delay time [s]')
            ax.set_xlabel('1/v [s/km]')
            ax.plot(x_s, fitData, color="red",label=f'd ={params[0]/6371}Re\n t_0={params[1]}\n r_corr = {r_corr}')
            xticks = np.linspace(0, deltaVs.max(), 6)
            xtick_labels = [f'{x:.2e}' for x in xticks]
            xtick_labels[0] = '0'
            ax.set_xticks(xticks,xtick_labels)
            plt.legend()
            plt.savefig(rf'C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\Science\AlfvenSingatureAnalysis\Particles\CrossCorrelationResults\STEB_{wDis}.png')
            Done(start_time)


        ########################################
        # --- Return results of the Analysis ---
        ########################################
        # return an array of: [STEB Number,Observation Time, Observation Altitude, Z_acc, Z_acc_error, correlationR]
        errorZ_avg = sum(errorZ) / len(errorZ)
        return [wDis,whenSTEBoccured_time,whenSTEBoccured_Alt,params[0]/Re, errorZ_avg, r_corr ]


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    STEBfitResults = []

    if len(wDispersions) == 1:
        results = AlfvenSignatureCrossCorrelation(wRocket, rocketFolderPath, justPrintFileNames, wDispersions[0])
        STEBfitResults.append(results)
    elif wDispersions == []:
        for i in range(len(dispersionAttributes.keyDispersionTimes)):
            results = AlfvenSignatureCrossCorrelation(wRocket, rocketFolderPath, justPrintFileNames, i+1)
            STEBfitResults.append(results)
    else:
        for i in range(len(wDispersions)):
            results = AlfvenSignatureCrossCorrelation(wRocket, rocketFolderPath, justPrintFileNames, wDispersions[i])
            STEBfitResults.append(results)


    if plotZsourceVsTime:

        # Load the data ESA... again
        data_dict = loadDictFromFile(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}{modifier}\*eepaa_fullcal*')[0])
        Energy = np.array(data_dict['Energy'][0])
        Pitch = data_dict['Pitch_Angle'][0]

        # Load the attitude data
        data_dict_attitude = loadDictFromFile(r'C:\Data\ACESII\attitude\high\ACESII_36359_Attitude_Solution.cdf')

        ###########################
        # --- Organize the data ---
        ###########################
        STEBfitResults = np.array(STEBfitResults)

        # FORMAT of STEB fit results:
        # [[STEB Number,Observation Time, Observation Altitude, Z_acc, Z_acc_error_avg, correlationR],
        #  [...],
        #  ]
        STEBtimes = STEBfitResults[:,1]
        STEB_geoLats = [data_dict_attitude['Lat_geom'][0][np.abs(data_dict_attitude['Epoch'][0] - tme).argmin()] for tme in STEBtimes]
        STEB_Zacc = (STEBfitResults[:,2]/1000)/Re + STEBfitResults[:,3]
        STEB_Zacc_avg = sum(STEB_Zacc)/len(STEB_Zacc)
        STEB_errorZacc_avg = STEBfitResults[:,4]

        # Find deltaE
        STEB_deltaE = np.array([Energy[np.abs(Energy - engyPair[1]).argmin()] - Energy[np.abs(Energy - engyPair[0]).argmin()] for engyPair in dispersionAttributes.keyDispersionDeltaE])

        # Find deltaT
        STEB_deltaT = np.array([pycdf.lib.datetime_to_tt2000(tmePair[1]) - pycdf.lib.datetime_to_tt2000(tmePair[0]) for tmePair in dispersionAttributes.keyDispersionDeltaT])/1E9

        # Find the labels to use for each STEB
        STEB_text_labels = []
        colorChoice = ['tab:orange', 'tab:purple', 'tab:cyan']
        # colorChoice = ['black', 'black', 'black']
        for type in dispersionAttributes.keyDispersionType:
            if 'isolat' in type.lower():
                STEB_text_labels.append(['I', colorChoice[0]])
            elif 'aurora' in type.lower():
                STEB_text_labels.append(['A', colorChoice[1]])
            elif 'edge' in type.lower():
                STEB_text_labels.append(['E', colorChoice[2]])

        STEB_text_labels = np.array(STEB_text_labels)

        # if not doing every steb, then reduce these
        if wDispersions != []:
            STEB_deltaT = STEB_deltaT[np.array(wDispersions) - 1]
            STEB_deltaE = STEB_deltaE[np.array(wDispersions) - 1]
            STEB_text_labels = STEB_text_labels[np.array(wDispersions) - 1]

        ##########################
        # --- PLOT THE RESULTS ---
        ##########################
        prgMsg('Outputting results of correlation analysis Zacc')
        fig, ax = plt.subplots(nrows=2)
        figure_width = 20
        figure_height = 10
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)

        # ---  toggles ---
        scatterSize = 250
        axisLabelFont = 25
        plotTextSize = 15

        # Raw Data
        if showErrorBars:
            ax[0].errorbar(STEBtimes, STEB_Zacc, yerr=STEB_errorZacc_avg, fmt='o', linewidth=2, capsize=6, color=scatColors)
        else:
            cmap = ax[0].scatter(x=STEBtimes, y=STEB_Zacc, s =scatterSize, c=STEB_text_labels[:,1])

        # plot the average z_acc
        ax[0].axhline(STEB_Zacc_avg, linestyle='--', color='black', linewidth=3)
        ax[0].text(dt.datetime(2022,11,20,17,25,45), STEB_Zacc_avg,f"Average: {round(STEB_Zacc_avg,2)} $R_E \pm ??$ ",fontsize=20,color='black', va='bottom')

        # Apply the text to each datapoint
        # for k in range(len(STEBtimes)):
        #     ax.text(STEBtimes[k], STEB_Zacc[k], STEB_text_labels[k][0], color=STEB_text_labels[k][1], fontsize='large', va='top', ha='left')

        # Plot Labels
        ax[0].set_ylabel('$Z_{acc}$ Altitude [$R_{E}$]', color='black', fontsize=axisLabelFont,labelpad=20)
        from matplotlib.dates import date2num,num2date
        startPoint, endPoint = date2num(min(STEBtimes)), date2num(max(STEBtimes))
        EpochTicks = num2date(np.linspace(startPoint,endPoint,8))
        offSetAwareAttitudeEpoch = date2num(data_dict_attitude['Epoch'][0])
        geomagLatTicks = [round(data_dict_attitude['Lat_geom'][0][ np.abs(offSetAwareAttitudeEpoch - date2num(tme)).argmin()],2) for tme in EpochTicks]
        ax[0].set_xticks(ticks=EpochTicks,labels=[ep.strftime('%H:%M.%S') for ep in EpochTicks])
        ax[0].minorticks_on()
        ax[0].tick_params(axis='both', which='major', labelsize=21,width=3,length=16)
        ax[0].tick_params(axis='both', which='minor', labelsize=21,width=3,length=8)
        ax[0].grid(True)
        ax[0].set_xticks(ticks=EpochTicks,labels=['' for thing in EpochTicks])

        for l in range(len(STEBtimes)):
            ax[0].text(x=STEBtimes[l], y=STEB_Zacc[l],s=f'{l+1}',va='center', ha='center',weight='bold',fontsize=plotTextSize)


        # --- create geomagnetic axis ---
        axGeoM = ax[0].twiny()
        axGeoM.scatter(STEBtimes,STEB_Zacc,alpha=0)
        axGeoM.set_xlabel('Geomagnetic Lattitude [deg]', fontsize=axisLabelFont, labelpad=20)

        # set the ticks for the geomagnetic top plot to match the epoch
        axGeoM.set_xticks(ticks=EpochTicks, labels=geomagLatTicks)

        axGeoM.minorticks_on()
        axGeoM.tick_params(axis='x', which='major', labelsize=21,width=3,length=16)
        axGeoM.tick_params(axis='x', which='minor', labelsize=21,width=3,length=8)

        # create my custom legend for Z_acc vs Epoch plot
        from matplotlib.lines import Line2D

        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor=colorChoice[0], label='Isolated', markersize=18),
            Line2D([0], [0], marker='o', color='w', markerfacecolor=colorChoice[1], label='Auroral', markersize=18),
            Line2D([0], [0], marker='o', color='w', markerfacecolor=colorChoice[2], label='Edge-Type', markersize=18)
        ]

        ax[0].legend(handles=legend_elements, loc='upper right', prop={'size': 15})


        # create delta E vs Epoch plot where I color the points in the ranges

        # determine the colors
        dT_colors = ['tab:blue','tab:green','tab:red']
        dT_color_labels = []
        for delT in STEB_deltaT:
            if delT <= 0.5:
                dT_color_labels.append(dT_colors[0])
            elif 0.5 < delT <= 1:
                dT_color_labels.append(dT_colors[1])
            elif delT > 1:
                dT_color_labels.append(dT_colors[2])

        ax[1].scatter(x=STEBtimes,y=STEB_deltaE,c=dT_color_labels,s=scatterSize)
        ax[1].grid(True)
        ax[1].set_xticks(ticks=EpochTicks, labels=[ep.strftime('%H:%M.%S') for ep in EpochTicks])
        ax[1].minorticks_on()
        ax[1].tick_params(axis='both', which='major', labelsize=21, width=3, length=16)
        ax[1].tick_params(axis='both', which='minor', labelsize=21, width=3, length=8)
        ax[1].set_ylabel('STEB $\Delta E$ [eV]',fontsize=axisLabelFont, labelpad=10)
        ax[1].set_xlabel('Epoch', color='black', fontsize=axisLabelFont, labelpad=10)
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor=dT_colors[0], label='$\Delta T \leq$ 0.5s', markersize=18),
            Line2D([0], [0], marker='o', color='w', markerfacecolor=dT_colors[1], label='0.5 < $\Delta T \leq$ 1s', markersize=18),
            Line2D([0], [0], marker='o', color='w', markerfacecolor=dT_colors[2], label='$\Delta T$ > 1s', markersize=18)
        ]

        for l in range(len(STEBtimes)):
            ax[1].text(x=STEBtimes[l], y=STEB_deltaE[l],s=f'{l+1}',va='center',ha='center',weight='bold',fontsize=plotTextSize)

        ax[1].legend(handles=legend_elements, loc='upper right', prop={'size': 15})

        # Show the figure
        # plt.subplots_adjust(wspace=0, hspace=1)
        plt.tight_layout()
        plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot5\STEB_CorrelationAnalysis_Results.png')
        Done(start_time)