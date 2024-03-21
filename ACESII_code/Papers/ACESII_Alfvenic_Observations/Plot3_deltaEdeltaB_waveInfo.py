# --- Plot3_detlaEdetlaB_waveInfo.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing field-aligned
# particle data along with electric and magnetic signatures. Also calculates the Local Alfven Speed
# Where the rocket is


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from numpy.fft import rfft, fftfreq
from ACESII_code.class_var_func import u0
from scipy.integrate import simpson

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
# figure_width = 7.5 # in inches
figure_width = 9.5 # in inches
figure_height = 11 # in inches
dpi = 500
PlotLabelSize = 10
PlotTitleSize = 10
PlotLegendSize = 10
PlotTickLabelSize = 10
freqLimit = 15
plotLineWidth = 1
TitlePadding = 10
EB_ratio_limits = [9E4, 1E8]
plotColors= ['tab:blue', 'tab:red', 'tab:orange', 'tab:green']
RegionNames = ['Quiet Region', 'Temporal Align', 'ILat Align']
kSets = [0, 1] # determines the E/B pairs: (1) B_e/E_r (2) B_r/-E_e or [0,1] for both
NORMALIZE = False

CalcPoyntingFluxEnergies = False


# --- TARGET REDUCTION VARIABLES ---
targetTimes = [[dt.datetime(2022, 11, 20, 17, 24, 31, 500), dt.datetime(2022, 11, 20, 17, 24, 38, 00)], # Quiet Time
               [dt.datetime(2022, 11, 20, 17, 24, 54, 250000), dt.datetime(2022, 11, 20, 17, 25, 3, 250000)], # Dispersed Time
               [dt.datetime(2022, 11, 20, 17, 25, 18, 000000), dt.datetime(2022, 11, 20, 17, 25, 27, 000000)]  # Auroral Time ( i.e. the FULL aurora)
               ]
# keeper - nonNormalized
# targetILats = [[71.255, 71.35], # Quiet Time
#                [71.54, 71.679], # temporally Time
#                [71.87, 71.945]] # spatially aligned

# trying to make it fit
targetILats = [[71.25, 71.35], # Quiet Time
               [71.54, 71.679], # temporally Time
               [71.859, 71.945]] # spatially aligned


targetVar = targetILats
targetVarName = 'ILat'

# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---

print(color.UNDERLINE + f'Plot3_WaveAnalysis' + color.END)
prgMsg('Loading Data')

inputFile_deltaB = glob('C:\Data\ACESII\L3\deltaB\low\*Field_Aligned*')[0] # get the deltaB data
inputFile_deltaE = glob('C:\Data\ACESII\L3\deltaE\low\*Field_Aligned*')[0] # get the deltaE data
inputFile_B = 'C:\Data\ACESII\L2\low\ACESII_36364_l2_RingCore_Field_Aligned.cdf'  # get the B data
inputFile_E = 'C:\Data\ACESII\L2\low\ACESII_36364_l2_E_Field_Field_Aligned_downsampled.cdf'  # get the E data
inputFile_poynting = glob('C:\Data\ACESII\science\PoyntingFlux\low\*Field_Aligned*')[0] # get the Poynting Flux data
inputFile_ASpeed = glob('C:\Data\ACESII\science\AlfvenSpeed_rkt\low\ACESII_36364_AlfvenSpeed_flight.cdf')[0] # get the Alfven Speed data
inputFile_Langmuir = 'C:\Data\ACESII\L3\Langmuir\low\ACESII_36364_langmuir_fixed.cdf'
inputFile_Bmag = 'C:\Data\ACESII\L1\low\ACESII_36364_l1_RingCore_rktFrm.cdf'
# --- Break up DataDicts into targetTime sections ---

# FORMAT: [[deltaB, deltaE, poynting, AlfvenSpeedRatio], ...]
dict_sets = []
sectionTimeRange = []

for i in range(len(targetVar)):
    data_dict_E = loadDictFromFile(inputFile_E,targetVar=[targetVar[i], targetVarName])
    data_dict_B = loadDictFromFile(inputFile_B,targetVar=[targetVar[i], targetVarName])
    data_dict_deltaB = deepcopy(loadDictFromFile(inputFile_deltaB,targetVar=[targetVar[i], targetVarName]))
    data_dict_deltaE = deepcopy(loadDictFromFile(inputFile_deltaE,targetVar=[targetVar[i], targetVarName]))
    data_dict_poynting = deepcopy(loadDictFromFile(inputFile_poynting,targetVar=[targetVar[i], targetVarName]))
    data_dict_langmuir = deepcopy(loadDictFromFile(inputFile_Langmuir,targetVar=[targetVar[i], targetVarName],wKeys_Load=['ni','Epoch','ILat']))

    # downsample the langmuir data
    indexVals = [np.abs(data_dict_langmuir['ILat'][0] - ILat).argmin() for ilt, ILat in enumerate(data_dict_B['ILat'][0])]
    data_dict_langmuir['ni'][0] = deepcopy(data_dict_langmuir['ni'][0][indexVals])
    data_dict_langmuir['Epoch'][0] = deepcopy(data_dict_langmuir['Epoch'][0][indexVals])
    data_dict_langmuir['ILat'][0] = deepcopy(data_dict_langmuir['ILat'][0][indexVals])
    sectionTimeRange.append([data_dict_deltaB['Epoch'][0][0], data_dict_deltaB['Epoch'][0][-1]])

    dict_sets.append([data_dict_deltaB,
                      data_dict_deltaE,
                      data_dict_poynting,
                      data_dict_langmuir,
                      data_dict_B,
                      data_dict_E
                      ])


    if CalcPoyntingFluxEnergies:
        print(f'\nCalculating Poynting Flux for {RegionNames[i]}')
        Sp = deepcopy(data_dict_poynting['S_p'][0])
        SILats = deepcopy(data_dict_poynting['ILat'][0])
        Sp_p = deepcopy(Sp)
        Sp_p[Sp_p <  0] = 0
        DownwardEnergy = simpson(Sp_p,x=SILats,dx=0.1)
        print(f'DownwardEnergy: {DownwardEnergy}')
        Sp_n = deepcopy(Sp)
        Sp_n[Sp_n > 0] = 0
        UpwardEnergy = simpson(Sp_n,x=SILats, dx=0.1)
        print(f'UpwardEnergy Area: {UpwardEnergy}')

Done(start_time)

#######################################
# --- PREPARE THE DATA FOR PLOTTING ---
#######################################
def Plot3_deltaEdeltaB_waveInfo(targetVar,dict_sets):

    # --- Calculate the FFT of deltaB, deltaE and their ratio ---
    # FORMAT: [[FFT_deltaBe, FFT_deltaBr,xf_deltaB, FFT_deltaEe, FFT_deltaEr,xf_DeltaE, FFT_ASpeed_ErBe, FFT_ASpeed_EeBr], ...]
    # Indicies:[          0,           1,        2,           3,           4,        5,               6,               7]
    data_dict_FFT = {
                'B0' : [],
                'Be_fft':[],
                'Br_fft':[],
                'Ee_fft':[],
                'Er_fft': [],
                'xf_deltaB':[],
                'xf_deltaE': [],
                'Be_norm_fft' : [],
                'Br_norm_fft' : [],
                'Ee_norm_fft' : [],
                'Er_norm_fft' : [],
                'ni' : [],
                'm_i_avg' : [],
                'VA_t' : [],
                'xf_B':[],
                'xf_E':[]}

    for i in range(len(targetVar)):
        print(i)

        data_dicts = dict_sets[i]

        # Calculate the normalized Quantities
        B_e = 1E-9 * data_dicts[4]['B_e'][0] # in tesla
        B_p = 1E-9 * data_dicts[4]['B_p'][0]
        B_r = 1E-9 * data_dicts[4]['B_r'][0]
        dB_e = data_dicts[0]['B_e'][0]  # in nanotesla
        dB_p = data_dicts[0]['B_p'][0]
        dB_r = data_dicts[0]['B_r'][0]
        B0 = 1E-9 * data_dicts[4]['Bmag'][0]

        E_e = data_dicts[5]['E_e'][0] # in V/m
        E_p = data_dicts[5]['E_p'][0]
        E_r = data_dicts[5]['E_r'][0]
        dE_e = data_dicts[1]['E_e'][0]  # in mV/m
        dE_p = data_dicts[1]['E_p'][0]
        dE_r = data_dicts[1]['E_r'][0]

        ni = (cm_to_m ** 3) * data_dicts[3]['ni'][0]
        m_i_avg = 2.45E-26
        VA_t = B0 / np.sqrt(u0 * ni * m_i_avg)

        # calculate the normalization
        Er_norm = E_r / (np.abs(VA_t) * np.abs(B0))
        Be_norm = B_e / (np.abs(B0))
        Ee_norm = E_e / (np.abs(VA_t) * np.abs(B0))
        Br_norm = B_r / (np.abs(B0))

        # calculate the FFTs

        # B_norm
        N, T = len(data_dicts[4]['Epoch'][0]), 1 / 128
        xf_B = fftfreq(N, T)[:N // 2]
        yf_Be = rfft(Be_norm)
        FFT_Be = 2.0 / N * np.abs(yf_Be[0:N // 2])
        yf_Br = rfft(Br_norm)
        FFT_Br = 2.0 / N * np.abs(yf_Br[0:N // 2])

        # deltaB
        N, T = len(data_dicts[0]['Epoch'][0]), 1 / 128
        xf_deltaB = fftfreq(N, T)[:N // 2]
        yf_dBe = rfft(dB_e)
        FFT_dBe = 2.0 / N * np.abs(yf_dBe[0:N // 2])
        yf_dBr = rfft(dB_r)
        FFT_dBr = 2.0 / N * np.abs(yf_dBr[0:N // 2])

        # E_norm
        N, T = len(data_dicts[1]['Epoch'][0]), 1 / 128
        xf_E = fftfreq(N, T)[:N // 2]
        yf_Ee = rfft(Ee_norm)
        FFT_Ee = 2.0 / N * np.abs(yf_Ee[0:N // 2])
        yf_Er = rfft(Er_norm)
        FFT_Er = 2.0 / N * np.abs(yf_Er[0:N // 2])

        # deltaE
        N, T = len(data_dicts[5]['Epoch'][0]), 1 / 128
        xf_deltaE = fftfreq(N, T)[:N // 2]
        yf_dEe = rfft(dE_e)
        FFT_dEe = 2.0 / N * np.abs(yf_dEe[0:N // 2])
        yf_dEr = rfft(dE_r)
        FFT_dEr = 2.0 / N * np.abs(yf_dEr[0:N // 2])

        # store everything
        data_dict_FFT['B0'].append(B0)
        data_dict_FFT['Be_fft'].append(FFT_dBe)
        data_dict_FFT['Br_fft'].append(FFT_dBr)
        data_dict_FFT['Ee_fft'].append(FFT_dEe)
        data_dict_FFT['Er_fft'].append(FFT_dEr)
        data_dict_FFT['xf_deltaB'].append(xf_deltaB)
        data_dict_FFT['xf_deltaE'].append(xf_deltaE)
        data_dict_FFT['Be_norm_fft'].append(FFT_Be)
        data_dict_FFT['Br_norm_fft'].append(FFT_Br)
        data_dict_FFT['Ee_norm_fft'].append(FFT_Ee)
        data_dict_FFT['Er_norm_fft'].append(FFT_Er)
        data_dict_FFT['xf_B'].append(xf_B)
        data_dict_FFT['xf_E'].append(xf_E)
        data_dict_FFT['ni'].append(ni)
        data_dict_FFT['m_i_avg'].append(m_i_avg)
        data_dict_FFT['VA_t'].append(VA_t)

    # ---------------------
    #######################
    # --- PLOT THE DATA ---
    #######################
    # ---------------------
    # prepare the data for plotting
    waveSetKeys = [
        ['B_e', 'E_r'],
        ['B_r', 'E_e']]

    waveSetLabels = [['$\delta B_{e}$', '$\delta E_{r}$'],
                     [' $\delta B_{r}$', '-$\delta E_{e}$']]

    if NORMALIZE:
        spectraFreqs = [data_dict_FFT['xf_B'], data_dict_FFT['xf_E']]

        spectraData = [
                       [data_dict_FFT['Be_norm_fft'], data_dict_FFT['Er_norm_fft']],
                       [data_dict_FFT['Br_norm_fft'], data_dict_FFT['Ee_norm_fft']]
                       ]
        spectraLabels = [['$B_{e}(t)/|B_{0}(t)|$', '$E_{r}(t)/|B_{0}(t)||V_{A}(t)|$'],
                         ['$B_{r}(t)/|B_{0}(t)|$', '$E_{e}(t)/|B_{0}(t)||V_{A}(t)|$']]

        alfvenData = [ [ np.array(spectraData[0][1][i])/np.array(spectraData[0][0][i])  for i in range(len(targetVar))],
                       [ np.array(spectraData[1][1][i])/np.array(spectraData[1][0][i]) for i in range(len(targetVar))],
                       ]

    else:
        spectraFreqs = [data_dict_FFT['xf_deltaB'], data_dict_FFT['xf_deltaE']]
        spectraData = [  [data_dict_FFT['Be_fft'], data_dict_FFT['Er_fft']],
                         [data_dict_FFT['Br_fft'], data_dict_FFT['Ee_fft']]
                         ]
        spectraLabels = [['$\delta B_{e}$ [nT]', '$\delta E_{r}$ [mV/m]'],
                         [' $\delta B_{r}$ [nT]', '-$\delta E_{e}$ [mV/m]']]
        alfvenData = [[(1E6)*np.array(spectraData[0][1][i]) / np.array(spectraData[0][0][i]) for i in range(len(targetVar))],
                      [(1E6)*np.array(spectraData[1][1][i]) / np.array(spectraData[1][0][i]) for i in range(len(targetVar))],
                      ]



    #############################
    ### DO THE PLOTTING LOOPS ###
    #############################

    for wKeySet in kSets:

        prgMsg('Plotting Data')
        fig, ax = plt.subplots(nrows=4, ncols=len(targetVar))
        fig.set_size_inches(figure_width, figure_height)

        wWaveSetKeys = waveSetKeys[wKeySet]
        wSpectraLabel = spectraLabels[wKeySet]
        wSpectraData = spectraData[wKeySet]
        wAlfvenData = alfvenData[wKeySet]
        wWaveSetLabels = waveSetLabels[wKeySet]


        ######## LOOP THROUGH ALL INVARIANT LATTITUDE REGIONS ########
        for i in range(len(targetVar)):

            data_dicts = dict_sets[i] # raw data

            # Set the title
            ax[0][i].set_title(RegionNames[i] + f'\n {sectionTimeRange[i][0].strftime("%H:%M:%S")} to {sectionTimeRange[i][1].strftime("%H:%M:%S")} UTC',
                               weight='bold',
                               fontsize=PlotTitleSize,
                               pad=TitlePadding)

            ####################
            # --- Everything ---
            #####################

            # set the grid
            for j in range(4):
                ax[j][i].grid(which='both', alpha=0.5)

            #################################
            # --- delta Be, delta Er plot ---
            #################################

            # delta Be
            ax[0][i].set_xlabel('ILat [deg]',fontsize=PlotLabelSize)
            ax[0][i].set_ylim(-8,8)
            ax[0][i].set_xmargin(0)
            ax[0][i].tick_params(axis='both', which='major', labelsize=PlotTickLabelSize)
            ax[0][i].tick_params(axis='y', which='minor', labelsize=PlotTickLabelSize,length=0,grid_alpha=0)
            ax[0][i].tick_params(axis='both', which='minor', labelsize=PlotTickLabelSize-2)
            ax[0][i].plot(data_dicts[0]['ILat'][0], data_dicts[0][wWaveSetKeys[0]][0], color=plotColors[0], linewidth=plotLineWidth, label=f'{wWaveSetLabels[0]} [nT]')

            # delta Er
            mod = -1 if wKeySet == 1 else 1
            ax[0][i].plot(data_dicts[1]['ILat'][0], mod*data_dicts[1][wWaveSetKeys[1]][0], color=plotColors[1], linewidth=plotLineWidth,label=f"{wWaveSetLabels[1]} [mV/m]")
            ax[0][i].set_xmargin(0)
            ax[0][i].minorticks_on()

            # fix the xticks
            newTicks = data_dicts[1]['ILat'][0][::int(len(data_dicts[1]['ILat'][0])/2)]
            newTicks = [round(tick,2) for tick in newTicks]
            newTicks.append(round(data_dicts[1]['ILat'][0][-1],2))
            newTickStr = [str(tick) for tick in newTicks]
            ax[0][i].set_xticks(newTicks, newTickStr)

            ############################
            # --- Poynting Flux plot ---
            ############################
            ax[1][i].plot(data_dicts[2]['ILat'][0], 1000*np.array(data_dicts[2]['S_p'][0]),plotColors[2], label='$\delta S_{p}$ [mW/$m^{2}$]',linewidth=plotLineWidth)
            ax[1][i].set_xlabel('ILat [deg]',fontsize=PlotLabelSize)
            ax[1][i].set_ylim(-2.2E-2, 2.8E-2)
            ax[1][i].set_xmargin(0)
            ax[1][i].invert_yaxis()
            ax[1][i].tick_params(axis='both', which='major', labelsize=PlotTickLabelSize)
            ax[1][i].tick_params(axis='both', which='minor', labelsize=PlotTickLabelSize-2)
            newTicks = data_dicts[1]['ILat'][0][:-1:int(len(data_dicts[1]['ILat'][0])/2)]
            newTicks = [round(tick, 2) for tick in newTicks]
            newTicks.append(round(data_dicts[1]['ILat'][0][-1], 2))
            newTickStr = [str(tick) for tick in newTicks]
            ax[1][i].set_xticks(newTicks, newTickStr)
            ax[1][i].minorticks_on()

            ######################
            # --- FFT E, FFT B ---
            ######################

            ln1 = ax[2][i].plot(spectraFreqs[0][i],wSpectraData[0][i], color=plotColors[0], label=f'{wSpectraLabel[0]}',linewidth=plotLineWidth)
            ax[2][i].set_xlabel('Frequency [Hz]',fontsize=PlotLabelSize)
            ax[2][i].set_xlim(0, freqLimit)
            ax[2][i].set_yscale('log')
            ax[2][i].set_ylim(1E-3,1E1)
            # ax[2][i].set_xscale('log')
            if NORMALIZE:
                ax[2][i].set_ylim(1E-7, 1E-2)
                # freqS = [val**(-1.67) for val in spectraFreqs[0][i] if val !=0]
                # ax[2][i].plot(spectraFreqs[0][i][:-1],freqS,color='black',linestyle='--',label='$f^{-1.67}$')
            else:
                # freqS = [val ** (-1.67) for val in spectraFreqs[0][i] if val != 0]
                # ax[2][i].plot(spectraFreqs[0][i][:-1], freqS, color='black', linestyle='--')
                for k in range(3):
                    ax[2][i].axvline(0.55*(k+1), color='green',linestyle='--',alpha=0.5)

            ax[2][i].tick_params(axis='both', which='major', labelsize=PlotTickLabelSize)
            ax[2][i].tick_params(axis='both', which='minor', labelsize=PlotTickLabelSize-2)

            axFFT_E = ax[2][i].twinx()
            ln2 = ax[2][i].plot(spectraFreqs[1][i], wSpectraData[1][i], color=plotColors[1], label=f'{wSpectraLabel[1]}',linewidth=plotLineWidth)
            axFFT_E.set_yscale('log')
            axFFT_E.set_ylim(1E-3, 1E1)
            axFFT_E.set_yticks([])




            ######################
            # --- ASpeed Ratio ---
            ######################
            ax[3][i].plot(spectraFreqs[0][i], wAlfvenData[i], color='black', label=f'', linewidth=0.75, linestyle='-', marker='o', markersize=1)
            # ax[3][i].plot(spectraFreqs[0][i], wAlfvenData[i], color='black', label=f'', linewidth=1.2, linestyle='-')


            if NORMALIZE:
                ax[3][i].axhline(1, linestyle='--', color='red', label='V$_{A}$ (Avg)')
                ax[3][i].set_ylim(0, 3)
            else:
                ax[3][i].axhline(sum(wAlfvenData[i])/len(wAlfvenData[i]),linestyle='--', color='red', label='V$_{A}$ (Avg) [m/s]')
                ax[3][i].set_ylim(EB_ratio_limits[0], EB_ratio_limits[1])
                ax[3][i].set_yscale('log')

                for k in range(3):
                    ax[3][i].axvline(0.55*(k+1), color='green',linestyle='--',alpha=0.5)


            ax[3][i].minorticks_on()
            ax[3][i].set_xlabel('Frequency [Hz]', fontsize=PlotLabelSize)
            ax[3][i].set_xlim(0, freqLimit)
            ax[3][i].tick_params(axis='both', which='major', labelsize=PlotTickLabelSize)
            ax[3][i].tick_params(axis='both', which='minor', labelsize=PlotTickLabelSize-2)


            # set the legend(s)
            lns = ln1 + ln2
            labs = [l.get_label() for l in lns]
            if i == 0:
                ax[0][i].set_ylabel(f'{wWaveSetLabels[0]} and {wWaveSetLabels[1]}' + "\n 0.4 to 20 Hz", fontsize=PlotLabelSize)
                ax[0][i].legend(loc='upper right', prop={'size': PlotLegendSize})
                ax[1][i].set_ylabel('Field-Aligned\n Poynting Flux', fontsize=PlotLabelSize)
                ax[1][i].legend(loc='upper right', prop={'size': PlotLegendSize})
                ax[2][i].set_ylabel('Amplitude Spectra', fontsize=PlotLabelSize)
                axFFT_E.legend(lns, labs, loc='upper right', prop={'size': PlotLegendSize})
                ax[3][i].legend(loc='upper right', prop={'size': PlotLegendSize})

                if NORMALIZE:
                    ax[3][i].set_ylabel('$E_{\perp}(t)/B_{\perp}(t)$', fontsize=PlotLabelSize)
                else:
                    ax[3][i].set_ylabel(f'{wWaveSetLabels[1]}/{wWaveSetLabels[0]}', fontsize=PlotLabelSize)


        ### OUTPUT THE DATA ###
        plt.tight_layout()
        fileOutName = "Plot3_WaveAnalysis_BeEr.png" if wWaveSetKeys[0] == 'B_e' and wWaveSetKeys[1] =='E_r' else "Plot3_WaveAnalysis_mBrEe.png"
        outputPath = rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot3\{fileOutName}'
        plt.savefig(outputPath, dpi=dpi)
        Done(start_time)



#################
# --- EXECUTE ---
#################
Plot3_deltaEdeltaB_waveInfo(targetVar,dict_sets)
