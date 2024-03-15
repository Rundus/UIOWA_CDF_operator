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

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
# figure_width = 7.5 # in inches
figure_width = 8.5 # in inches
figure_height = 11 # in inches
dpi = 300
PlotLabelSize = 13
PlotTitleSize = 11.5
PlotLegendSize = 10
PlotTickSize = 10
freqLimit = 15
plotLineWidth = 1
TitlePadding = 15
EB_ratio_limits = [9E4,1E8]
plotColors= ['tab:blue','tab:red', 'tab:orange','tab:green']
RegionNames = ['Quiet Region', 'Temporal Align', 'ILat Align']
kSets = [0,1] # determines the E/B pairs: (1) B_e/E_r (2) B_r/-E_e or [0,1] for both
highFlyerShennanigans = False
E_Field_scale_Factor = 0.6



# --- TARGET REDUCTION VARIABLES ---
targetTimes = [[dt.datetime(2022,11,20,17,24,31,500), dt.datetime(2022,11,20,17, 24, 38, 00)], # Quiet Time
               [dt.datetime(2022,11,20,17,24,54,250000), dt.datetime(2022,11,20,17, 25, 3, 250000)], # Dispersed Time
               [dt.datetime(2022,11,20,17,25,18,000000), dt.datetime(2022,11,20,17, 25, 27, 000000)]  # Auroral Time ( i.e. the FULL aurora)
               ]
# keeper
# targetILats = [[71.25, 71.34], # Quiet Time
#                [71.581, 71.68], # Dispersed Time
#                [71.831, 71.922]]

# trying to make it fit
# targetILats = [[71.25, 71.34], # Quiet Time
#                [71.581, 71.68], # Dispersed Time
#                [71.85, 71.922]]

targetILats = [[71.2, 71.34], # Quiet Time
               [71.49, 71.68], # Dispersed Time
               [71.8, 71.94]]


# HF sshennanigans
# targetILats = [[71.25, 71.33],  # Quiet Time
#                    [71.4, 71.693],  # Dispersed Time
#                    [71.9, 72.06]]

targetVar = targetILats
targetVarName = 'ILat'

# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---

print(color.UNDERLINE + f'Plot3_WaveAnalysis' + color.END)
prgMsg('Loading Data')

if highFlyerShennanigans:
    inputFile_deltaB = glob('C:\Data\ACESII\L3\deltaB\high\*Field_Aligned*')[0]  # get the deltaB data
    inputFile_deltaE = glob('C:\Data\ACESII\L3\deltaE\low\*Field_Aligned*')[0]  # get the deltaE data
    inputFile_poynting = glob('C:\Data\ACESII\science\PoyntingFlux\low\*Field_Aligned*')[0]  # get the Poynting Flux data
    inputFile_Langmuir = 'C:\Data\ACESII\L3\Langmuir\high\ACESII_36359_langmuir_fixed.cdf'
    inputFile_Bmag = 'C:\Data\ACESII\L1\high\ACESII_36359_l1_RingCore_rktFrm.cdf'
else:
    inputFile_deltaB = glob('C:\Data\ACESII\L3\deltaB\low\*Field_Aligned*')[0] # get the deltaB data
    inputFile_deltaE = glob('C:\Data\ACESII\L3\deltaE\low\*Field_Aligned*')[0] # get the deltaE data
    inputFile_poynting = glob('C:\Data\ACESII\science\PoyntingFlux\low\*Field_Aligned*')[0] # get the Poynting Flux data
    inputFile_ASpeed = glob('C:\Data\ACESII\science\AlfvenSpeed_rkt\low\ACESII_36364_AlfvenSpeed_flight.cdf')[0] # get the Alfven Speed data
    inputFile_Langmuir = 'C:\Data\ACESII\L3\Langmuir\low\ACESII_36364_langmuir_fixed.cdf'
    inputFile_Bmag = 'C:\Data\ACESII\L1\low\ACESII_36364_l1_RingCore_rktFrm.cdf'
# --- Break up DataDicts into targetTime sections ---

# FORMAT: [[deltaB, deltaE, poynting, AlfvenSpeedRatio], ...]
dict_sets = []
sectionTimeRange = []

for i in range(len(targetVar)):
    data_dict_deltaB = deepcopy(loadDictFromFile(inputFile_deltaB,targetVar=[targetVar[i],targetVarName]))
    data_dict_deltaE = deepcopy(loadDictFromFile(inputFile_deltaE,targetVar=[targetVar[i],targetVarName]))
    data_dict_poynting = deepcopy(loadDictFromFile(inputFile_poynting,targetVar=[targetVar[i],targetVarName]))
    data_dict_langmuir = deepcopy(loadDictFromFile(inputFile_Langmuir,targetVar=[targetVar[i],targetVarName],wKeys_Load=['ni','ni_error','Epoch','ILat']))
    data_dict_Bmag = deepcopy(loadDictFromFile(inputFile_Bmag,targetVar=[targetVar[i],targetVarName]))
    sectionTimeRange.append([data_dict_deltaB['Epoch'][0][0],data_dict_deltaB['Epoch'][0][-1]])

    # downsample the langmuir data



    if highFlyerShennanigans:
        data_dict_deltaE['E_r'][0] = deepcopy(E_Field_scale_Factor * data_dict_deltaB['B_e'][0])
        data_dict_deltaE['E_e'][0] = -1*deepcopy(E_Field_scale_Factor * data_dict_deltaB['B_r'][0])
        data_dict_deltaE['E_p'][0] = deepcopy(E_Field_scale_Factor * data_dict_deltaB['B_p'][0])
        data_dict_deltaE['Epoch'][0] = deepcopy(data_dict_deltaB['Epoch'][0])
        data_dict_deltaE['ILat'][0] = deepcopy(data_dict_deltaB['ILat'][0])

        Evec = np.array([  [data_dict_deltaE['E_r'][0][i],data_dict_deltaE['E_e'][0][i],data_dict_deltaE['E_p'][0][i]]   for i in range(len(data_dict_deltaE['Epoch'][0]))])
        Bvec = np.array([  [data_dict_deltaB['B_r'][0][i],data_dict_deltaB['B_e'][0][i],data_dict_deltaB['B_p'][0][i]]   for i in range(len(data_dict_deltaB['Epoch'][0]))])
        S = (1E-9 * 1E-3)*np.array([np.cross(Evec[i],Bvec[i])/u0 for i in range(len(Evec))])
        data_dict_poynting['S_r'][0] = S[:, 0]
        data_dict_poynting['S_e'][0] = S[:, 1]
        data_dict_poynting['S_p'][0] = S[:, 2]
        data_dict_poynting['ILat'][0] = deepcopy(data_dict_deltaB['ILat'][0])
        data_dict_poynting['Epoch'][0] = deepcopy(data_dict_deltaB['Epoch'][0])

    dict_sets.append([data_dict_deltaB,
                      data_dict_deltaE,
                      data_dict_poynting,
                      data_dict_langmuir,
                      data_dict_Bmag
                      ])
Done(start_time)


#######################################
# --- PREPARE THE DATA FOR PLOTTING ---
#######################################
def Plot3_deltaEdeltaB_waveInfo(targetVar,dict_sets,highFlyerShennanigans):

    # --- Calculate the FFT of deltaB, deltaE and their ratio ---
    # FORMAT: [[FFT_deltaBe, FFT_deltaBr,xf_deltaB, FFT_deltaEe, FFT_deltaEr,xf_DeltaE, FFT_ASpeed_ErBe, FFT_ASpeed_EeBr], ...]
    # Indicies:[          0,           1,        2,           3,           4,        5,               6,               7]
    FFT_sets = []

    for i in range(len(targetVar)):
        print(i)

        data_dicts = dict_sets[i]
        scale = 1 # 1.363636

        # deltaB
        wDict = 0
        N, T = len(data_dicts[wDict]['Epoch'][0]), 1 / 128
        xf_deltaB = fftfreq(N, T)[:N // 2]
        yf_Be = rfft(data_dicts[wDict]['B_e'][0])
        FFT_Be = 2.0 / N * np.abs(yf_Be[0:N // 2]) if i in [0,1] else scale*2.0 / N * np.abs(yf_Be[0:N // 2])
        yf_Br = rfft(data_dicts[wDict]['B_r'][0])
        FFT_Br = 2.0 / N * np.abs(yf_Br[0:N // 2]) if i in [0,1] else scale*2.0 / N * np.abs(yf_Br[0:N // 2])

        # deltaE
        wDict = 1
        N, T = len(data_dicts[wDict]['Epoch'][0]), 1 / 128
        xf_deltaE = fftfreq(N, T)[:N // 2]
        yf_Ee = rfft(data_dicts[wDict]['E_e'][0])
        FFT_Ee = 2.0 / N * np.abs(yf_Ee[0:N // 2])
        yf_Er = rfft(data_dicts[wDict]['E_r'][0])
        FFT_Er = 2.0 / N * np.abs(yf_Er[0:N // 2])

        # ASpeed
        ErBe_ratio = ((1E-3*FFT_Er)/(FFT_Be*1E-9))
        mEeBr_ratio = ((1E-3*FFT_Ee)/(FFT_Br*1E-9))

        # Average MHD Alfven Speed
        ni_avg = (cm_to_m ** 3) * sum(data_dicts[3]['ni'][0]) / len(data_dicts[3]['ni'][0])
        ni_error_avg =  (cm_to_m ** 3) *sum(data_dicts[3]['ni_error'][0]) / len(data_dicts[3]['ni_error'][0])
        Bmag_avg = (1E-9) * sum(data_dicts[4]['Bmag'][0]) / len(data_dicts[4]['Bmag'][0])
        # m_i_avg = IonMasses[1]
        m_i_avg = 2.45E-26
        MHD_Alfven_avg = Bmag_avg / np.sqrt(u0 * ni_avg * m_i_avg)


        # normalization temp plotting
        Be = data_dicts[0]['B_e'][0]
        Br = data_dicts[0]['B_r'][0]
        Bp = data_dicts[0]['B_p'][0]

        Ee = data_dicts[1]['E_e'][0]
        Er = data_dicts[1]['E_r'][0]
        Ep = data_dicts[1]['E_p'][0]

        ni = data_dicts[3]['ni'][0]
        Bmag = 1E-9 * data_dicts[4]['Bmag'][0]

        print(len(Be),len(Ee),len(ni),len(Bmag))

        fig, ax = plt.subplots(3)




        # --- error analysis ---
        # deltaB_error = 0.5E-9
        # deltaN_error = ni_error_avg
        # print(f'ni_avg: {ni_avg}')
        # print(f'ni_error_avg: {ni_error_avg}')
        # print(f'ni_error/n_i: {ni_error_avg/ni_avg}')
        # print(f'Alfv_Speed_avg: {MHD_Alfven_avg}')
        # dVdB = (u0*ni_avg*m_i_avg)**(-0.5)
        # dVdni= 0.5*(ni_avg**(-1.5)) * (Bmag_avg/np.sqrt(u0*m_i_avg))
        # MHD_Alfven_error = np.sqrt((deltaB_error*dVdB)**2 + (deltaN_error*dVdni)**2)
        # print(f'Alfv_Speed_error: {MHD_Alfven_error}')
        # print(f'Speed_error/Speed: {MHD_Alfven_error/MHD_Alfven_avg}')

        # store everything
        FFT_sets.append([FFT_Be, FFT_Br, xf_deltaB, FFT_Ee, FFT_Er, xf_deltaE, ErBe_ratio, mEeBr_ratio, MHD_Alfven_avg])

    # ---------------------
    #######################
    # --- PLOT THE DATA ---
    #######################
    # ---------------------

    for wKeySet in kSets:

        prgMsg('Plotting Data')

        # for each of the three regions, plot the following:
        # (1) deltaBe, deltaEr (Different yaxis)
        # (2) Poynting FLux
        # (3) FFT E, FFT B
        # (4) ASpeed Ratio


        fig, ax = plt.subplots(nrows=4, ncols=len(targetVar))
        # plt.subplots_adjust(hspace=0.35, wspace=0.32)
        # fig.suptitle('Low Flyer E and B Data',fontsize=30)
        fig.set_size_inches(figure_width,figure_height)

        keySet = [
            [['B_e','E_r'],['$\delta B_{e}$','$\delta E_{r}$']],
            [['B_r','E_e'],[' $\delta B_{r}$','-$\delta E_{e}$']],
        ]

        for i in range(len(targetVar)):

            data_dicts = dict_sets[i] # raw data
            FFT_data = FFT_sets[i]

            # Set the title
            ax[0][i].set_title(RegionNames[i] + f'\n {sectionTimeRange[i][0].strftime("%H:%M:%S")} to {sectionTimeRange[i][1].strftime("%H:%M:%S")} (UTC)',
                               weight='bold',
                               fontsize=PlotTitleSize,
                               pad=TitlePadding)

            ####################
            # --- Everything ---
            ######################

            # set the grid
            for j in range(4):
                ax[j][i].grid(which='both', alpha=0.25)

            #################################
            # --- delta Be, delta Er plot ---
            #################################
            zorder = [[0,1],[0,1],[0,1]]

            # delta Be
            ax[0][i].set_xlabel('ILat [deg]',fontsize=PlotLabelSize)
            ax[0][i].set_ylim(-13, 13)
            ax[0][i].set_xmargin(0)
            ax[0][i].tick_params(axis='both', which='major', labelsize=PlotTickSize)
            ax[0][i].tick_params(axis='both', which='minor', labelsize=PlotTickSize-2)
            if i == 0:
                ax[0][i].set_ylabel(f'{keySet[wKeySet][0][0]} and {keySet[wKeySet][0][1]}' + "\n 0.4 to 20 Hz",fontsize=PlotLabelSize)
            ln1 = ax[0][i].plot(data_dicts[0]['ILat'][0], data_dicts[0][keySet[wKeySet][0][0]][0], color=plotColors[0], linewidth=plotLineWidth, label=f'{keySet[wKeySet][1][0]} [nT]', zorder=zorder[i][0])

            # delta Er
            plotThisEData = data_dicts[1][keySet[wKeySet][0][1]][0] if keySet[wKeySet][0][1] == 'E_r' else -1 * np.array( data_dicts[1][keySet[wKeySet][0][1]][0])
            ax[0][i].plot(data_dicts[1]['ILat'][0], plotThisEData, color=plotColors[1], linewidth=plotLineWidth,label=f"{keySet[wKeySet][1][1]} [mV/m]", zorder=zorder[i][1])
            ax[0][i].set_xmargin(0)

            #set the legend
            if i == 0:
                ax[0][i].legend(loc='upper right', prop={'size':PlotLegendSize})

            # fix the xticks
            newTicks = data_dicts[1]['ILat'][0][::int(len(data_dicts[1]['ILat'][0])/2)]
            newTicks = [round(tick,2) for tick in newTicks]
            newTickStr = [str(tick) for tick in newTicks]
            ax[0][i].set_xticks(newTicks,newTickStr)



            ############################
            # --- Poynting Flux plot ---
            ############################
            ax[1][i].plot(data_dicts[2]['ILat'][0], 1000*np.array(data_dicts[2]['S_p'][0]),plotColors[2], label='$\delta S_{p}$ [mW/$m^{2}$]',linewidth=plotLineWidth)
            if i == 0:
                ax[1][i].set_ylabel('Field-Aligned Poynting Flux',fontsize=PlotLabelSize)
            ax[1][i].set_xlabel('ILat [deg]',fontsize=PlotLabelSize)

            ax[1][i].set_ylim(-2E-2, 3E-2)

            ax[1][i].set_xmargin(0)
            if i == 0:
                ax[1][i].legend(loc='upper right',prop={'size': PlotLegendSize})
            ax[1][i].invert_yaxis()
            ax[1][i].tick_params(axis='both', which='major', labelsize=PlotTickSize)
            ax[1][i].tick_params(axis='both', which='minor', labelsize=PlotTickSize-2)
            newTicks = data_dicts[2]['ILat'][0][:-1:int(len(data_dicts[1]['ILat'][0])/2)]
            newTicks = [round(tick, 2) for tick in newTicks]
            newTickStr = [str(tick) for tick in newTicks]
            ax[1][i].set_xticks(newTicks, newTickStr)

            ######################
            # --- FFT E, FFT B ---
            ######################
            plotThisFFT_B_data = FFT_data[0] if keySet[wKeySet][0][0] == 'B_e' and keySet[wKeySet][0][1]=='E_r' else FFT_data[1]
            plotThisFFT_E_data = FFT_data[3] if keySet[wKeySet][0][0] == 'B_r' and keySet[wKeySet][0][1]=='E_e' else FFT_data[4]
            ln1 = ax[2][i].plot(FFT_data[2],plotThisFFT_B_data , color=plotColors[0], label=f'{keySet[wKeySet][1][0]} [nT]',linewidth=plotLineWidth)
            if i == 0:
                ax[2][i].set_ylabel('Amplitude Spectra',fontsize=PlotLabelSize)
            ax[2][i].set_xlabel('Frequency [Hz]',fontsize=PlotLabelSize)
            ax[2][i].set_xlim(0, freqLimit)
            ax[2][i].set_yscale('log')
            if not highFlyerShennanigans:
                ax[2][i].set_ylim(1E-3, 1E1)
            ax[2][i].tick_params(axis='both', which='major', labelsize=PlotTickSize)
            ax[2][i].tick_params(axis='both', which='minor', labelsize=PlotTickSize-2)
            ax[2][i].minorticks_on()

            axFFT_E = ax[2][i].twinx()
            ln2 = ax[2][i].plot(FFT_data[2], plotThisFFT_E_data, color=plotColors[1],label=f"{keySet[wKeySet][1][1]} [mV/m]", linewidth=plotLineWidth)
            axFFT_E.set_yscale('log')
            if not highFlyerShennanigans:
                axFFT_E.set_ylim(1E-3, 1E1)
            axFFT_E.set_yticks([])

            # set the legend
            lns = ln1 + ln2
            labs = [l.get_label() for l in lns]
            if i == 0:
                axFFT_E.legend(lns, labs, loc='upper right',prop={'size': PlotLegendSize})

            ######################
            # --- ASpeed Ratio ---
            ######################

            if keySet[wKeySet][0][0] == 'B_e' and keySet[wKeySet][0][1] == 'E_r':
                FFT_data_plot_this = FFT_data[6]
            elif keySet[wKeySet][0][0] == 'B_r' and keySet[wKeySet][0][1] == 'E_e':
                FFT_data_plot_this = FFT_data[7]

            ax[3][i].plot(FFT_data[2], FFT_data_plot_this,color='black',label=f'{keySet[wKeySet][1][1]}/{keySet[wKeySet][1][0]}',linewidth=0.5,linestyle='-',marker='o',markersize=1)

            if i == 0:
                ax[3][i].set_ylabel('E/B Ratio [m/s]',fontsize=PlotLabelSize)

            # plot the Alfven Speed Calculated from Langmuir Data
            for k in range(3):
                ax[3][i].axvline(0.55*(k+1),linestyle='--', color='green')

            ax[3][i].axhline(FFT_data[8],linestyle='--', color='red', label='V$_{A}$ (MHD)')
            ax[3][i].minorticks_on()
            ax[3][i].set_xlabel('Frequency [Hz]', fontsize=PlotLabelSize)
            ax[3][i].set_xlim(0, freqLimit)
            if not highFlyerShennanigans:
                ax[3][i].set_ylim(EB_ratio_limits[0], EB_ratio_limits[1])
            ax[3][i].set_yscale('log')
            ax[3][i].tick_params(axis='both', which='major', labelsize=PlotTickSize)
            ax[3][i].tick_params(axis='both', which='minor', labelsize=PlotTickSize-2)
            if i == 0:
                ax[3][i].legend(loc='upper right',prop={'size': PlotLegendSize})


            # # Plot reflection function
            # freq = np.linspace(min(FFT_data[2]),max(FFT_data[2]),100)
            # print(FFT_data[8])
            # def IAR(SimgaP,freq, Alt):
            #     VA = FFT_data[8]
            #     SigmaA = 1 /(u0*VA)
            #     R = (SigmaA - SimgaP)/(SigmaA + SimgaP)
            #     return VA * (1 - R*np.cos(-2*freq*Alt/VA)) / (1+R*np.cos(-2*freq*Alt/VA))
            #
            # IAR_Ratio = [IAR(10,val,380E3) for val in freq]
            # ax[3][i].plot(freq,IAR_Ratio, color='orange')


        plt.tight_layout()
        fileOutName = "Plot3_WaveAnalysis_BeEr.png" if keySet[wKeySet][0][0] == 'B_e' and keySet[wKeySet][0][1]=='E_r' else "Plot3_WaveAnalysis_mBrEe.png"
        outputPath = rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot3\HF_shennanigans\{fileOutName}' if highFlyerShennanigans else rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot3\{fileOutName}'
        plt.savefig(outputPath,dpi=dpi)
        Done(start_time)



#################
# --- EXECUTE ---
#################
Plot3_deltaEdeltaB_waveInfo(targetVar,dict_sets,highFlyerShennanigans)
