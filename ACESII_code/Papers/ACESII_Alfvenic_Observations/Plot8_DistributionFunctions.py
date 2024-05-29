# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing field-aligned
# particle data along with electric and magnetic signatures

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt

from ACESII_code.myImports import *
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from ACESII_code.class_var_func import EpochTo_T0_Rocket, q0,m_e, Re
import matplotlib.gridspec as gridspec
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
print(color.UNDERLINE + f'Plot8_Conjugacy' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# --- Physics Toggles ---
datasetReduction_TargetTime = [dt.datetime(2022,11, 20, 17, 24, 23, 000000), dt.datetime(2022,11,20,17,25,10,000000)]
targetVar = [datasetReduction_TargetTime, 'Epoch']

plasmaSheet_targetTimes = [dt.datetime(2022,11, 20, 17, 24, 21, 000000), dt.datetime(2022,11,20,17,24,48,000000)]
plasmaSheet_TargetTimes_data = [  [dt.datetime(2022, 11, 20, 17, 24, 24, 000000), dt.datetime(2022,11,20,17,24,30,500000)],
                             [dt.datetime(2022, 11, 20, 17, 24, 39, 500000), dt.datetime(2022,11,20,17,24,45,000000)]
                             ]

invertedV_targetTimes = [dt.datetime(2022,11, 20, 17, 24, 56, 000000), dt.datetime(2022,11,20,17,25,3,000000)]
invertedV_TargetTimes_data = [[dt.datetime(2022,11, 20, 17, 25, 1, 396000), dt.datetime(2022,11,20,17,25,1,596000)],
                         [dt.datetime(2022,11, 20, 17, 25, 1, 462000), dt.datetime(2022,11,20,17,25,1,512000)]
                         ]

# --- Plot toggles - General ---
figure_width = 10 # in inches
figure_height =18 # in inches
Title_FontSize = 20
Label_FontSize = 20
Label_Padding = 8
Tick_FontSize = 12
Tick_Length = 1
Tick_Width = 1
Tick_FontSize_minor = 10
Tick_Length_minor = 1
Tick_Width_minor = 1
Plot_LineWidth = 0.5
plot_MarkerSize = 14
legend_fontSize = 15
dpi = 800

# --- Cbar ---
mycmap = apl_rainbow_black0_cmap()
cbarMin, cbarMax = 1E-18, 1E-14
cbarTickLabelSize = 14
# cmap.set_bad(color=(1, 1, 1))

# --- Distribution Fit toggles---
fit_figure_width = 10 # in inches
fit_figure_height =10 # in inches
Plot_DistributionFit = True # outputs a sets of slices in pitch angle vs distrubution function
DistributionFit_wPtchSlice = 2

# --- Dist Overview toggles ---
show_OverviewPlot = False
OverviewPlot_ptchRange = [2] # determines which pitch to include in the parallel distribution functions

# --- Comparison Plot toggles ---
show_ComparisonPlot = False
ComparisonPlot_distLimits = [1E-19, 1E-13]
ComparisonPlot_xAxis_energyLimits = [10, 1E4]
ComparisonPlot_ptchRange = [1]
ComparisonPlot_countNoiseLevel = 1

# --- Backmapped toggles ---
show_BackmappedPlot = False
BackmappedPlot_altsMapped = [400, 2000, 5000] # altitude [in kilometers] to determine the distribution functions at higher altitudes
BackmappedPlot_ptchRange = [i for i in range(11)] # for the
BackmappedPlot_X_Velocity_limits, BackmappedPlot_Y_Velocity_limit = [-0.5, 1.8], [-0.5, 2]


##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
prgMsg('Loading Data')

rocketAttrs, b, c = ACES_mission_dicts()

# EEPAA Distribution Data
inputEEPAA_dist_high = glob('C:\Data\ACESII\L3\DistFunc\high\*ACESII_36359_distFunc_eepaa_fullCal*')[0]
data_dict_dist = loadDictFromFile(inputFilePath=inputEEPAA_dist_high, targetVar=targetVar, wKeys_Reduce=['Distribution_Function', 'Epoch'])
inputEEPAA_counts_high = glob('C:\Data\ACESII\L1\high\*eepaa_fullCal*')[0]
data_dict_counts = loadDictFromFile(inputFilePath=inputEEPAA_counts_high, targetVar=targetVar, wKeys_Reduce=['eepaa', 'Epoch'])

# Define the data
counts = deepcopy(data_dict_counts['eepaa'][0])
distFunc = deepcopy(data_dict_dist['Distribution_Function'][0])
Epoch = deepcopy(data_dict_dist['Epoch'][0])
Energy = deepcopy(data_dict_dist['Energy'][0])
Pitch = deepcopy(data_dict_dist['Pitch_Angle'][0])
Done(start_time)

##############################
# --- --- --- --- --- --- ----
# --- CALCULATED VARIABLES ---
# --- --- --- --- --- --- ----
##############################
prgMsg('Preparing Data')

###### 1-Count THRESHOLD LINE ######
distFunc_NoiseCount = np.zeros(shape=(len(Energy)))
geo_factor = rocketAttrs.geometric_factor[0]
count_interval = 0.8992E-3

for engy in range(len(Energy)):
    deltaT = (count_interval) - (ComparisonPlot_countNoiseLevel * rocketAttrs.deadtime[0])
    fluxVal = (ComparisonPlot_countNoiseLevel) / (geo_factor[0] * deltaT)
    distFunc_NoiseCount[engy] = ((cm_to_m * m_e / (q0 * Energy[engy])) ** 2) * (fluxVal / 2)

# Time since launch
LaunchDateTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 000000))
Epoch_timeSince = EpochTo_T0_Rocket(InputEpoch=Epoch, T0=LaunchDateTime)
plasmaSheet_timeSince = EpochTo_T0_Rocket(InputEpoch=plasmaSheet_targetTimes, T0=LaunchDateTime)
invertedV_timeSince = EpochTo_T0_Rocket(InputEpoch=invertedV_targetTimes, T0=LaunchDateTime)

if show_OverviewPlot:

    ###### averaged distribution function over pitch angle ######
    omniDist_high = np.zeros(shape=(len(distFunc), len(Energy)))
    for tme in range(len(Epoch)):
        for engy in range(len(Energy)):
            sumVal = 0

            for ptch in OverviewPlot_ptchRange:
                val = distFunc[tme, ptch, engy]
                if val > 0:
                    sumVal += val

            # Average the Omni-flux by the number of bins. ONLY include bins 10deg - 170 since they have full coverage
            omniDist_high[tme][engy] = sumVal / len(range(2, 18 + 1))

if show_ComparisonPlot:


    ###### PARALLEL DISTRIBUTION FUNCTION VS ENERGY ######
    paraDist = [[], []]
    paraEngy = [[], []]
    regionTimes = [plasmaSheet_targetTimes, invertedV_targetTimes]

    for regionIdx in range(len(paraDist)):
        rTime = regionTimes[regionIdx]
        low_idx, high_idx = np.abs(Epoch - rTime[0]).argmin(), np.abs(Epoch - rTime[1]).argmin()
        tempDist = []
        tempEngy = []
        tempPtch = []

        for ptch in ComparisonPlot_ptchRange:
            data = distFunc[low_idx:high_idx, ptch, :].T
            data[np.where(data == rocketAttrs.epoch_fillVal)] = 0

            average = np.array([np.average(arr[np.nonzero(arr)])    for arr in data]) # average over time for each energy bin
            average = np.array([val if np.isnan(val) == False else 0 for val in average ]) # remove the points where there was no distribution function for a specific enregy

            # remove all the zero values from average and from energy
            zeroFinder = np.where(average!=0)


            tempDist.append(average[zeroFinder])
            tempEngy.append(Energy[zeroFinder])
            tempPtch.append(ptch)

        if len(tempPtch) > 1:
            print('GOTTA WRITE SOMETHING')

        else:
            paraDist[regionIdx] = tempDist[0]
            paraEngy[regionIdx] = tempEngy[0]

    Done(start_time)

if show_BackmappedPlot:
    ###### BACK-MAPPING COUNTS TO DISTRIBUTION FUNCTION ######
    prgMsg('Mapping Data')
    # Time since launch
    LaunchDateTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 000000))
    Epoch_timeSince = EpochTo_T0_Rocket(InputEpoch=Epoch, T0=LaunchDateTime)
    plasmaSheet_timeSince = EpochTo_T0_Rocket(InputEpoch=plasmaSheet_targetTimes, T0=LaunchDateTime)
    invertedV_timeSince = EpochTo_T0_Rocket(InputEpoch=invertedV_targetTimes, T0=LaunchDateTime)
    regionTimes = [plasmaSheet_targetTimes,invertedV_targetTimes]
    mappedCounts_output = [[], []]

    # loop over altitudes
    for alt in BackmappedPlot_altsMapped:

        # loop over plasma sheet or inverted-V altitudes
        for regionIdx in range(len(regionTimes)):

            rTime = regionTimes[regionIdx]
            low_idx, high_idx = np.abs(Epoch - rTime[0]).argmin(), np.abs(Epoch - rTime[1]).argmin()
            countsData = counts[low_idx:high_idx]
            mappedCounts = [[[[] for engy in range(len(Energy))] for ptch in range(len(Pitch))] for tme in range(len(counts[low_idx:high_idx]))]

            # loop over all values in the data and re-bin them to a new pitch angle position based on altitude
            sizes = [len(countsData), len(countsData[0]), len(countsData[0][0])]
            ranges = [range(sizes[0]), range(sizes[1]), range(sizes[2])]

            for tme, ptch,engy in itertools.product(*ranges):

                currPitch = Pitch[ptch]
                mappedPitch = np.degrees(np.arcsin(  np.sin(np.radians(currPitch))*np.power((Re + 400)/(Re + alt), 1.5)))

                # determine where new mapped pitch should go
                newPtchIdx = np.abs(Pitch - mappedPitch).argmin()
                mappedCounts[tme][newPtchIdx][engy].append(countsData[tme][ptch][engy])

            # reduce the lists by summing the total counts at each point
            for tme, ptch, engy in itertools.product(*ranges):

                val = np.array(mappedCounts[tme][ptch][engy])

                # remove fillvals from data
                val = val[np.where(val > 0)]

                if len(val) == 0:
                    mappedCounts[tme][ptch][engy] = 0
                else:
                    mappedCounts[tme][ptch][engy] = sum(val)

            # Calculate the Distribution Function for each point
            geo_factor = rocketAttrs.geometric_factor[0]
            count_interval = 0.8992E-3

            for tme, ptch, engy in itertools.product(*ranges):
                deltaT = (count_interval) - (mappedCounts[tme][ptch][engy] * rocketAttrs.deadtime[0])
                fluxVal = (mappedCounts[tme][ptch][engy]) / (geo_factor[0] * deltaT)
                mappedCounts[tme][ptch][engy] = ((cm_to_m * m_e / (q0 * Energy[engy])) ** 2) * (fluxVal / 2)

            mappedCounts_output[regionIdx].append(np.array(mappedCounts))

    mapped_paraDist = [[[], []] for alt in BackmappedPlot_altsMapped]
    mapped_paraEngy = [[[], []] for alt in BackmappedPlot_altsMapped]
    mapped_paraPtch = [[[], []] for alt in BackmappedPlot_altsMapped]

    for altidx, alt in enumerate(BackmappedPlot_altsMapped):

        for regionIdx in range(len(mapped_paraDist[0])):
            rTime = regionTimes[regionIdx]
            low_idx, high_idx = np.abs(Epoch - rTime[0]).argmin(), np.abs(Epoch - rTime[1]).argmin()
            tempDist = []
            tempEngy = []
            tempPtch = []

            for ptch in BackmappedPlot_ptchRange:

                data = mappedCounts_output[regionIdx][altidx][:, ptch, :].T

                average = np.array([np.average(arr[np.nonzero(arr)]) for arr in data]) # average over time for each energy bin
                average = np.array([val if np.isnan(val) == False else 0 for val in average ]) # remove the points where there was no distribution function for a specific enregy

                # remove all the zero values from average and from energy
                # zeroFinder = np.where(average!=0)
                # tempDist.append(average[zeroFinder])
                # tempEngy.append(Energy[zeroFinder])
                tempDist.append(average)
                tempEngy.append(Energy)
                tempPtch.append(ptch)

            mapped_paraDist[altidx][regionIdx] = tempDist
            mapped_paraEngy[altidx][regionIdx] = tempEngy
            mapped_paraPtch[altidx][regionIdx] = tempPtch
    Done(start_time)


    ###### CALCUALTE V_PERP/V_PARALLEL ######
    prgMsg('Calculating Vperp and Vparallel')
    countsTemp = deepcopy(counts[0][0:len(BackmappedPlot_ptchRange)])

    for regionIdx, tmes in enumerate(regionTimes):

        Vperp = deepcopy(countsTemp)
        Vpara = deepcopy(countsTemp)

        for ptch in BackmappedPlot_ptchRange:
            for engy in range(len(Energy)):
                Emag = np.sqrt(2 * q0 * Energy[engy] / (m_e))
                Vperp[ptch][engy] = np.sin(np.radians(Pitch[ptch])) * Emag
                Vpara[ptch][engy] = np.cos(np.radians(Pitch[ptch])) * Emag

    Vpara, Vperp = np.array(Vpara) / (10000 * 1000), np.array(Vperp) / (10000 * 1000)
    Done(start_time)

###############################
# --- --- --- --- --- --- --- -
# --- FIT THE DISTRIBUTIONS ---
# --- --- --- --- --- --- --- -
###############################
from numpy import pi, exp, sqrt
from math import gamma

if Plot_DistributionFit:

    ##################################
    # --- DEFINE THE FIT FUNCTIONS ---
    ##################################
    def MaxwellianDist(x, n, T):  # Fits the NATURAL LOG of a Maxwellian for a Uniform/Thermalized Plasma
        return 4 * pi * n * ((m_e / (2 * pi * q0 * T)) ** 1.5) * (2 * x / m_e) * exp(-x / (T))

    def KappaDist(x, kappa, n, T):
        w = (2 * q0 * T / m_e) * ((kappa - 3 / 2) / kappa)
        term1 = n / (2 *(pi**1.5)* (w**1.5) )
        term2 = gamma(kappa+1) / ((kappa**1.5)*gamma(kappa + 0.5))
        term3 = (1 + (x**2) /(kappa * w))**(-kappa - 1)
        func = term1 * term2*term3

        # term1 = n / (power(pi, 3/2)*power(w,3))
        # term2 = gamma(kappa+1)/(power(kappa,3/2)*gamma(kappa-1/2))
        # term3 = (1 + power(x,2)/(kappa*power(w,2)) )**(-kappa-1)

        # term1 = n * gamma(kappa + 1)
        # term2 = (1 / (((sqrt(pi) * (sqrt(q0 * T * (2 * kappa - 3) / (kappa * m_e)))) ** 3) * ((kappa ** 3 / 2) * gamma(kappa - 0.5))))
        # term3 = ((1 + 2 * q0 * x / (m_e * kappa * ((sqrt(q0 * T * (2 * kappa - 3) / (kappa * m_e))) ** 2))) ** (-kappa - 1))
        return func

    def collectDistData(datasetTimes, datasetPitchRange):

        ###### PARALLEL DISTRIBUTION FUNCTION VS ENERGY ######
        distributionVals= []
        energyVals = []
        pitchVals = []
        low_idx, high_idx = np.abs(Epoch - datasetTimes[0]).argmin(), np.abs(Epoch - datasetTimes[1]).argmin()

        for ptch in datasetPitchRange:
            data = distFunc[low_idx:high_idx, ptch, :].T
            data[np.where(data == rocketAttrs.epoch_fillVal)] = 0
            average = np.array([np.average(arr[np.nonzero(arr)]) for arr in data])  # average over time for each energy bin
            average = np.array([val if np.isnan(val) == False else 0 for val in average])  # remove the points where there was no distribution function for a specific enregy

            # remove all the zero values from average and from energy
            zeroFinder = np.where(average != 0)
            distributionVals.append(average[zeroFinder])
            energyVals.append(Energy[zeroFinder])
            pitchVals.append(ptch)

        return distributionVals, energyVals, pitchVals

    ##############################
    # --- COLLECT THE FIT DATA ---
    ##############################
    fitPtches = [2, 3, 4] # pitch angles 10deg, 20deg, 30deg
    InvertedVData = [collectDistData(invertedV_TargetTimes_data[i], fitPtches) for i in range(len(invertedV_TargetTimes_data))]
    PlasmaSheetData = [collectDistData(plasmaSheet_TargetTimes_data[i], fitPtches) for i in range(len(plasmaSheet_TargetTimes_data))]

    # guess = [1.509 1.3E6, 500]
    # NOTE: 1 cm^3 = 1E6 m^3. Magnetospheric particles from plasma sheet are around 1 cm^-3 at 500eV
    # guess = [1.51, 1E6, 200]  # observed plasma at dispersive region is 0.5E5 cm^-3 BUT this doesn't make sense to use as the kappa fit since the kappa fit comes from MUCH less dense populations above
    # boundVals = [[1.5000000001, 2], [1E4, 1E14], [0.001, 1E4]]  # kappa, n, T
    # bounds = tuple([[boundVals[i][0] for i in range(len(boundVals))], [boundVals[i][1] for i in range(len(boundVals))]])
    # params, cov = curve_fit(KappaDist, paraEngy[0], paraDist[0], p0=guess, maxfev=int(1E9), bounds=bounds)
    # xData_fit, yData_fit = np.array(paraEngy[0]), np.array([KappaDist(val, *params) for val in paraEngy[0]])

    # --- PLOT THE DISTRIBUTION FIT---
    dataTimes = [plasmaSheet_TargetTimes_data, invertedV_TargetTimes_data]
    dataDists = [PlasmaSheetData,InvertedVData]
    wRegionIdx = 0 # determines which region we're interested in
    wTimeIdx = 0
    wPtchIdx = 0
    rLabels = ['Plasma Sheet', 'Inverted-V']

    for wRegionIdx in [0, 1]:
        for wTimeIdx in range(len(dataTimes[wRegionIdx])):
            for wPtchIdx in range(len(fitPtches)):
                fitTimes = dataTimes[wRegionIdx][wTimeIdx] # [timestart,timeend]
                fitDist = dataDists[wRegionIdx][wTimeIdx][0][wPtchIdx]
                fitEngy = dataDists[wRegionIdx][wTimeIdx][1][wPtchIdx]
                fitPtch = dataDists[wRegionIdx][wTimeIdx][2][wPtchIdx]

                fig, ax = plt.subplots(2)
                fig.set_size_inches(fit_figure_width, fit_figure_height)

                # general overplot
                low_idx, high_idx = np.abs(Epoch - fitTimes[0]).argmin(), np.abs(Epoch - fitTimes[1]).argmin()
                fitEpoch = Epoch[low_idx:high_idx+1]
                generaldata = distFunc[low_idx:high_idx+1, fitPtch, :].T
                ax[0].set_title(rLabels[wRegionIdx]+f'\n {fitTimes[0].strftime("%H:%M:%S")} to {fitTimes[1].strftime("%H:%M:%S")}' +f'\nPitch = {Pitch[fitPtch]}', fontsize=Title_FontSize)
                ax[0].pcolormesh(fitEpoch,Energy,generaldata, cmap=mycmap, vmin=cbarMin, vmax=cbarMax, norm='log')
                ax[0].set_yscale('log')
                ax[0].set_ylabel('Energy [eV]', fontsize=Label_FontSize)

                ax[1].plot(fitEngy, fitDist, color='black', marker='.', ms=plot_MarkerSize-7)
                ax[1].plot(Energy, distFunc_NoiseCount, label=f'{ComparisonPlot_countNoiseLevel}-count level', color='black')
                ax[1].set_yscale('log')
                ax[1].set_ylim(ComparisonPlot_distLimits[0], ComparisonPlot_distLimits[1])
                ax[1].set_xlim(ComparisonPlot_xAxis_energyLimits[0], ComparisonPlot_xAxis_energyLimits[1])
                ax[1].set_xscale('symlog')

                for i in range(2):
                    ax[i].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
                    ax[i].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize_minor, length=Tick_Length_minor, width=Tick_Width_minor)
                plt.tight_layout()
                plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\fitData\Plot8_distribution_fit_base{wRegionIdx}{wTimeIdx}{wPtchIdx}.png',dpi=dpi)

else:
    ############################
    # --- --- --- --- --- --- --
    # --- START THE PLOTTING ---
    # --- --- --- --- --- --- --
    ############################
    prgMsg('Plotting Figure')

    # DETERMINE THE PLOT DIMENSIONS
    fig = plt.figure()
    fig.set_size_inches(figure_width, figure_height)
    gs0 = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[2/10, 2/10, 6/10], hspace=0.12)

    # --- DISTRIBUTION DATA OVERVIEW---
    if show_OverviewPlot:
        gs00 = gs0[0].subgridspec(1,1)
        axDistOverview = fig.add_subplot(gs00[0])
        axDistOverview.pcolormesh(Epoch_timeSince, Energy, omniDist_high.T, cmap=mycmap, vmin=cbarMin, vmax=cbarMax, norm='log')
        axDistOverview.set_ylabel('Energy [eV]', fontsize=Label_FontSize)
        axDistOverview.tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
        axDistOverview.tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize_minor, length=Tick_Length_minor, width=Tick_Width_minor)
        axDistOverview.set_yscale('log')
        axDistOverview.axvspan(invertedV_timeSince[0], invertedV_timeSince[1], color='gray', alpha=0.25, zorder=1)
        axDistOverview.axvspan(plasmaSheet_timeSince[0], plasmaSheet_timeSince[1], color='gray', alpha=0.25, zorder=1)

        # --- cbar ---
        # cax = fig.add_axes([0.91, 0.288, 0.02, 0.592])
        # cbar = plt.colorbar(cmap, cax= cax)
        # cbar.set_label('Omni-Dir. diff E. Flux \n' + '[cm$^{-2}$str$^{-1}$eV/eV]', rotation=-90, labelpad=20, fontsize=General_LabelFontSize)
        # cbar.ax.minorticks_on()
        # cbar.ax.tick_params(labelsize=cbarTickLabelSize + 5)

    if show_ComparisonPlot:
        # --- COMPARISON DISTRIBUTIONS ---
        gs01 = gs0[1].subgridspec(1, 1)
        ax_Compare = fig.add_subplot(gs01[0, 0])
        ax_Compare.scatter(paraEngy[0], paraDist[0], label='Plasma Sheet Ambient', color='blue')
        ax_Compare.scatter(paraEngy[1], paraDist[1], label='Inverted-V', color='red')
        ax_Compare.plot(Energy, distFunc_NoiseCount, label=f'{ComparisonPlot_countNoiseLevel}-count level',color='black')
        ax_Compare.set_yscale('log')
        ax_Compare.set_ylim(ComparisonPlot_distLimits[0], ComparisonPlot_distLimits[1])
        ax_Compare.set_xlim(ComparisonPlot_xAxis_energyLimits[0], ComparisonPlot_xAxis_energyLimits[1])
        ax_Compare.set_xscale('symlog')
        ax_Compare.grid()
        ax_Compare.set_ylabel(r'Parallel Distribution Function',fontsize=Label_FontSize)
        ax_Compare.legend(prop={'size': legend_fontSize})

    if show_BackmappedPlot:
        # --- BACK-MAPPED DISTRIBUTIONS ---
        gs02 = gs0[2].subgridspec(3, 2)
        subAxes = [[fig.add_subplot(gs02[j, i]) for j in range(3)] for i in range(2)]
        backmapped_labels = ['P.S.', 'Inverted-V']
        for rIdx in range(2):
            for z, alt in enumerate(BackmappedPlot_altsMapped):
                subAxes[rIdx][z].pcolormesh(Vperp, Vpara, mapped_paraDist[z][rIdx], cmap=mycmap, shading='nearest', norm='log', vmin=cbarMin, vmax=cbarMax)
                subAxes[rIdx][z].scatter([0], [0], alpha=0, label=f'{backmapped_labels[rIdx]} Alt {alt} [km]')
                subAxes[rIdx][z].set_xlim(BackmappedPlot_X_Velocity_limits[0], BackmappedPlot_X_Velocity_limits[1])
                subAxes[rIdx][z].set_ylim(BackmappedPlot_Y_Velocity_limit[0], BackmappedPlot_Y_Velocity_limit[1])
                subAxes[rIdx][z].invert_yaxis()
                subAxes[rIdx][z].legend()


    plt.tight_layout()
    # fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)  # remove the space between plots
    plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot8_DistributionFunc_base.png', dpi=dpi)
    Done(start_time)

