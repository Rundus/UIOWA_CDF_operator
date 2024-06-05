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
import numpy as np

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
datasetReduction_TargetTime = [dt.datetime(2022,11, 20, 17, 24, 50, 000000), dt.datetime(2022,11,20,17,25,15,000000)]
targetVar = [datasetReduction_TargetTime, 'Epoch']

plasmaSheet_targetTimes = [dt.datetime(2022,11, 20, 17, 24, 21, 000000), dt.datetime(2022,11,20,17,24,48,000000)]
plasmaSheet_TargetTimes_data = [
                            [dt.datetime(2022, 11, 20, 17, 26, 22, 500000), dt.datetime(2022,11,20,17,26,28,500000)],
                            [dt.datetime(2022, 11, 20, 17, 27, 0, 000000), dt.datetime(2022,11,20,17,27,10,000000)]
                             ]

invertedV_targetTimes = [dt.datetime(2022,11, 20, 17, 24, 56, 000000), dt.datetime(2022,11,20,17,25,3,000000)]
# invertedV_TargetTimes_data = [[dt.datetime(2022,11, 20, 17, 25, 1, 396000), dt.datetime(2022,11,20,17,25,1,712000)]]
invertedV_TargetTimes_data = [[dt.datetime(2022,11, 20, 17, 25, 1, 300000), dt.datetime(2022,11,20,17,25,1,800000)]]

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
dpi = 200

# --- Cbar ---
mycmap = apl_rainbow_black0_cmap()
cbarMin, cbarMax = 1E-18, 1E-14
cbarTickLabelSize = 14
# cmap.set_bad(color=(1, 1, 1))


##################################
Fit_FITDATA = True

# --- Distribution Fit (AVERAGE) toggles---
Fit_AverageDistributionFit = False # outputs a sets of slices in pitch angle vs distrubution function
fit_figure_width = 10 # in inches
fit_figure_height =10 # in inches
DistributionFit_wPtchSlice = 2

# --- Distribution Fit - Density, Temperature and Potential ---
Fit_DistributionDensityTempPotental = True
wPitchToFit = 2


# --- Dist Overview toggles ---
show_OverviewPlot = False
OverviewPlot_ptchRange = [2] # determines which pitch to include in the parallel distribution functions

# --- Comparison Plot toggles ---
show_ComparisonPlot = False
ComparisonPlot_distLimits = [1E-21, 1E-9]
ComparisonPlot_xAxis_energyLimits = [1E-2, 1E4]
ComparisonPlot_ptchRange = [1]
ComparisonPlot_countNoiseLevel = 4

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
inputEEPAA_diffFlux_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0]
data_dict_diffFlux = loadDictFromFile(inputEEPAA_diffFlux_high, targetVar=targetVar, wKeys_Reduce=['Differential_Energy_Flux','Differential_Number_Flux', 'Epoch'])

# Define the data
counts = deepcopy(data_dict_counts['eepaa'][0])
distFunc = deepcopy(data_dict_dist['Distribution_Function'][0])
diffEFlux = deepcopy(data_dict_diffFlux['Differential_Energy_Flux'][0])
diffNFlux = deepcopy(data_dict_diffFlux['Differential_Number_Flux'][0])
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
    regionTimes = [plasmaSheet_targetTimes, invertedV_targetTimes]
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


if Fit_FITDATA:

    if Fit_AverageDistributionFit:

        from numpy import pi, exp, sqrt
        from math import gamma

        ##################################
        # --- DEFINE THE FIT FUNCTIONS ---
        ##################################
        def MaxwellianDist(x, n, T):  # Fits the NATURAL LOG of a Maxwellian for a Uniform/Thermalized Plasma
            return 4 * pi * n * ((m_e / (2 * pi * q0 * T)) ** 1.5) * (2 * x / m_e) * exp(-x / (T))

        def KappaDist(x, kappa, n, T):
            w = (2 * q0 * T / m_e) * ((kappa - 3 / 2) / kappa)
            term1 = n / (2 *(pi**1.5)* (w**1.5) )
            term2 = gamma(kappa+1) / ((kappa**1.5)*gamma(kappa + 0.5))
            term3 = (1 + (2*x*q0) /(m_e*kappa * w))**(-kappa - 1)
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
        fitPtches = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                     18]  # pitch angles 10deg, 20deg, 30deg
        InvertedVData = [collectDistData(invertedV_TargetTimes_data[i], fitPtches) for i in
                         range(len(invertedV_TargetTimes_data))]
        PlasmaSheetData = [collectDistData(plasmaSheet_TargetTimes_data[i], fitPtches) for i in
                           range(len(plasmaSheet_TargetTimes_data))]

        # --- PLOT THE DISTRIBUTION FIT---
        dataTimes = [plasmaSheet_TargetTimes_data, invertedV_TargetTimes_data]
        dataDists = [PlasmaSheetData, InvertedVData]
        wRegionIdx = 0  # determines which region we're interested in
        wTimeIdx = 0
        wPtchIdx = 0
        rLabels = ['Plasma Sheet', 'Inverted-V']

        from matplotlib.widgets import Button, Slider

        init_kappa = 1.572
        init_Temp = 485
        init_n = 9.56E6

        for wRegionIdx in [0, 1]:
            for wTimeIdx in range(len(dataTimes[wRegionIdx])):
                for wPtchIdx in range(len(fitPtches)):
                    fitTimes = dataTimes[wRegionIdx][wTimeIdx]  # [timestart,timeend]
                    fitDist = dataDists[wRegionIdx][wTimeIdx][0][wPtchIdx]
                    fitEngy = dataDists[wRegionIdx][wTimeIdx][1][wPtchIdx]
                    fitPtch = dataDists[wRegionIdx][wTimeIdx][2][wPtchIdx]

                    # NOTE: 1 cm^3 = 1E6 m^3. Magnetospheric particles from plasma sheet are around 1 cm^-3 at 500eV
                    # guess = [1.55, 20000000, 100]  # observed plasma at dispersive region is 0.5E5 cm^-3 BUT this doesn't make sense to use as the kappa fit since the kappa fit comes from MUCH less dense populations above
                    # boundVals = [[1.500001, 2], [1E4,1E6], [10, 100]]  # kappa, n, T
                    # bounds = tuple([[boundVals[i][0] for i in range(len(boundVals))], [boundVals[i][1] for i in range(len(boundVals))]])
                    # # params, cov = curve_fit(KappaDist, fitEngy, fitDist, maxfev=int(1E9), bounds=bounds, p0=guess)
                    # params, cov = curve_fit(KappaDist, fitEngy, fitDist, maxfev=int(1E9), bounds=bounds)
                    # # params, cov = curve_fit(KappaDist, fitEngy, fitDist, maxfev=int(1E9), p0=guess)
                    # xData_fit, yData_fit = np.array(fitEngy), np.array([KappaDist(val, *params) for val in fitEngy])

                    # Calculate the ChiSquare value
                    nu = 3 - 1
                    # ChiSquared = (1/nu) * sum([ (yData_fit[i] - fitDist[i])**2 / (fitDist[i]**2) for i in range(len(xData_fit))])

                    fig, ax = plt.subplots(2)
                    fig.set_size_inches(fit_figure_width, fit_figure_height)

                    # general overplot
                    low_idx, high_idx = np.abs(Epoch - fitTimes[0]).argmin(), np.abs(Epoch - fitTimes[1]).argmin()
                    fitEpoch = Epoch[low_idx:high_idx + 1]
                    generaldata = distFunc[low_idx:high_idx + 1, fitPtch, :].T
                    ax[0].set_title(rLabels[
                                        wRegionIdx] + f'\n {fitTimes[0].strftime("%H:%M:%S")} to {fitTimes[1].strftime("%H:%M:%S")}' + f'\nPitch = {Pitch[fitPtch]}',
                                    fontsize=Title_FontSize)
                    ax[0].pcolormesh(fitEpoch, Energy, generaldata, cmap=mycmap, vmin=cbarMin, vmax=cbarMax,
                                     norm='log')
                    ax[0].set_yscale('log')
                    ax[0].set_ylabel('Energy [eV]', fontsize=Label_FontSize)

                    ax[1].plot(Energy, distFunc_NoiseCount, label=f'{ComparisonPlot_countNoiseLevel}-count level',
                               color='red')
                    ax[1].plot(fitEngy, fitDist, color='black', marker='.', ms=plot_MarkerSize - 7)
                    fittedData, = ax[1].plot(fitEngy, KappaDist(fitEngy, init_kappa, init_n, init_Temp), lw=2)
                    # ax[1].plot(xData_fit, yData_fit, color='green', alpha=0.8, label=rf'$\kappa = {params[0]}$' + f'\nn={params[1]}' + f'\nT={params[2]}' + '\n $\chi^{2}$ =' +f'{ChiSquared}')
                    ax[1].set_yscale('log')
                    ax[1].set_ylim(ComparisonPlot_distLimits[0], ComparisonPlot_distLimits[1])
                    ax[1].set_xlim(ComparisonPlot_xAxis_energyLimits[0], ComparisonPlot_xAxis_energyLimits[1])
                    ax[1].set_xscale('symlog')
                    ax[1].legend(fontsize=legend_fontSize)

                    # --- Create the Sliders ---

                    # Kappa
                    axKappa = fig.add_axes([0.75, 0.25, 0.0225, 0.63])
                    kappa_slider = Slider(
                        ax=axKappa,
                        label="Kappa",
                        valmin=1.5000001,
                        valmax=3,
                        valinit=init_kappa,
                        orientation="vertical"
                    )

                    # Density
                    axDensity = fig.add_axes([0.85, 0.25, 0.0225, 0.63])
                    density_slider = Slider(
                        ax=axDensity,
                        label="density",
                        valmin=1E3,
                        valmax=1E8,
                        valinit=init_n,
                        orientation="vertical"
                    )

                    # temperature
                    axTemp = fig.add_axes([0.95, 0.25, 0.0225, 0.63])
                    Temp_slider = Slider(
                        ax=axTemp,
                        label="Temp",
                        valmin=1,
                        valmax=1000,
                        valinit=init_Temp,
                        orientation="vertical"
                    )

                    def update(val):
                        fittedData.set_ydata(
                            KappaDist(fitEngy, kappa_slider.val, density_slider.val, Temp_slider.val))
                        fig.canvas.draw_idle()

                    kappa_slider.on_changed(update)
                    density_slider.on_changed(update)
                    Temp_slider.on_changed(update)

                    for i in range(2):
                        ax[i].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize,
                                          length=Tick_Length, width=Tick_Width)
                        ax[i].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize_minor,
                                          length=Tick_Length_minor, width=Tick_Width_minor)
                    plt.tight_layout()
                    # plt.show()
                    # plt.close()
                    plt.savefig(
                        rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\fitData\Plot8_distribution_fit_base{wRegionIdx}{wTimeIdx}{wPtchIdx}.png',
                        dpi=dpi)

    if Fit_DistributionDensityTempPotental:


        def diffNFlux_for_mappedMaxwellian(x,n,T,beta,V,alpha):
            Vpara_sqrd = (2*x*np.power(np.cos(np.radians(alpha)),2)/m_e  ) - 2*V/m_e + (1 - 1/beta)*(2*x/m_e)*(np.power(np.sin(np.radians(alpha)),2))
            Vperp_sqrd = ((2*x)/(beta*m_e))*np.power(np.sin(np.radians(alpha)), 2)

            return (2*x)*((q0/m_e)**2) * (1E2*n) * np.power(m_e/(2*np.pi*q0*T), 3/2) * np.exp((-m_e/(2*T))*(Vpara_sqrd + Vperp_sqrd))

        # define my function at the specific pitch angle of interest
        from functools import partial
        fitFuncAtPitch = partial(diffNFlux_for_mappedMaxwellian, alpha= Pitch[wPitchToFit])

        ##############################
        # --- COLLECT THE FIT DATA ---
        ##############################

        for timeset in invertedV_TargetTimes_data:
            # collect the data
            low_idx, high_idx = np.abs(Epoch - timeset[0]).argmin(), np.abs(Epoch - timeset[1]).argmin()
            EpochFitData = Epoch[low_idx:high_idx + 1]
            fitData = diffNFlux[low_idx:high_idx+1, wPitchToFit, :]

            # for each slice in time, loop over the data and identify the peak differentialNumberFlux (This corresponds to the
            # peak energy of the inverted-V since the location of the maximum number flux tells you what energy the low-energy BULk got accelerated to)
            # Note: The peak in the number flux is very likely the maximum value AFTER 100 eV, just find this point

            for tmeIdx in range(len(EpochFitData)):

                # Determine the peak point based on a treshold limit
                threshEngy = 200
                EngyIdx = np.abs(Energy - threshEngy).argmin()
                peakDiffNVal = fitData[tmeIdx][:EngyIdx].max()
                peakDiffNVal_index = np.argmax(fitData[tmeIdx][:EngyIdx])

                # get the subset of data to fit to and fit it. Only include data with non-zero points
                xData_fit = np.array(Energy[:peakDiffNVal_index+1])
                yData_fit = np.array(fitData[tmeIdx][:peakDiffNVal_index+1])
                nonZeroIndicies = np.where(yData_fit!=0)[0]
                xData_fit = xData_fit[nonZeroIndicies]
                yData_fit = yData_fit[nonZeroIndicies]

                deviation = 0.18
                guess = [1.55, 20000000, 100]  # observed plasma at dispersive region is 0.5E5 cm^-3 BUT this doesn't make sense to use as the kappa fit since the kappa fit comes from MUCH less dense populations above
                boundVals = [[0.0001, 30000], # n [cm^-3]
                             [10,500], # T [eV]
                             [1, 10], # beta [unitless]
                             [(1-deviation)*Energy[peakDiffNVal_index], (1+deviation)*Energy[peakDiffNVal_index]]]  # V [eV]
                bounds = tuple([[boundVals[i][0] for i in range(len(boundVals))], [boundVals[i][1] for i in range(len(boundVals))]])
                params, cov = curve_fit(fitFuncAtPitch,xData_fit,yData_fit,maxfev=int(1E9), bounds=bounds)

                print(cov)

                pairs = [f'({x},{y})' for x,y in zip(xData_fit,yData_fit)]
                print(Energy[peakDiffNVal_index])
                for pairVal in pairs:
                    print(pairVal)

                fittedX = np.linspace(xData_fit.min(), xData_fit.max(), 100)
                fittedY = fitFuncAtPitch(fittedX,*params)

                # --- Calculate ChiSquare ---

                for i in range(len(xData_fit)):
                    print(fitFuncAtPitch(xData_fit[i],*params),yData_fit[i],xData_fit[i],(fitFuncAtPitch(xData_fit[i],*params) - yData_fit[i])**2)

                ChiSquare = (1/2)*sum([(fitFuncAtPitch(xData_fit[i],*params) - yData_fit[i])**2 / (yData_fit[i]**2) for i in range(len(xData_fit))])


                # Plot the result to see
                fig, ax = plt.subplots(2)
                fig.suptitle(f'Pitch Angle = {Pitch[wPitchToFit]} \n {EpochFitData[tmeIdx]} UTC')
                ax[0].pcolormesh(EpochFitData,Energy,fitData.T, vmin=1E4, vmax=1E7,cmap='turbo', norm='log')
                ax[0].set_yscale('log')
                ax[0].set_ylabel('Energy [eV]')
                ax[0].set_xlabel('Time')
                ax[0].axvline(EpochFitData[tmeIdx],color='black', linestyle='--')
                ax[0].set_ylim(28,1000)
                ax[1].plot(Energy, fitData[tmeIdx][:],'-o')
                ax[1].plot(fittedX, fittedY, color='purple', label=f'n = {params[0]}\n T = {params[1]} \n' + rf'$\beta$={params[2]}' + f'\n V = {params[3]}\n' + r'$\chi^{2}$='+f'{ChiSquare}')
                ax[1].legend()
                ax[1].axvline(Energy[peakDiffNVal_index],color='red')
                ax[1].set_yscale('log')
                ax[1].set_xscale('log')
                ax[1].set_xlabel('Energy [eV]')
                ax[1].set_ylabel('diffNFlux')
                ax[1].set_xlim(28,1E4)
                ax[1].set_ylim(1E4, 1E7)


                plt.show()
                plt.close()






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

