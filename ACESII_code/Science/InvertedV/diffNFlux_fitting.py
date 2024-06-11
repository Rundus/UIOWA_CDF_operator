# --- diffNFlux_fitting.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: using the method outline in Kaeppler's thesis, we can fit inverted-V distributions
# to get estimate the magnetospheric temperature, density and electrostatic potential that accelerated
# our particles


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from ACESII_code.class_var_func import EpochTo_T0_Rocket, q0,m_e
from ACESII_code.Science.InvertedV.Evans_class_var_funcs import calc_diffNFlux,dist_Maxwellian
from functools import partial
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---


print(color.UNDERLINE + f'diffNFlux_fitting' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# --- Physics Toggles ---
invertedV_targetTimes = [dt.datetime(2022,11, 20, 17, 24, 56, 000000), dt.datetime(2022,11,20,17,25,3,000000)]
invertedV_TargetTimes_data = [[dt.datetime(2022,11, 20, 17, 25, 1, 300000), dt.datetime(2022,11,20,17,25,1,800000)]]

# --- Plot toggles - General ---
figure_width = 10 # in inches
figure_height =8 # in inches
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

# ---  Density, Temperature and Potential ---
invertedV_fitDensityTempPotential = False
wPitchToFit = 2
countNoiseLevel = 1


##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
prgMsg('Loading Data')
from ACESII_code.Science.InvertedV.Evans_class_var_funcs import loadDiffNFluxData
diffNFlux,Epoch,Energy,Pitch = loadDiffNFluxData()
Done(start_time)

###############################
# --- --- --- --- --- --- --- -
# --- FIT THE DISTRIBUTIONS ---
# --- --- --- --- --- --- --- -
###############################


# define my function at the specific pitch angle of interest

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


        pairs = [f'({x},{y})' for x,y in zip(xData_fit,yData_fit)]


        fittedX = np.linspace(xData_fit.min(), xData_fit.max(), 100)
        fittedY = fitFuncAtPitch(fittedX,*params)

        # --- Calculate ChiSquare ---

        for i in range(len(xData_fit)):
            print(fitFuncAtPitch(xData_fit[i],*params),yData_fit[i],xData_fit[i],(fitFuncAtPitch(xData_fit[i],*params) - yData_fit[i])**2)

        ChiSquare = (1/2)*sum([(fitFuncAtPitch(xData_fit[i],*params) - yData_fit[i])**2 / (yData_fit[i]**2) for i in range(len(xData_fit))])


        # Plot the result to see
        fig, ax = plt.subplots(2)
        fig.set_size_inches(figure_width, figure_height)
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
