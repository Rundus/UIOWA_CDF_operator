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
from scipy.interpolate import griddata
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---


print(color.UNDERLINE + f'diffNFlux_fitting' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# --- Physics Toggles ---
datasetReduction_TargetTime = [dt.datetime(2022,11, 20, 17, 24, 50, 000000), dt.datetime(2022,11,20,17,25,15,000000)]
targetVar = [datasetReduction_TargetTime, 'Epoch']

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

# --- Beta fit via model params ---
invertedV_betafit = True
outputCDF = True


##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
prgMsg('Loading Data')

rocketAttrs, b, c = ACES_mission_dicts()

# EEPAA Distribution Data
inputEEPAA_diffFlux_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0]
data_dict_diffFlux = loadDictFromFile(inputEEPAA_diffFlux_high, targetVar=targetVar, wKeys_Reduce=['Differential_Energy_Flux','Differential_Number_Flux', 'Epoch'])

# Define the data
diffEFlux = deepcopy(data_dict_diffFlux['Differential_Energy_Flux'][0])
diffNFlux = deepcopy(data_dict_diffFlux['Differential_Number_Flux'][0])
Epoch = deepcopy(data_dict_diffFlux['Epoch'][0])
Energy = deepcopy(data_dict_diffFlux['Energy'][0])
Pitch = deepcopy(data_dict_diffFlux['Pitch_Angle'][0])
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
    deltaT = (count_interval) - (countNoiseLevel * rocketAttrs.deadtime[0])
    fluxVal = (countNoiseLevel) / (geo_factor[0] * deltaT)
    distFunc_NoiseCount[engy] = ((cm_to_m * m_e / (q0 * Energy[engy])) ** 2) * (fluxVal / 2)

# Time since launch
LaunchDateTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 000000))
Epoch_timeSince = EpochTo_T0_Rocket(InputEpoch=Epoch, T0=LaunchDateTime)
invertedV_timeSince = EpochTo_T0_Rocket(InputEpoch=invertedV_targetTimes, T0=LaunchDateTime)
Done(start_time)
###############################
# --- --- --- --- --- --- --- -
# --- FIT THE DISTRIBUTIONS ---
# --- --- --- --- --- --- --- -
###############################
from ACESII_code.Science.InvertedV.fitFunctions import diffNFlux_for_mappedMaxwellian,calc_diffNFlux,dist_Maxwellian




if invertedV_fitDensityTempPotential:



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

if invertedV_betafit:

    # --- toggles ---
    beta_model = [5,10,15,20,25]
    paramSet = 2

    # note: these values come from the pitch = 10deg fit
    modelParams = [['2022-11-20 17:25:01:312207', 1.25, 1150, 3400],# kappler
                   ['2022-11-20 17:25:01:312207', 1.5, 800, 2000],  # Evans 1974
                   ['2022-11-20 17:25:01:312207', 3.799273, 108.904, 235.2], # our data
                   ['2022-11-20 17:25:01:362211', 5.365, 146.9, 201],
                   ['2022-11-20 17:25:01:412211', 2.957, 124.599, 275],
                   ['2022-11-20 17:25:01:462208', 3.3039, 115.86, 235.4]] # format: time, density [cm^-3], temperature [eV], [potential eV]

    Vpar_model = []
    Vperp_model = []
    distFunc_model = []
    Epoch_model =[dt.datetime.strptime(st[0], "%Y-%m-%d %H:%M:%S:%f") for st in modelParams]

    #################################
    # --- GENERATE THE MODEL DATA ---
    #################################
    prgMsg('Creating Model Data')
    for betaChoice in beta_model:

        # betaChoice = 7  # unitless
        model_T = modelParams[paramSet][2]  # in eV
        model_n = modelParams[paramSet][1]  # in cm^-3
        model_V0 = modelParams[paramSet][3]

        # --- Define a grid a velocities (static) ---
        N = 1000
        Vperp_gridVals = np.linspace(-np.sqrt(2*Energy.max()*q0/m_e), np.sqrt(2*Energy.max()*q0/m_e), N)
        Vpara_gridVals = np.linspace(0, np.sqrt(2*Energy.max()*q0/m_e), N)
        VperpGrid, VparaGrid = np.meshgrid(Vperp_gridVals, Vpara_gridVals)
        distGrid = dist_Maxwellian(VperpGrid, VparaGrid, n=model_n, T=model_T)
        diffNFluxGrid = calc_diffNFlux(Vperp_gridVals,Vpara_gridVals,distGrid)

        # --- Determine the Accelerated Velocities ---
        Vperp_gridVals_Accel = Vperp_gridVals
        Vpar_gridVals_Accel = np.array([np.sqrt(val**2 + 2*model_V0*q0/m_e) for val in Vpara_gridVals])
        VperpGrid_Accel,VparGrid_Accel = np.meshgrid(Vperp_gridVals_Accel,Vpar_gridVals_Accel)
        diffNFluxGrid_Accel = calc_diffNFlux(Vperp_gridVals_Accel,Vpar_gridVals_Accel,distGrid)

        # --- Determine the new velocities at different beta ---
        VperpArray_magsph = VperpGrid_Accel.flatten()
        VparaArray_magsph = VparGrid_Accel.flatten()
        distFuncArray= distGrid.flatten()
        VperpArray_iono = np.array([ np.sqrt(betaChoice)*val for val in VperpArray_magsph])
        VparaArray_iono_sqrd = np.array([ Vpar_magsph**2 + (1 - betaChoice)*(Vper_magsph**2) for Vper_magsph, Vpar_magsph in zip(VperpArray_magsph,VparaArray_magsph)])
        VparaArray_iono = np.array([np.sqrt(val) if val >=0 else -1*np.sqrt(np.abs(val)) for val in VparaArray_iono_sqrd ])

        if outputCDF:
            Vpar_model.append(VparaArray_iono)
            Vperp_model.append(VperpArray_iono)
            distFunc_model.append(distFuncArray)

        # make the grids
        VperpGrid_iono = VperpArray_iono.reshape(N,N)
        VparGrid_iono = VparaArray_iono.reshape(N,N)
        diffNFluxGrid_iono = calc_diffNFlux(VperpGrid_iono,VparGrid_iono,distGrid)

        # --- Determine the Distribution at a beta ratio ---

        # --- Plot it ---
        fig, ax = plt.subplots(3,2)
        fig.set_size_inches(figure_width, figure_height)

        titles = ['Plasma Sheet Model','Accelerated', 'Observed Ionosphere Model']
        fig.suptitle(rf'$\beta$ = {betaChoice}')
        for k in range(2):

            if k == 0:
                grids = [[VperpGrid, VparaGrid, distGrid], [VperpGrid_Accel, VparGrid_Accel, distGrid], [VperpGrid_iono, VparGrid_iono, distGrid]]
                vmin, vmax = 1E-22, 1E-14
                cbarLabel = 'Distribution Function'

            else:
                grids = [[VperpGrid, VparaGrid, diffNFluxGrid], [VperpGrid_Accel, VparGrid_Accel, diffNFluxGrid_Accel], [VperpGrid_iono, VparGrid_iono, diffNFluxGrid_iono]]
                vmin, vmax = 1E4, 1E7
                cbarLabel = 'diff_N_Flux'


            for i in [0, 1, 2 ]:
                cmap = ax[i,k].pcolormesh(grids[i][0]/(1E7), grids[i][1]/(1E7), grids[i][2], cmap=mycmap, norm='log', vmin=vmin, vmax=vmax)
                cbar = plt.colorbar(cmap, ax=ax[i,k])
                cbar.set_label(cbarLabel)
                ax[i,k].set_ylabel('Vpara')
                ax[i,k].set_xlabel('Vperp')
                ax[i,k].set_ylim(0, 3)
                ax[i,k].set_xlim(-3, 3)
                ax[i,k].invert_yaxis()
                ax[i,k].set_title(titles[i])

                if i in [1,2]:
                    ax[i,k].axhline(np.sqrt(2*model_V0*q0/m_e)/(1000*10000),color='red', label='$V_{0}$'+f'= {model_V0} eV')

                    # plot the 110 deg line
                    pitchBarVal = 110
                    EmagVal = (2*10000*q0/m_e)
                    ax[i,k].plot([0, EmagVal*np.sin(np.radians(pitchBarVal))],[0, EmagVal*np.cos(np.radians(pitchBarVal))], color='black',label='$110^{\circ}$')
                    ax[i,k].plot([0,-1*EmagVal * np.sin(np.radians(pitchBarVal))], [0, EmagVal * np.cos(np.radians(pitchBarVal))], color='black')
                    ax[i,k].legend()
        plt.tight_layout()
        plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\PlotExtra\BetaFit\BetaFit_{betaChoice}.png',dpi=dpi)
        plt.close()
    Done(start_time)

    if outputCDF:

        prgMsg('Outputting Data')

        exampleAttrs = {'LABLAXIS': None, 'DEPEND_0': None,
                    'DEPEND_1': None,
                    'DEPEND_2': None,
                    'FILLVAL': rocketAttrs.epoch_fillVal,
                    'FORMAT': 'E12.2',
                    'UNITS': None,
                    'VALIDMIN': None,
                    'VALIDMAX': None,
                    'VAR_TYPE': 'data', 'SCALETYP': 'linear'}


        # ACESI setup
        # ACESI_EEPAA_path = r'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\PlotExtra\BetaFit\ACESI_EEPAA.cdf'
        # data_dict_ACESI = loadDictFromFile(ACESI_EEPAA_path)
        # EnergyBins = np.flip(np.array(data_dict_ACESI['energy_cal'][0][0]))
        # PitchBins = np.array([-180 + i * 15 for i in range(24 + 1)])

        # ACESII setup
        EnergyBins = Energy
        PitchBins = np.array([-180 + i * 10 for i in range(36 + 1)])

        data_dict_output = {'Pitch_Angle':[PitchBins, deepcopy(exampleAttrs)],
                            'Energy':[EnergyBins,deepcopy(exampleAttrs)],
                            'Epoch':[ np.array([pycdf.lib.datetime_to_tt2000(Epoch_model[paramSet])]),deepcopy(exampleAttrs)]}

        data_dict_output['Pitch_Angle'][1]['UNITS'] = 'deg'
        data_dict_output['Energy'][1]['UNITS'] = 'eV'


        for betaIdx in range(len(distFunc_model)):

            diffNFlux_model = [  [ [[] for engy in range(len(EnergyBins))] for ptch in range(len(PitchBins))] for tme in range(len([0]))]

            calcEnergies = [0.5*m_e*(perp**2 + par**2)/q0 for perp,par in zip(Vperp_model[betaIdx],Vpar_model[betaIdx])]
            calcPitch = [(180/np.pi)*np.arctan2(perp, par) for perp,par in zip(Vperp_model[betaIdx],Vpar_model[betaIdx])]

            # assign the values to the now calculated diffNFlux n
            for i in range(len(calcEnergies)):

                engyIdx = np.abs(EnergyBins-calcEnergies[i]).argmin()
                ptchIdx = np.abs(PitchBins-calcPitch[i]).argmin()

                diffNFlux_model[0][ptchIdx][engyIdx].append(   (2 * calcEnergies[i]) * ((q0 / m_e) ** 2)*distFunc_model[betaIdx][i])

            # flatten the values in the diffnFlux new array
            for tme in range(len([0])):
                for ptch in range(len(PitchBins)):
                    for engy in range(len(EnergyBins)):
                        try:
                            diffNFlux_model[tme][ptch][engy] = sum(diffNFlux_model[tme][ptch][engy])/len(diffNFlux_model[tme][ptch][engy])
                            # diffNFlux_model[tme][ptch][engy] = sum(diffNFlux_model[tme][ptch][engy])
                        except:
                            diffNFlux_model[tme][ptch][engy] = sum(diffNFlux_model[tme][ptch][engy])

            # write it all out to datafile
            data_dict_output = {**data_dict_output, **{f'diffNFlux_beta{beta_model[betaIdx]}':[np.array(diffNFlux_model),deepcopy(exampleAttrs)]}}
            data_dict_output[f'diffNFlux_beta{beta_model[betaIdx]}'][1]['UNITS']='cm^-2 str^-1 s^-1 eV^-1'
            data_dict_output[f'diffNFlux_beta{beta_model[betaIdx]}'][1]['DEPEND_0'] = 'Epoch'
            data_dict_output[f'diffNFlux_beta{beta_model[betaIdx]}'][1]['DEPEND_1'] = 'Pitch_Angle'
            data_dict_output[f'diffNFlux_beta{beta_model[betaIdx]}'][1]['DEPEND_2'] = 'Energy'


            # --- now create some comparison plots between the real data and your model ---

            # get the real data at the timeslice
            diffNFluxSlice = diffNFlux[np.abs(Epoch - Epoch_model[paramSet]  ).argmin()]
            compareThesePitches = [2,3, 4,5, 6,7, 8,10]

            for ptchIdxval in compareThesePitches:

                realData = diffNFluxSlice[ptchIdxval][:]
                closestPitch = np.abs(PitchBins-Pitch[ptchIdxval]).argmin()
                modelData = diffNFlux_model[0][closestPitch][:]

                fig, ax = plt.subplots()
                fig.suptitle(rf'$\beta$={beta_model[betaIdx]}' + f'\n Pitch = {Pitch[ptchIdxval]}' + f'\n{Epoch_model[paramSet]}')
                ax.plot(Energy,realData,color='black',label='real Data',marker='.')
                ax.plot(EnergyBins,modelData,color='red',label='model',marker='.')
                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_ylim(1E4,1E7)
                ax.set_xlim(1E1,1E4)
                ax.legend()
                plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\PlotExtra\BetaFit\betaComparisonPlots\Beta_{beta_model[betaIdx]}_pitch_{Pitch[ptchIdxval]}.png')
                plt.close()

        outputPath = rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\PlotExtra\BetaFit\BetaFit_data.cdf'
        outputCDFdata(outputPath, data_dict_output, instrNam='EEPAA')
        Done(start_time)