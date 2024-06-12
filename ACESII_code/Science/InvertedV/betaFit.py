# --- betaFit.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: using the data from diffNFlux_fitting we can generate distributions at various altitudes
# to see which height our data most closely matches to


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import matplotlib.pyplot as plt
import numpy as np
import itertools
import os
from ACESII_code.myImports import *
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from ACESII_code.class_var_func import q0,m_e
from ACESII_code.Science.InvertedV.Evans_class_var_funcs import velocitySpace_to_PitchEnergySpace, loadDiffNFluxData,calc_DistributionMapping
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---



#################
# --- TOGGLES ---
#################


# --- Model Distribution Toggles ---
useKaepplerData = False
N = 1000 # velcoity space grid density
modifyInitialBeam = False
beamPitchThreshold = 50
# beta_model = [3,4,5,6,7,8,9,10] # If I wanted distance values between 1-2 R_E, then beta should be between 3 to 24
beta_model = [3] # If I wanted distance values between 1-2 R_E, then beta should be between 3 to 24

# --- mapped distribution to different betas ---
Plot_modelDistribution = False

# --- beam widening fits for beta ---
Plot_beamWidthBetaFit = True
betaIdxToPlot = 0 # Even though I'm trying to find beta, I still might want to plot it in order to compare

# --- beta pitch slice comparison plots ---
Plot_BetaSliceComparison = False

# --- output data as EEPAA CDF ---
outputCDF = False






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
mycmap = apl_rainbow_black0_cmap()





##########################
# --- LOADING THE DATA ---
##########################
prgMsg('Loading Data')
# note: these values come from the pitch = 10deg fit
if useKaepplerData:
    # ACESI setup
    ACESI_EEPAA_path = r'C:\Data\ACESII\science\invertedV\ACESI_EEPAA.cdf'
    data_dict_ACESI = loadDictFromFile(ACESI_EEPAA_path)
    diffEFlux = data_dict_ACESI['diff_flux'][0]
    Epoch, Energy, Pitch = data_dict_ACESI['Epoch'][0], data_dict_ACESI['energy_cal'][0][0], data_dict_ACESI['PA_bin'][0][0]
    diffNFlux = deepcopy(diffEFlux)
    Energy = np.array(Energy)

    # Convert Kaeppler's diffEFlux to diffNFlux
    ranges = [range(len(diffEFlux)), range(len(diffEFlux[0])), range(len(diffEFlux[0][0]))]
    for tme, ptch, engy in itertools.product(*ranges):
        val = diffEFlux[tme][ptch][engy]
        if val >= 0:
            diffNFlux[tme][ptch][engy] = val/Energy[engy]

    EnergyBins = Energy
    PitchBins = np.array([-180 + i * 15 for i in range(24 + 1)])
    compareThesePitches = [2,4,6,8,10]
    paramSet = 0
else:
    diffNFlux, Epoch, Energy, Pitch = loadDiffNFluxData()
    EnergyBins = Energy
    compareThesePitches = [2,4,6,8,10]

    PitchBins = np.array([-180 + i * 10 for i in range(36 + 1)])
    paramSet = 4


modelParams = [['2009-01-29 09:54:57:673000', 1.25, 1150, 3400],# kappler
               ['2022-11-20 17:25:01:312207', 1.5, 800, 2000],  # Evans 1974
               ['2022-11-20 17:25:01:312207', 3.799273, 108.904, 235.2], # our data
               ['2022-11-20 17:25:01:412211', 2.957, 124.599, 275],
               ['2022-11-20 17:25:01:612000', 4.748186, 137.043, 172.64] # OUR DATA: the GOOD SLICE
               ] # format: time, density [cm^-3], temperature [eV], [potential eV]

Done(start_time)


#################################
# --- --- --- --- --- --- --- ---
# --- GENERATE THE MODEL DATA ---
# --- --- --- --- --- --- --- ---
#################################
# Choose the input beam paramaters
model_T = modelParams[paramSet][2]  # in eV
model_n = modelParams[paramSet][1]  # in cm^-3
model_V0 = modelParams[paramSet][3]

# Containers to store data
VparaGrids_iono_model = []
VperpGrids_iono_model = []
distFunc_model = []
diffNFluxGrids_iono_model=[]
Epoch_model =[dt.datetime.strptime(st[0], "%Y-%m-%d %H:%M:%S:%f") for st in modelParams]

for betaChoice in beta_model:

    # --- Define a grid of intial velocities ---
    Vperp_gridVals = np.linspace(-np.sqrt(2*Energy.max()*q0/m_e), np.sqrt(2*Energy.max()*q0/m_e), N)
    Vpara_gridVals = np.linspace(0, np.sqrt(2*Energy.max()*q0/m_e), N)

    # --- Calculate the grids of the accelerate and mirror_mapped distributions ---
    prgMsg(rf'Creating Model Data for beta = {betaChoice}')
    distGrid,VperpGrid, VparaGrid, diffNFluxGrid, VperpGrid_Accel, VparaGrid_Accel, diffNFluxGrid_Accel, VperpGrid_iono, VparaGrid_iono, diffNFluxGrid_iono = calc_DistributionMapping(Vperp_gridVals=Vperp_gridVals,
                                                                                                                                                                              Vpara_gridVals=Vpara_gridVals,
                                                                                                                                                                              model_T=model_T,
                                                                                                                                                                              model_n=model_n,
                                                                                                                                                                              model_V0=model_V0,
                                                                                                                                                                              beta=betaChoice,
                                                                                                                                                                             modifyInitialBeam=modifyInitialBeam,
                                                                                                                                                                              beamPitchThreshold=beamPitchThreshold)
    # Store the data for various beta values
    VparaGrids_iono_model.append(VparaGrid_iono)
    VperpGrids_iono_model.append(VperpGrid_iono)
    distFunc_model.append(distGrid)
    diffNFluxGrids_iono_model.append(diffNFluxGrid_iono)

    Done(start_time)






#####################################
# --- Distribution Function Plots ---
#####################################
if Plot_modelDistribution:

    # remove the old files in the folder
    for file in glob(r'C:\Data\ACESII\science\invertedV\betaMappingPlots\*.png*'):
        os.remove(file)


    for betaIdx,betaChoice in enumerate(beta_model):

        prgMsg('Creating the Plot')

        # --- Plot it ---
        fig, ax = plt.subplots(3, 2)
        fig.set_size_inches(figure_width, figure_height)

        titles = ['Plasma Sheet Model','Accelerated', 'Observed Ionosphere Model']
        fig.suptitle(rf'$\beta$ = {betaChoice}')
        for k in range(2):

            if k == 0:
                grids = [[VperpGrid, VparaGrid, distGrid], [VperpGrid_Accel, VparaGrid_Accel, distGrid], [VperpGrids_iono_model[betaIdx], VparaGrids_iono_model[betaIdx], distGrid]]
                vmin, vmax = 1E-22, 1E-14
                cbarLabel = 'Distribution Function'

            else:
                grids = [[VperpGrid, VparaGrid, diffNFluxGrid], [VperpGrid_Accel, VparaGrid_Accel, diffNFluxGrid_Accel], [VperpGrids_iono_model[betaIdx], VparaGrids_iono_model[betaIdx], diffNFluxGrids_iono_model[betaIdx]]]
                vmin, vmax = 1E2, 1E8
                cbarLabel = 'diff_N_Flux'

            for i in [0, 1, 2]:
                cmap = ax[i, k].pcolormesh(grids[i][0]/(1E7), grids[i][1]/(1E7), grids[i][2], cmap=mycmap, norm='log', vmin=vmin, vmax=vmax)
                cbar = plt.colorbar(cmap, ax=ax[i,k])
                cbar.set_label(cbarLabel)
                ax[i, k].set_ylabel('Vpara')
                ax[i, k].set_xlabel('Vperp')
                ax[i, k].set_ylim(0, 8)
                ax[i, k].set_xlim(-8, 8)
                ax[i, k].invert_yaxis()
                ax[i, k].set_title(titles[i])

                if i in [1, 2]:
                    # plot the 110 deg line
                    ax[i, k].axhline(np.sqrt(2*model_V0*q0/m_e)/(1000*10000),color='red', label='$V_{0}$'+f'= {model_V0} eV')
                    ax[i, k].legend()
        plt.tight_layout()
        plt.savefig(rf'C:\Data\ACESII\science\invertedV\betaMappingPlots\BetaFit_{betaChoice}.png')
        plt.close()
        Done(start_time)



###########################################
# --- BETA PITCH SLICE COMPARISON PLOTS ---
###########################################
if Plot_BetaSliceComparison:
    prgMsg('Creating Beta Slice Comparison Plots')
    # remove the old files in the folder
    for file in glob(r'C:\Data\ACESII\science\invertedV\betaComparisonPlots\*.png*'):
        os.remove(file)


    for betaIdx in range(len(distFunc_model)):

        diffNFlux_model_pE, EnergyGrid, PitchGrid = velocitySpace_to_PitchEnergySpace(EnergyBins=EnergyBins,
                                                                             PitchBins=PitchBins,
                                                                             VperpGrid=VperpGrids_iono_model[betaIdx],
                                                                             VparaGrid=VparaGrids_iono_model[betaIdx],
                                                                             ZGrid=diffNFluxGrids_iono_model[betaIdx])


        # --- now create some comparison plots between the real data and your model ---

        # get the real data at the timeslice
        diffNFluxSlice = diffNFlux[np.abs(Epoch - Epoch_model[paramSet]).argmin()]

        for ptchIdxval in compareThesePitches:

            realData = diffNFluxSlice[ptchIdxval][:]
            closestPitch = np.abs(PitchBins-Pitch[ptchIdxval]).argmin()
            modelData = diffNFlux_model_pE[closestPitch][:]

            fig, ax = plt.subplots()
            fig.suptitle(rf'$\beta$={beta_model[betaIdx]}' + f'\n Pitch = {Pitch[ptchIdxval]}' + f'\n{Epoch_model[paramSet]}')
            ax.plot(Energy,realData,color='black',label='real Data',marker='.')
            ax.plot(EnergyBins,modelData,color='red',label='model',marker='.')
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_ylim(1E4,1E7)
            ax.set_xlim(1E1,1E4)
            ax.legend()
            plt.savefig(rf'C:\Data\ACESII\science\invertedV\betaComparisonPlots\Beta_{beta_model[betaIdx]}_pitch_{Pitch[ptchIdxval]}.png')
            plt.close()

    outputPath = rf'C:\Data\ACESII\science\invertedV\BetaFit_data.cdf'
    Done(start_time)


#########################################
# --- BETA BEAM PITCH WIDENNING PLOTS ---
#########################################
if Plot_beamWidthBetaFit:
    prgMsg('Creating Beta Slice Comparison Plots')

    # remove the old files in the folder
    for file in glob(r'C:\Data\ACESII\science\invertedV\BeamWidth\*.png*'):
        os.remove(file)

    # Get the model data (non-mirrored) in Pitch-Energy Space
    diffNFlux_model_pE, EnergyGrid, PitchGrid = velocitySpace_to_PitchEnergySpace(EnergyBins=EnergyBins,
                                                                         PitchBins=PitchBins,
                                                                         VperpGrid=VperpGrid_Accel,
                                                                         VparaGrid=VparaGrid_Accel,
                                                                         ZGrid=diffNFluxGrid_Accel)

    diffNFlux_model_mirrored_pE, EnergyGrid_mirrored, PitchGrid_mirrored = velocitySpace_to_PitchEnergySpace(EnergyBins=EnergyBins,
                                                                                  PitchBins=PitchBins,
                                                                                  VperpGrid=VperpGrids_iono_model[betaIdxToPlot],
                                                                                  VparaGrid=VparaGrids_iono_model[betaIdxToPlot],
                                                                                  ZGrid=diffNFluxGrids_iono_model[betaIdxToPlot])

    # get the maximum angle for each pitch angle in the model beam
    beamPitchModel = [arr.max() for arr in diffNFlux_model_pE]
    BeamEnergyModel = [EnergyGrid[0] for arr in EnergyGrid]
    print(beamPitchModel)
    print(BeamEnergyModel)


    # get the real data at the timeslice
    diffNFluxSlice = np.array(diffNFlux[np.abs(Epoch - Epoch_model[paramSet]).argmin()])

    # define alpha(E)_max for the data (the maximum beam angle for a given energy)
    BeamEnergyChoice = [725.16, 621.29, 532.3, 456.06, 390.74, 334.77, 286.82, 245.74, 210.54, 180.39, 154.55, 132.41, 113.45, 97.2, 83.28, 71.35, 61.13, 52.37, 44.87, 38.44, 32.94, 28.22]
    BeamPitchChoice = [90-2.5 for i in range(len(BeamEnergyChoice))]

    # calculate the


    # --- --- --- --- --- ----
    # Make the Beam-Width Plot
    # --- --- --- --- --- ----
    fig, ax = plt.subplots(4)
    fig.set_size_inches(figure_width, 2*figure_height)

    # Accelerated Model Data
    cmapObj=ax[0].pcolormesh(PitchGrid, EnergyGrid, np.array(diffNFlux_model_pE),vmin=1E4, vmax=5E6, norm='log',cmap=mycmap)
    ax[0].set_yscale('log')
    ax[0].set_ylim(100, 2000)
    ax[0].set_xlim(0, 90)
    cbar = plt.colorbar(cmapObj)

    # Mirrored Model Data
    cmapObj = ax[1].pcolormesh(PitchGrid_mirrored, EnergyGrid_mirrored, np.array(diffNFlux_model_mirrored_pE), vmin=1E4, vmax=5E6, norm='log', cmap=mycmap)
    ax[1].set_yscale('log')
    ax[1].set_ylim(100, 2000)
    ax[1].set_xlim(0, 90)
    cbar = plt.colorbar(cmapObj)

    # Real Data + Beam Choice
    cmapObj=ax[2].pcolormesh(Pitch, Energy, diffNFluxSlice.T, vmin=1E4, vmax=5E6, norm='log',cmap=mycmap)
    ax[2].scatter(BeamPitchChoice, BeamEnergyChoice, color='black')
    ax[2].set_yscale('log')
    ax[2].set_ylim(100, 2000)
    ax[2].set_xlim(0, 90)
    cbar = plt.colorbar(cmapObj)

    # Beta fit from Beam Choice
    ax[3].plot([1, 2, 3], [1, 2, 3])
    ax[3].set_ylabel(r'$\sin^{2}(\alpha_{I})$',fontsize=Label_FontSize)
    ax[3].set_xlabel(r'$\sin^{2}(\alpha_{Model})$',fontsize=Label_FontSize)
    plt.tight_layout()
    plt.savefig(rf'C:\Data\ACESII\science\invertedV\BeamWidth\BeamWidthAnalysis.png')


    Done(start_time)


############################
# --- OUTPUT DATA AS CDF ---
############################
if outputCDF:

    prgMsg('Outputting Data')
    rocketAttrs, b, c = ACES_mission_dicts()
    exampleAttrs = {'LABLAXIS': None, 'DEPEND_0': None,
                'DEPEND_1': None,
                'DEPEND_2': None,
                'FILLVAL': rocketAttrs.epoch_fillVal,
                'FORMAT': 'E12.2',
                'UNITS': None,
                'VALIDMIN': None,
                'VALIDMAX': None,
                'VAR_TYPE': 'data', 'SCALETYP': 'linear'}


    data_dict_output = {'Pitch_Angle':[PitchBins, deepcopy(exampleAttrs)],
                        'Energy':[EnergyBins,deepcopy(exampleAttrs)],
                        'Epoch':[ np.array([pycdf.lib.datetime_to_tt2000(Epoch_model[paramSet])]),deepcopy(exampleAttrs)]}

    data_dict_output['Pitch_Angle'][1]['UNITS'] = 'deg'
    data_dict_output['Energy'][1]['UNITS'] = 'eV'


    for betaIdx in range(len(distFunc_model)):
        diffNFlux_model_pE, EnergyGrid, PitchGrid = velocitySpace_to_PitchEnergySpace(EnergyBins=EnergyBins,
                                                                                      PitchBins=PitchBins,
                                                                                      VperpGrid=VperpGrids_iono_model[betaIdx],
                                                                                      VparaGrid=VparaGrids_iono_model[ betaIdx],
                                                                                      ZGrid=diffNFluxGrids_iono_model[ betaIdx])

        # write it all out to datafile
        data_dict_output = {**data_dict_output, **{f'diffNFlux_beta{beta_model[betaIdx]}':[np.array([diffNFlux_model_pE]),deepcopy(exampleAttrs)]}}
        data_dict_output[f'diffNFlux_beta{beta_model[betaIdx]}'][1]['UNITS']='cm^-2 str^-1 s^-1 eV^-1'
        data_dict_output[f'diffNFlux_beta{beta_model[betaIdx]}'][1]['DEPEND_0'] = 'Epoch'
        data_dict_output[f'diffNFlux_beta{beta_model[betaIdx]}'][1]['DEPEND_1'] = 'Pitch_Angle'
        data_dict_output[f'diffNFlux_beta{beta_model[betaIdx]}'][1]['DEPEND_2'] = 'Energy'

    outputCDFdata(outputPath, data_dict_output, instrNam='EEPAA')
    Done(start_time)