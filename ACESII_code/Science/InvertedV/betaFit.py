# --- betaFit.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: using the data from diffNFlux_fitting we can generate distributions at various altitudes
# to see which height our data most closely matches to


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


##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
prgMsg('Loading Data')
from ACESII_code.Science.InvertedV.Evans_class_var_funcs import loadDiffNFluxData
diffNFlux,Epoch,Energy,Pitch = loadDiffNFluxData()
Done(start_time)

##########################
# --- BETA FIT TOGGLES ---
##########################
# mapped distribution to different betas
beta_model = [1,2,4]
paramSet = 2
N = 1000 # velcoity space grid density
pitchBarVal = 110 # pitch angle to plot on graphs

# output data as if EEPAA detected it
outputCDF = True
EnergyBins = Energy
PitchBins = np.array([-180 + i * 10 for i in range(36 + 1)]) #
compareThesePitches = [2,3, 4,5, 6,7, 8,10]


# note: these values come from the pitch = 10deg fit
# modelParams = [['2022-11-20 17:25:01:312207', 1.25, 1150, 3400],# kappler
#                ['2022-11-20 17:25:01:312207', 1.5, 800, 2000],  # Evans 1974
#                ['2022-11-20 17:25:01:312207', 3.799273, 108.904, 235.2], # our data
#                ['2022-11-20 17:25:01:362211', 5.365, 146.9, 201],
#                ['2022-11-20 17:25:01:412211', 2.957, 124.599, 275],
#                ['2022-11-20 17:25:01:462208', 3.3039, 115.86, 235.4]] # format: time, density [cm^-3], temperature [eV], [potential eV]

modelParams = [['2022-11-20 17:25:01:312207', 1.25, 1150, 3400],# kappler
               ['2022-11-20 17:25:01:312207', 1.5, 800, 2000],  # Evans 1974
               ['2022-11-20 17:25:01:412211', 2.957, 124.599, 275]] # OUR DATA: the GOOD SLICE

#################################
# --- GENERATE THE MODEL DATA ---
#################################
prgMsg('Creating Model Data')

Vpar_model = []
Vperp_model = []
distFunc_model = []
Epoch_model =[dt.datetime.strptime(st[0], "%Y-%m-%d %H:%M:%S:%f") for st in modelParams]

for betaChoice in beta_model:

    # betaChoice = 7  # unitless
    model_T = modelParams[paramSet][2]  # in eV
    model_n = modelParams[paramSet][1]  # in cm^-3
    model_V0 = modelParams[paramSet][3]

    # --- Define a grid a velocities (static) ---
    Vperp_gridVals = np.linspace(-np.sqrt(2*Energy.max()*q0/m_e), np.sqrt(2*Energy.max()*q0/m_e), N)
    Vpara_gridVals = np.linspace(0, np.sqrt(2*Energy.max()*q0/m_e), N)
    VperpGrid, VparaGrid = np.meshgrid(Vperp_gridVals, Vpara_gridVals)
    distGrid = dist_Maxwellian(VperpGrid, VparaGrid, n=model_n, T=model_T)
    diffNFluxGrid = calc_diffNFlux(VperpGrid,VparaGrid,distGrid)

    # --- Determine the Accelerated Velocities ---
    Vperp_gridVals_Accel = Vperp_gridVals
    Vpar_gridVals_Accel = np.array([np.sqrt(val**2 + 2*model_V0*q0/m_e) for val in Vpara_gridVals])
    VperpGrid_Accel,VparGrid_Accel = np.meshgrid(Vperp_gridVals_Accel,Vpar_gridVals_Accel)
    diffNFluxGrid_Accel = calc_diffNFlux(VperpGrid_Accel,VparGrid_Accel,distGrid)

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
            vmin, vmax = 1E4, 4E6
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
                EmagVal = (2*10000*q0/m_e)
                ax[i,k].plot([0, EmagVal*np.sin(np.radians(pitchBarVal))],[0, EmagVal*np.cos(np.radians(pitchBarVal))], color='black',label='$110^{\circ}$')
                ax[i,k].plot([0,-1*EmagVal * np.sin(np.radians(pitchBarVal))], [0, EmagVal * np.cos(np.radians(pitchBarVal))], color='black')
                ax[i,k].legend()
    plt.tight_layout()
    plt.savefig(rf'C:\Data\ACESII\science\invertedV\BetaFit_{betaChoice}.png',dpi=dpi)
    plt.close()
Done(start_time)





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


    # ACESI setup
    # ACESI_EEPAA_path = r'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\PlotExtra\BetaFit\ACESI_EEPAA.cdf'
    # data_dict_ACESI = loadDictFromFile(ACESI_EEPAA_path)
    # EnergyBins = np.flip(np.array(data_dict_ACESI['energy_cal'][0][0]))
    # PitchBins = np.array([-180 + i * 15 for i in range(24 + 1)])

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
            plt.savefig(rf'C:\Data\ACESII\science\invertedV\betaComparisonPlots\Beta_{beta_model[betaIdx]}_pitch_{Pitch[ptchIdxval]}.png')
            plt.close()

    outputPath = rf'C:\Data\ACESII\science\invertedV\BetaFit_data.cdf'
    outputCDFdata(outputPath, data_dict_output, instrNam='EEPAA')
    Done(start_time)