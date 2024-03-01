# --- alfven_Eperp_Generator.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:
# Generates the entire E-Field for all time in the simulation and returns a variable that looks like:
# [
#   [Ez(x=0,t=0),Ez(x=1,t=0),Ez(x=2,t=0)...., Ez(x=len(Alt),t=0)],
#   [Ez(x=0,t=1),Ez(x=1,t=1),Ez(x=2,t=1)...., Ez(x=len(Alt),t=1)]
#   ,...]
import numpy as np

# --- imports ---
from ACESII_code.myImports import *
from ACESII_code.class_var_func import L2_ACES_Quick
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.simToggles import m_to_km, R_REF, GenToggles,EToggles
from ACESII_code.class_var_func import loadDictFromFile
start_time = time.time()

##################
# --- PLOTTING ---
##################
plot_Epara = False

################
# --- OUTPUT ---
################
outputData = True

# get the Eperp and plasma environment Profiles
data_dict_plasEvrn = loadDictFromFile(f'{GenToggles.simFolderPath}\plasmaEnvironment\plasmaEnvironment.cdf')
data_dict_Eperp = loadDictFromFile(rf'{GenToggles.simFolderPath}\alfvenWave\Eperp\alfvenWave_Eperp.cdf')


def alfvenEparaGenerator(outputData, **kwargs):
    plotBool = kwargs.get('showPlot', False)

    def EparaProfile(altRange, timeRange, plotBool):

        # get the profiles and flip them so we begin with high altitudes
        altRange = altRange
        lambdaPerp, kperp, skinDepth = data_dict_plasEvrn['lambdaPerp'][0],data_dict_plasEvrn['kperp'][0], data_dict_plasEvrn['skinDepth'][0]
        Eperp = data_dict_Eperp['Eperp'][0]

        # create the X dimension
        simXRange = np.linspace(0, EToggles.lambdaPerp0, EToggles.lambdaPerp_Rez)

        # create a meshgrid for determining dE_perp /dz
        Epara = []
        for tme, timeVal in tqdm(enumerate(timeRange)):

            spaceGrid = np.zeros(shape=(len(altRange),len(simXRange)))

            for x in range(len(simXRange)):
                for z in range(len(altRange)-1):
                    Eperp_n = Eperp[tme][z][x]
                    Eperp_n1 = Eperp[tme][z+1][x]
                    EperpGradVal = (Eperp_n1 - Eperp_n) / (altRange[z+1] - altRange[z])
                    spaceGrid[z][x] = kperp[z]*(skinDepth[z]**2) * EperpGradVal / (1 + (kperp[z]*skinDepth[z]**2))  # calculate the gradient in

            Epara.append(spaceGrid)

        Epara = np.array(Epara)

        # plotting
        if plotBool:
            from ACESII_code.class_var_func import prgMsg,Done
            import time
            start_time = time.time()
            prgMsg('Creating Epara Plots')
            import matplotlib.pyplot as plt
            plt.rcParams.update({'font.size': 22})
            for timeIndexChoice in range(len(timeRange)):
                fig,ax = plt.subplots()
                fig.set_figwidth(15)
                fig.set_figheight(15)
                ax.set_title('$E_{\parallel}$ Propogation vs Altitude \n' + f't = {timeRange[timeIndexChoice]}'+r'  $\tau_{0}$=' +f'{EToggles.tau0} s' +'\n' + '$\lambda_{\perp}$ =' + f'{EToggles.lambdaPerp0} [m],  ' + '$\omega_{wave}$ =' + f'{EToggles.waveFreq_Hz} [Hz]')
                cmap = ax.pcolormesh(simXRange/EToggles.lambdaPerp0,altRange/R_REF, Epara[timeIndexChoice]*1E6, cmap='bwr', vmin=-1*1E6,vmax=1E6)
                ax.set_xlabel('X Distance [$\lambda_{\perp}$]')
                ax.set_ylabel('Z Distance [$R_{E}$]')
                ax.grid(True)
                cbar = plt.colorbar(cmap)
                cbar.set_label('$|E_{\parallel}$| [uV/m]')
                plt.savefig(rf'{GenToggles.simFolderPath}\alfvenWave\Epara\plots\Epara_t{timeIndexChoice}.png')
                plt.close()
            Done(start_time)

        # create the output variable
        return np.array(Epara), simXRange

    # get all the variables
    Epara, simXRange = EparaProfile(altRange=GenToggles.simAlt, timeRange=GenToggles.simTime, plotBool=plotBool)

    ################
    # --- OUTPUT ---
    ################
    if outputData:
        prgMsg('Writing out Epara Data')

        # --- ACES II Flight/Integration Data ---
        wRocket = 4
        rocketAttrs, b, c = ACES_mission_dicts()
        globalAttrsMod = rocketAttrs.globalAttributes[wRocket - 4]
        globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
        ModelData = L2_ACES_Quick(wRocket - 4)

        # --- Construct the Data Dict ---
        exampleVar = {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                      'FORMAT': 'I5', 'UNITS': 'm', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                      'SCALETYP': 'linear', 'LABLAXIS': 'simAlt'}

        data_dict = {'Epara': [Epara, deepcopy(exampleVar)],
                     'simXRange': [simXRange, deepcopy(exampleVar)],
                     'simTime': [GenToggles.simTime, deepcopy(exampleVar)],
                     'simAlt': [GenToggles.simAlt, deepcopy(exampleVar)]}

        data_dict['Epara'][1]['UNITS'] = 'V/m'
        data_dict['Epara'][1]['LABLAXIS'] = 'Epara'
        data_dict['Epara'][1]['DEPEND_0'] = 'simTime'
        data_dict['Epara'][1]['DEPEND_1'] = 'simAlt'
        data_dict['Epara'][1]['DEPEND_2'] = 'simXRange'

        data_dict['simXRange'][1]['UNITS'] = 'm'
        data_dict['simXRange'][1]['LABLAXIS'] = 'simXRange'
        data_dict['simXRange'][1]['DEPEND_0'] = 'simXRange'

        data_dict['simTime'][1]['UNITS'] = 'seconds'
        data_dict['simTime'][1]['LABLAXIS'] = 'simTime'
        data_dict['simTime'][1]['DEPEND_0'] = 'simTime'

        data_dict['simAlt'][1]['UNITS'] = 'm'
        data_dict['simAlt'][1]['LABLAXIS'] = 'simAlt'
        data_dict['simAlt'][1]['DEPEND_0'] = 'simAlt'

        # --- output the data ---
        outputPath = rf'{GenToggles.simFolderPath}\alfvenWave\Epara\alfvenWave_Epara.cdf'
        outputCDFdata(outputPath, data_dict, ModelData, globalAttrsMod, 'simulation')
        Done(start_time)



#################
# --- EXECUTE ---
#################
alfvenEparaGenerator(outputData=outputData,showPlot=plot_Epara)