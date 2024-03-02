# --- alfven_Eperp_Generator.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:
# Generates the entire E-Field for all time in the simulation and returns a variable that looks like:
# [
#   [Ez(x=0,t=0),Ez(x=1,t=0),Ez(x=2,t=0)...., Ez(x=len(Alt),t=0)],
#   [Ez(x=0,t=1),Ez(x=1,t=1),Ez(x=2,t=1)...., Ez(x=len(Alt),t=1)]
#   ,...]

# --- imports ---
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.simToggles import m_to_km, R_REF, GenToggles,EToggles
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.plasmaEnvironment.plasmaEnvironment_Generator import generatePlasmaEnvironment
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.geomagneticField.geomagneticField_Generator import generateGeomagneticField
import time
import numpy as np
from itertools import product
from copy import deepcopy
from ACESII_code.class_var_func import prgMsg,Done, outputCDFdata, loadDictFromFile
start_time = time.time()

##################
# --- PLOTTING ---
##################
plot_Eperp = False

################
# --- OUTPUT ---
################
outputData = True


# --- Re-run the plasma environment and load the data ---
regenerateEnvironment = True
if regenerateEnvironment:
    generateGeomagneticField(outputData=True)
    generatePlasmaEnvironment(outputData=True)

data_dict_plasEvrn = loadDictFromFile(f'{GenToggles.simFolderPath}\plasmaEnvironment\plasmaEnvironment.cdf')


def alfvenEperpGenerator(outputData, **kwargs):


    # --- Eperp Makder ---
    def Eperp_generator(x, z, t, Vel, tau, KperpVal):

        # # the backside of the wave
        # if GenToggles.simAltHigh  > z > GenToggles.simAltHigh - EToggles.Z0_wave - Vel*(t-tau):
        #     EperpVal = 0
        #
        # # the middle of wave
        # elif EToggles.Z0_wave - Vel*(t-tau) > z > EToggles.Z0_wave -  Vel*t:
        #     amplitude = EToggles.Eperp0*np.sin(KperpVal*x/2)
        #     EperpVal = amplitude* (1 - np.cos((z - Vel*t) * (2*np.pi / (Vel*tau)) ) )
        #
        # # the front part of the wave
        # elif EToggles.Z0_wave -  Vel*t > z >  GenToggles.simAltLow:
        #     EperpVal = 0
        # else:
        #     EperpVal = 0

        # the middle of wave
        if EToggles.Z0_wave - Vel * (t - tau) > z > EToggles.Z0_wave - Vel * t:
            amplitude = EToggles.Eperp0 * np.sin(KperpVal * x / 2)
            EperpVal = amplitude * (1 - np.cos(((z - EToggles.Z0_wave) + Vel * t) * (2 * np.pi / (Vel * tau))))
        else:
            EperpVal = 0

        return EperpVal

    def EperpProfile(altRange, timeRange, **kwargs):
        plotBool = kwargs.get('showPlot', False)

        Eperp = []

        # get the profiles and flip them so we begin with high altitudes
        altRange = altRange
        lambdaPerp, kperp = data_dict_plasEvrn['lambdaPerp'][0],data_dict_plasEvrn['kperp'][0]
        alfSpdMHD = data_dict_plasEvrn['alfSpdMHD'][0]
        alfSpdInertial = data_dict_plasEvrn['alfSpdInertial'][0]
        # speed = [R_REF for i in range(len(altRange))]
        speed = alfSpdInertial

        # create the X dimension
        simXRange = np.linspace(0, EToggles.lambdaPerp0, EToggles.lambdaPerp_Rez)

        # create a meshgrid
        for tme, timeVal in enumerate(timeRange):

            spaceGrid = np.zeros(shape=(len(altRange),len(simXRange)))

            for z,x in product(*[range(len(altRange)),range(len(simXRange))]):
                zVal = altRange[z]
                xVal = simXRange[x]
                spaceGrid[z][x] = Eperp_generator(x=xVal,
                                                  z=zVal,
                                                  t=timeVal,
                                                  Vel=speed[z],
                                                  tau=EToggles.tau0,
                                                  KperpVal=kperp[z])
            Eperp.append(spaceGrid)

        # plotting
        if plotBool:
            from ACESII_code.class_var_func import prgMsg,Done
            import time
            start_time = time.time()
            prgMsg('Creating Eperp Plots')
            import matplotlib.pyplot as plt
            plt.rcParams.update({'font.size': 22})
            for timeIndexChoice in range(len(timeRange)):
                fig,ax = plt.subplots()
                fig.set_figwidth(15)
                fig.set_figheight(15)
                ax.set_title('$E_{\perp}$ Propogation vs Altitude \n' + f't = {timeRange[timeIndexChoice]}'+r'  $\tau_{0}$=' +f'{EToggles.tau0} s' +'\n' + '$\lambda_{\perp}$ =' + f'{EToggles.lambdaPerp0} [m],  ' + '$\omega_{wave}$ =' + f'{EToggles.waveFreq_Hz} [Hz]')
                cmap = ax.pcolormesh(simXRange/EToggles.lambdaPerp0,altRange/R_REF,Eperp[timeIndexChoice]*m_to_km, cmap='turbo', vmin=0*-1*EToggles.Eperp0*m_to_km,vmax=EToggles.Eperp0*m_to_km)
                ax.set_xlabel('X Distance [$\lambda_{\perp}$]')
                ax.set_ylabel('Z Distance [$R_{E}$]')
                ax.grid(True)
                cbar = plt.colorbar(cmap)
                cbar.set_label('$|E_{\perp}$| [mV/m]')
                plt.savefig(rf'{GenToggles.simFolderPath}\alfvenWave\Eperp\plots\Eperp_t{timeIndexChoice}.png')
                plt.close()
            Done(start_time)

        # create the output variable
        return np.array(Eperp), simXRange

    ################
    # --- OUTPUT ---
    ################
    if outputData:
        prgMsg('Writing out Eperp Data')

        # get all the variables
        Eperp, simXRange = EperpProfile(altRange=GenToggles.simAlt, timeRange=GenToggles.simTime, **kwargs)

        # --- Construct the Data Dict ---
        exampleVar = {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808,
                      'FORMAT': 'I5', 'UNITS': 'm', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                      'SCALETYP': 'linear', 'LABLAXIS': 'simAlt'}


        data_dict = {'Eperp': [Eperp, {'DEPEND_0': 'simTime', 'DEPEND_1': 'simAlt', 'DEPEND_2': 'simXRange', 'UNITS': 'V/m', 'LABLAXIS': 'Eperp'}],
                     'simXRange': [simXRange, {'DEPEND_0': 'simAlt', 'UNITS': 'T', 'LABLAXIS': 'simXRange'}],
                     'simTime': [GenToggles.simTime, {'DEPEND_0': 'simAlt', 'UNITS': 'm', 'LABLAXIS': 'simTime'}],
                     'simAlt': [GenToggles.simAlt, {'DEPEND_0': 'simAlt', 'UNITS': 'm', 'LABLAXIS': 'simAlt'}]}

        # update the data dict attrs
        for key, val in data_dict.items():
            newAttrs = deepcopy(exampleVar)

            for subKey, subVal in data_dict[key][1].items():
                newAttrs[subKey] = subVal

            data_dict[key][1] = newAttrs

        # --- output the data ---
        outputPath = rf'{GenToggles.simOutputPath}\Eperp\Eperp.cdf'
        outputCDFdata(outputPath, data_dict)
        Done(start_time)



#################
# --- EXECUTE ---
#################
alfvenEperpGenerator(outputData=outputData,showPlot=plot_Eperp)