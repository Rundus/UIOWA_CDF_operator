# --- alfven_Eperp_Generator.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:
# Generates the entire E-Field for all time in the simulation and returns a variable that looks like:
# [
#   [Ez(x=0,t=0),Ez(x=1,t=0),Ez(x=2,t=0)...., Ez(x=len(Alt),t=0)],
#   [Ez(x=0,t=1),Ez(x=1,t=1),Ez(x=2,t=1)...., Ez(x=len(Alt),t=1)]
#   ,...]

# --- imports ---
from ACESII_code.myImports import *
from ACESII_code.class_var_func import L2_ACES_Quick
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.simToggles import m_to_km, R_REF, GenToggles,EToggles
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.plasmaEnvironment.plasmaEnvironment_Generator import generatePlasmaEnvironment
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.geomagneticField.geomagneticField_Generator import generateGeomagneticField
from ACESII_code.class_var_func import loadDictFromFile
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

            for z,x in itertools.product(*[range(len(altRange)),range(len(simXRange))]):
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

        # --- ACES II Flight/Integration Data ---
        wRocket = 4
        rocketAttrs, b, c = ACES_mission_dicts()
        globalAttrsMod = rocketAttrs.globalAttributes[wRocket - 4]
        globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
        ModelData = L2_ACES_Quick(wRocket - 4)

        # get all the variables
        Eperp, simXRange = EperpProfile(altRange=GenToggles.simAlt, timeRange=GenToggles.simTime, **kwargs)

        # --- Construct the Data Dict ---
        exampleVar = {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                      'FORMAT': 'I5', 'UNITS': 'm', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                      'SCALETYP': 'linear', 'LABLAXIS': 'simAlt'}

        data_dict = {'Eperp': [Eperp, deepcopy(exampleVar)],
                     'simXRange': [simXRange, deepcopy(exampleVar)],
                     'simTime': [GenToggles.simTime, deepcopy(exampleVar)],
                     'simAlt': [GenToggles.simAlt, deepcopy(exampleVar)]}

        data_dict['Eperp'][1]['UNITS'] = 'V/m'
        data_dict['Eperp'][1]['LABLAXIS'] = 'Eperp'
        data_dict['Eperp'][1]['DEPEND_0'] = 'simTime'
        data_dict['Eperp'][1]['DEPEND_1'] = 'simAlt'
        data_dict['Eperp'][1]['DEPEND_2'] = 'simXRange'

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
        outputPath = rf'{GenToggles.simFolderPath}\alfvenWave\Eperp\alfvenWave_Eperp.cdf'
        outputCDFdata(outputPath, data_dict, ModelData, globalAttrsMod, 'simulation')
        Done(start_time)



#################
# --- EXECUTE ---
#################
alfvenEperpGenerator(outputData=outputData,showPlot=plot_Eperp)