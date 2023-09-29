# --- PoyntingFlux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Determine the PoyntingFLux of the data using E-Field and B-Field Measurements.
# For the low flyer, it ONLY accepts despun data, high flyer has its own case



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from ACESII_code.myImports import *

start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False
printMagFiles = True
printElecFiles = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

modifier = ''
inputPath_modifier_elec = 'l1'
inputPath_modifier_mag = 'l2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science/PoyntingFlux' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# reduce the size of the data
reduceData = True
targetTimes = [pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 24, 26, 000)),
               pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 25, 15, 000))]

plotSPoynting = False
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.interpolate import CubicSpline
from ACESII_code.class_var_func import butter_filter, u0, IonMasses

def PoyntingFlux(wRocket, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles_elec = glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wflyer]}{modifier}\*E_Field*')
    inputFiles_mag = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}{modifier}\*RingCore*')

    fileoutName = f'ACESII_{rocketID}_PoyntingFlux.cdf'

    if justPrintFileNames:
        if printMagFiles:
            for i, file in enumerate(inputFiles_mag):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_mag[i], round(getsize(file) / (10 ** 6), 1)))
        elif printElecFiles:
            for i, file in enumerate(inputFiles_elec):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_elec[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Calculating Poynting flux for ACESII {rocketID}' + color.END)

        # --- get the data from the mag file ---
        prgMsg(f'Loading data from mag Files')
        data_dict_mag = loadDictFromFile(inputFiles_mag[0], {})
        data_dict_mag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag['Epoch'][0]])
        indicies = [np.abs(data_dict_mag['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_mag['Epoch'][0] - targetTimes[1]).argmin()]

        for key, val in data_dict_mag.items():
            data_dict_mag[key][0] = data_dict_mag[key][0][indicies[0]:indicies[1]]

        # create vector variable and convert to tesla
        dB = np.array([np.array([data_dict_mag['dB_east'][0][i], data_dict_mag['dB_north'][0][i], data_dict_mag['dB_up'][0][i]])*(1E-9) for i in range(len(data_dict_mag['Epoch'][0]))])

        Done(start_time)

        if wRocket == 4:

            # collect the Magnitude of B from L1 spun data
            prgMsg('Getting Bmag')
            inputFileBmag = glob('C:\Data\ACESII\L1\high\*RingCore_rktFrm*')
            data_dict_Bmag = loadDictFromFile(inputFileBmag[0], {})
            data_dict_Bmag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_Bmag['Epoch'][0]])
            Done(start_time)

            prgMsg('Getting Plasma Density')
            inputFileBmag = glob('C:\Data\ACESII\science\Langmuir\high\*Temp&Density*')
            data_dict_density = loadDictFromFile(inputFileBmag[0], {})
            data_dict_density['fixed_Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_density['fixed_Epoch'][0]])
            Done(start_time)

            # reduce the datasets:
            prgMsg('Reducing Data')

            # Bmag
            indicies = [np.abs(data_dict_Bmag['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_Bmag['Epoch'][0] - targetTimes[1]).argmin()]
            for key, val in data_dict_Bmag.items():
                data_dict_Bmag[key][0] = deepcopy(data_dict_Bmag[key][0][indicies[0]:indicies[1]])

            # Density
            indicies = [np.abs(data_dict_density['fixed_Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_density['fixed_Epoch'][0] - targetTimes[1]).argmin()]
            for key, val in data_dict_density.items():
                data_dict_density[key][0] = deepcopy(data_dict_density[key][0][indicies[0]:indicies[1]])

            Done(start_time)

            #################################
            # --- CALCULATE POYNTING FLUX ---
            #################################
            prgMsg('Calculating Poynting Flux using EigenFunction')

            # down sample the density data onto the magnetometer data
            indiciesDownsampled = [np.abs(data_dict_density['fixed_Epoch'][0] - data_dict_mag['Epoch'][0][i]).argmin() for i in range(len(data_dict_mag['Epoch'][0]))]
            plasmaDensity = np.array([data_dict_density['fixed_ni_density'][0][index] for index in indiciesDownsampled])

            # calculate Alfven Velocity
            AlfvenVelocity = np.array([((1E-9)*data_dict_Bmag['Bmag'][0][i])/np.sqrt(u0*plasmaDensity[i]*IonMasses[0]) for i in range(len(plasmaDensity))])


            # calculate Alfven eigenfunction E


            # calculate Poyning Flux
            S = np.array([(AlfvenVelocity[i]/(2*u0))*np.array([dB[i][0]**2, dB[i][1]**2, dB[i][2]**2]) for i in range(len(data_dict_mag['Epoch'][0])) ])

            Done(start_time)

            if plotSPoynting:
                prgMsg('Plotting HF Poynting Flux')
                Epoch = [ pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_mag['Epoch'][0]]
                fig, ax = plt.subplots(3)
                fig.suptitle('Alfven EigenFunction')
                compNames = ['$S_{East}$', '$S_{North}$', '$S_{up}$']
                for i in range(3):
                    ax[i].plot(Epoch, S[:, i])
                    ax[i].set_ylabel(compNames[i] + '[$W/m^{2}$]')
                    ax[i].set_ylim(-0.01, 0.15)

                plt.show()
                Done(start_time)

        elif wRocket == 5:

            # --- get the data from the electric file ---
            prgMsg(f'Loading data from Electric Field Files')
            data_dict_elec = loadDictFromFile(inputFiles_elec[0], {})
            data_dict_elec['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_elec['Epoch'][0]])

            # --- reduce each dataset ---
            indicies = [np.abs(data_dict_elec['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_elec['Epoch'][0] - targetTimes[1]).argmin()]
            for key, val in data_dict_elec.items():
                data_dict_elec[key][0] = deepcopy(data_dict_elec[key][0][indicies[0]:indicies[1]])

            E_Field = np.array([ [data_dict_elec['E_east'][0][i], data_dict_elec['E_north'][0][i],data_dict_elec['E_up'][0][i]] for i in range(len(data_dict_elec['Epoch'][0]))])
            Done(start_time)

            ####################################################
            # --- INTERPOLATE MAG/ATTITUDE ONTO EFI TIMEBASE ---
            ####################################################

            prgMsg('Interpolating Mag Data')
            spline_dict_mag = {}
            dataInterp_dict_mag = {key: [] for key, val in data_dict_mag.items()}
            for key, newDataList in data_dict_mag.items():
                if 'Epoch'.lower() not in key.lower():
                    # --- cubic interpolation ---
                    splCub = CubicSpline(data_dict_mag['Epoch'][0], data_dict_mag[key][0])

                    # store the spline information for later
                    spline_dict_mag = {**spline_dict_mag, **{key: splCub}}

                    # evaluate the interpolation at all the epoch_mag points
                    dataInterp_dict_mag[key] = np.array([splCub(timeVal) for timeVal in data_dict_elec['Epoch'][0]])

            dataInterp_dict_mag['Epoch'] = np.array(data_dict_elec['Epoch'][0])

            # create a magnetic vector variable
            B_interp = np.array([ [dataInterp_dict_mag['B_east'][i],dataInterp_dict_mag['B_north'][i],dataInterp_dict_mag['B_up'][i]] for i in range(len(dataInterp_dict_mag['Epoch']))])

            Done(start_time)

            #################################
            # --- CALCULATE POYNTING FLUX ---
            #################################
            prgMsg('Calculating Poynting Flux')
            S = np.array([ u0*(np.cross(E_Field[i], B_interp[i])) for i in range(len(E_Field))])
            if plotSPoynting:
                Epoch = data_dict_elec['Epoch'][0]
                fig, ax = plt.subplots(3)
                fig.suptitle('B_filtered')
                ax[0].plot(Epoch, S[:, 0])
                ax[0].set_ylabel('$S_{East}$')
                ax[1].plot(Epoch, S[:, 1])
                ax[1].set_ylabel('$S_{North}$')
                ax[2].plot(Epoch, S[:, 2])
                ax[2].set_ylabel('$S_{Up}$')
                plt.show()


        prgMsg('Preparing Data')
        # --- prepare data for output ---
        if wRocket == 4:
            data_for_output = np.array([[S[i][0], S[i][1], S[i][2], np.linalg.norm(S[i])] for i in range(len(dB))])
        elif wRocket == 5:
            data_for_output = np.array([[S[i][0],S[i][1],S[i][2],np.linalg.norm(S[i])] for i in range(len(E_Field))])

        data_dict = deepcopy(data_dict_mag)
        comps = ['dB_east', 'dB_north', 'dB_up', 'dBmag']
        newComps = ['S_east', 'S_north', 'S_up', 'S_tot']

        # --- Magnetic Components ---
        # get the attributes of the old components and replace them
        for i, key in enumerate(comps):
            newAttrs = deepcopy(data_dict[key][1])
            newAttrs['LABLAXIS'] = newComps[i]

            # remove the old key
            del data_dict[key]

            # append the new key
            data_dict = {**data_dict, **{newComps[i]: [data_for_output[:, i], newAttrs]}}

        if wRocket == 5:
            data_dict['Epoch'][0] = data_dict_elec['Epoch'][0]

        Done(start_time)



        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('Creating output file')

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, 'PoyntingFlux')

            Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no electric field .cdf files in the specified directory' + color.END)
elif len(glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no B-field .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        PoyntingFlux(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        PoyntingFlux(wRocket, rocketFolderPath, justPrintFileNames, wflyer)
