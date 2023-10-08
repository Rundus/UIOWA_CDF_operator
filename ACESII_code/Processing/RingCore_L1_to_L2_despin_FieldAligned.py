# # --- RingCore_L1_to_L2_despin_FieldAligned.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: Just interpolate the attitude solution data and apply the wallops DCM to the magnetometer data, getting it in ENU coordinates.
# Apply it over the entire journey of the flight



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_attitude = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier_despin = 'l2' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# --- --- --- Reduce Data --- --- ---
reduceTimes = [ [dt.datetime(2022, 11, 20, 17, 22, 00, 000000), dt.datetime(2022, 11, 20, 17, 30, 00, 000000)],
                [dt.datetime(2022, 11, 20, 17, 24, 00, 000000), dt.datetime(2022, 11, 20, 17, 28, 00, 000000)]]

# --- --- --- MODEL COMPARE --- --- ---
plotMODELcompare = False
wModel = 1 # 0 --> IGRF, 1--> CHAOS7
outputData = False

# --- --- --- Subtract Model --- --- ---
SECTION_subtractCHAOSModel = False
plotSubtractedModel = False
outputSubtractedModel = False

# --- --- --- Field Aligned --- --- ---
SECTION_convertToFieldAligned = True
outputFieldAligned = True


fitResults = {
            'Bx': {'Spin Amp': 25.42873940404161, 'Spin Freq': 0.6463295881639182, 'Spin Phase': 91.9759995936283,
                   'Cone Amp': 625.8772357084948, 'Cone Freq': 0.05294818121871208, 'Cone Phase': -138.77308595997619,
                   'Offset': -44919.748937299344},
            'By': {'Spin Amp': 7.378420193701481, 'Spin Freq': 0.6442248190622027, 'Spin Phase': 109.20255873087793,
                   'Cone Amp': 1380.5616077430786, 'Cone Freq': 0.02700105226961604, 'Cone Phase': 109.87799606103452,
                   'Offset': -139.74554466082876},
            'Bz': {'Spin Amp': 8.095746809541962, 'Spin Freq': 0.6442537451458561, 'Spin Phase': 19.11852573798773,
                   'Cone Amp': 1257.0313161879794, 'Cone Freq': 0.026874206798816504, 'Cone Phase': -69.78175516947503,
                   'Offset': 32.456720919269245}
        }

coneFreq = sum([fitResults['By']['Cone Freq'], fitResults['Bz']['Cone Freq']]) / 2
spinFreq = sum([fitResults['Bz']['Spin Freq'], fitResults['By']['Spin Freq'], fitResults['Bz']['Spin Freq']]) / 3 if wRocket == 4 else 0.55


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from pyIGRF import igrf_value
from ACESII_code.myImports import *
from scipy.interpolate import CubicSpline
from ACESII_code.class_var_func import CHAOS


def RingCore_L1_to_L2_Despin(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L0_ACES_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*RingCore_rktFrm*')
    inputFiles_attitude = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[wflyer]}{modifier}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]

    fileoutName_despin = f'ACESII_{rocketID}_l2_RingCore_DeSpun'


    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'DeSpining RingCore Data' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the Magnetometer file ---
        prgMsg(f'Loading data from {inputPath_modifier} RingCore Files')
        data_dict_mag = loadDictFromFile(inputFiles[wFile],{})
        Done(start_time)

        # --- get the data from the attitude file ---
        prgMsg(f'Loading data from {inputPath_modifier_attitude} Files')
        data_dict_attitude = loadDictFromFile(inputFiles_attitude[0], {})
        Done(start_time)

        ########################
        # --- Reduce dataset ---
        ########################
        prgMsg('Reducing Dataset')

        targetTimes = reduceTimes[wRocket - 4]

        # --- apply reduction to mag data ---
        lowCutoff, highCutoff = np.abs(data_dict_mag['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_mag['Epoch'][0] - targetTimes[1]).argmin()

        for key, val in data_dict_mag.items():
            data_dict_mag[key][0] = np.array(data_dict_mag[key][0][lowCutoff:highCutoff])

        data_dict_mag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag['Epoch'][0]])
        Epoch_seconds = np.array([(tme - data_dict_mag['Epoch'][0][0]) / 1E9 for tme in data_dict_mag['Epoch'][0]])
        Epoch_dt = np.array([pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_mag['Epoch'][0]])

        # --- apply reduction to attitude data ---
        lowCutoff, highCutoff = np.abs(data_dict_attitude['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_attitude['Epoch'][0] - targetTimes[1]).argmin()

        for key, val in data_dict_attitude.items():
            data_dict_attitude[key][0] = np.array(data_dict_attitude[key][0][lowCutoff:highCutoff])

        data_dict_attitude['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude['Epoch'][0]])
        Done(start_time)

        ############################################################
        # --- interpolate attitude data up to magnetometer epoch ---
        ############################################################
        prgMsg('Interpolating Attitude Data')

        # results from RinCore_L1_DetermineAdjustments_Before_Despin that suggest a linear offset of the attitude data in time
        # offsetResults_intercept = [104944444.44444445, 117000000]
        # offsetResults_intercept = [0, 117000000]
        offsetResults_intercept = [0, 0]
        # offsetResults_VectorScalar = np.array([[0, 0, -234.47368421052633],[0,0,-355.47368421052633]])
        offsetResults_VectorScalar = np.array([[0, 0, 0], [0, 0, 0]])

        # Do the interpolation
        Epoch_attitude_loop = np.array([int(tme + offsetResults_intercept[wRocket - 4]) for tme in data_dict_attitude['Epoch'][0]])
        dataKeys = ['Epoch', 'Alt', 'Latgd', 'Long', 'Y_Az', 'a11', 'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33']
        dataKeysVal = [deepcopy(data_dict_mag['Epoch'][0]), [], [], [], [], [], [], [], [], [], [], [], [], []]
        attitudeData = [deepcopy(data_dict_attitude[key][0]) for key in dataKeys]  # a list to contain the attitude only the data that I care about
        dataInterp_dict_attitude = {key: value for key, value in zip(dataKeys, dataKeysVal)}

        counter = 0
        for key, newDataList in dataInterp_dict_attitude.items():
            if key != 'Epoch':
                # --- cubic interpolation ---
                splCub = CubicSpline(Epoch_attitude_loop, attitudeData[counter])

                # evaluate the interpolation at all the epoch_mag points
                dataInterp_dict_attitude[key] = np.array([splCub(timeVal) for timeVal in dataInterp_dict_attitude['Epoch']])
            counter += 1

        # convert altitude to km
        dataInterp_dict_attitude['Alt'] = dataInterp_dict_attitude['Alt']/1000

        Done(start_time)

        #######################################
        # --- REVERSE ROTATE TO REMOVE SPIN ---
        #######################################

        prgMsg('Applying Despin')

        # define the Yaw,Pitch,Roll angles to use over the timeseries. USE DEGREES since DCM does a radian conversion
        B_rkt = np.array([[data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]] for i in range(len(Epoch_seconds))])

        DCMmat = np.array([
            [[dataInterp_dict_attitude['a11'][i], dataInterp_dict_attitude['a12'][i], dataInterp_dict_attitude['a13'][i]],
             [dataInterp_dict_attitude['a21'][i], dataInterp_dict_attitude['a22'][i], dataInterp_dict_attitude['a23'][i]],
             [dataInterp_dict_attitude['a31'][i], dataInterp_dict_attitude['a32'][i], dataInterp_dict_attitude['a33'][i]]]
            for i in range(len(Epoch_seconds))
        ])

        # prepare data for further processing
        data_for_output = np.array([np.matmul(DCMmat[i], B_rkt[i]) for i in range(len(Epoch_seconds))]) + offsetResults_VectorScalar[wRocket - 4]

        Done(start_time)

        ##########################################
        # --- COMPARE DESPIN TO IGRF OR CHAOS7 ---
        ##########################################

        if plotMODELcompare:

            modelNames = ['IGRF', 'CHAOS7']


            if wModel == 0:
                prgMsg('Comparing to IGRF Model')

                # --- get IGRF ENU ---
                date = 2022 + 323 / 365  # Corresponds to 11/20/2022
                ### IGRF info ###
                # [3] North Comp (+ N | - S)
                # [4] East Comp (+ E | - W)
                # [5] Vertical Comp (+ D | - U)
                # [6] Total Field

                IGRF = np.array([igrf_value(dataInterp_dict_attitude['Latgd'][i], dataInterp_dict_attitude['Long'][i], dataInterp_dict_attitude['Alt'][i] / 1000, date) for i in range(len(dataInterp_dict_attitude['Epoch']))])
                B_model = np.array([[vec[4], vec[3], -1 * vec[5]] for vec in IGRF])

            elif wModel:
                prgMsg('Loading CHAOS model')
                Epoch_in_dt = np.array([pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_mag['Epoch'][0]])
                B_model = CHAOS(lat=dataInterp_dict_attitude['Latgd'],
                                long=dataInterp_dict_attitude['Long'],
                                alt=dataInterp_dict_attitude['Alt'],
                                times=Epoch_in_dt)

            Done(start_time)


            # compare ENU IGRF to my unspun rocket data
            fig, ax = plt.subplots(3)
            fig.suptitle(f'ACESII {rocketID} {modelNames[wModel]} vs Attitude DCM comparison')

            mapping = [0, 1, 2]

            # East (B_rkt_Y)
            B_rkt_east_plot, = ax[0].plot(Epoch_dt, data_for_output[:, mapping[0]],label='B_rkt_East')
            ax[0].plot(Epoch_dt, B_model[:, 0],label='Model East')
            ax[0].set_ylabel(f'B_East [nT]')
            ax[0].legend()
            # ax[0].set_ylim(-10000, 10000)

            # North (B_rkt_Z)
            i = 1
            B_rkt_north_plot, = ax[1].plot(Epoch_dt, data_for_output[:, mapping[1]],label='B_rkt_North')
            ax[1].plot(Epoch_dt, B_model[:, 1], label='Model North')
            ax[1].set_ylabel(f'B_North [nT]')
            ax[1].legend()


            # Up (B_rkt_X)
            B_rkt_up_plot, = ax[2].plot(Epoch_dt, data_for_output[:, mapping[2]],label='B_rkt_Up')
            ax[2].plot(Epoch_dt, B_model[:, 2], label='Model Up')
            ax[2].set_ylabel(f'B_Up [nT]')
            ax[2].set_xlabel('Epoch')
            ax[2].legend()

            plt.show()

            Done(start_time)

        if outputData:
            prgMsg('Creating despin output file')

            # create the output data_dict
            data_dict = deepcopy(data_dict_mag)
            comps = ['Bx', 'By', 'Bz', 'Bmag']
            newComps = ['B_East', 'B_North', 'B_Up', 'B_Mag']
            data_for_output_despin = np.array([[data_for_output[i][0],data_for_output[i][1],data_for_output[i][2], np.linalg.norm(data_for_output[i])] for i in range(len(data_for_output))])

            # --- Magnetic Components ---
            # get the attributes of the old components and replace them
            for i, key in enumerate(comps):
                newAttrs = deepcopy(data_dict[key][1])
                newAttrs['LABLAXIS'] = newComps[i]

                # remove the old key
                del data_dict[key]

                # append the new key
                data_dict = {**data_dict, **{newComps[i]: [data_for_output_despin[:, i], newAttrs]}}

            outputPath = f'{rocketFolderPath}{outputPath_modifier_despin}\{fliers[wflyer]}\\{fileoutName_despin}_ENU.cdf'

            outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, 'RingCore')

            Done(start_time)

        if SECTION_subtractCHAOSModel:


            prgMsg('Loading CHAOS model')
            Epoch_in_dt = np.array([pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_mag['Epoch'][0]])
            B_model = CHAOS(lat=dataInterp_dict_attitude['Latgd'],
                            long=dataInterp_dict_attitude['Long'],
                            alt=dataInterp_dict_attitude['Alt'],
                            times=Epoch_in_dt)
            Done(start_time)

            prgMsg('Subtracting Model')
            data_for_output = data_for_output - B_model
            Done(start_time)

            if plotSubtractedModel:

                prgMsg('Plotting Subtraction')

                fig, ax = plt.subplots(3)
                compNames = ['B_East', 'B_North', 'B_Up']

                for i in range(3):
                    ax[i].plot(Epoch_dt, data_for_output[:, i])
                    ax[i].set_ylabel(compNames[i])

                ax[2].set_xlabel('Epoch')
                plt.show()

                Done(start_time)

            if outputSubtractedModel:
                prgMsg('Creating model subtracted output file')

                # create the output data_dict
                data_dict = deepcopy(data_dict_mag)
                comps = ['Bx', 'By', 'Bz', 'Bmag']
                newComps = ['B_East', 'B_North', 'B_Up', 'B_Mag']
                data_for_output_despin = np.array([[data_for_output[i][0], data_for_output[i][1], data_for_output[i][2], np.linalg.norm(data_for_output[i])] for i in range(len(data_for_output))])

                # --- Magnetic Components ---
                # get the attributes of the old components and replace them
                for i, key in enumerate(comps):
                    newAttrs = deepcopy(data_dict[key][1])
                    newAttrs['LABLAXIS'] = newComps[i]

                    # remove the old key
                    del data_dict[key]

                    # append the new key
                    data_dict = {**data_dict, **{newComps[i]: [data_for_output_despin[:, i], newAttrs]}}

                outputPath = f'{rocketFolderPath}{outputPath_modifier_despin}\{fliers[wflyer]}\\{fileoutName_despin}_subModel.cdf'

                outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, 'RingCore')

                Done(start_time)

        if SECTION_convertToFieldAligned:

            # #############################################################
            # --- Convert ENU coordinates to Field Aligned coordinates ---
            # #############################################################
            prgMsg('Converting to Field Aligned Coordinates')

            # Get the Data
            B_rkt_ENU = np.array([np.matmul(DCMmat[i], B_rkt[i]) for i in range(len(Epoch_seconds))]) + offsetResults_VectorScalar[wRocket - 4]
            Epoch_in_dt = np.array([pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_mag['Epoch'][0]])
            B_model = CHAOS(lat=dataInterp_dict_attitude['Latgd'],
                            long=dataInterp_dict_attitude['Long'],
                            alt=dataInterp_dict_attitude['Alt'],
                            times=Epoch_in_dt) # CHAOS in ENU coordinates


            # --- Convert B-Data to GEO (ECEF) XYZ coordinates ---
            from ACESII_code.class_var_func import ENUtoECEF
            ENUtoGEOmatrix = np.array([ENUtoECEF(Lat=dataInterp_dict_attitude['Latgd'][i], Long=dataInterp_dict_attitude['Long'][i]) for i in range(len(Epoch_in_dt))])


            B_rkt_GEO = np.array([np.matmul(ENUtoGEOmatrix[i], B_rkt_ENU[i]) for i in range(len(Epoch_in_dt))])
            B_CHAOS_GEO = np.array([np.matmul(ENUtoGEOmatrix[i], B_model[i]) for i in range(len(Epoch_in_dt))])

            # --- determine the Payload's Position Vector in GEO (ECEF) coordinate XYZ ---
            R_REF = 6371.2 # earth Radius in km
            Radius = dataInterp_dict_attitude['Alt'] + R_REF
            coLatRad = [np.radians(90-lat) for lat in dataInterp_dict_attitude['Latgd']]
            LongRad = [np.radians(long) for long in dataInterp_dict_attitude['Long']]
            Rsc = np.array([
                [Radius[i]*np.sin(coLatRad[i])*np.cos(LongRad[i]),
                 Radius[i]*np.sin(coLatRad[i])*np.sin(LongRad[i]),
                 Radius[i]*np.cos(coLatRad[i])] for i in range(len(Epoch_in_dt))])


            # --- calculate Field Aligned unit vectors over the duration of the flight ---

            # pHat comes from the CHAOS model direction of B in GEO
            pHat = np.array([B_CHAOS_GEO[i]/np.linalg.norm(B_CHAOS_GEO[i]) for i in range(len(Epoch_in_dt))])

            # e-hat comes from the cross of pHat and the Rocket's radius vector (in geomagnetic coordinates)
            eHat = np.array([np.cross(pHat[i], Rsc[i])/np.linalg.norm(np.cross(pHat[i], Rsc[i])) for i in range(len(Epoch_in_dt))])

            # rHat comes from the cross of eHat and pHat
            rHat = np.array([np.cross(eHat[i], pHat[i]) for i in range(len(Epoch_in_dt))])

            # form the transformation matrix FROM GEO TO FIELD ALIGNED
            DCM_FA = np.array([[eHat[i], pHat[i], rHat[i]] for i in range(len(Epoch_in_dt))])

            # --- Transform B_ENU to B_EPR ---
            B_rkt_EPR = np.array([np.matmul(DCM_FA[i], B_rkt_GEO[i]) for i in range(len(Epoch_in_dt))])

            data_for_output = deepcopy(B_rkt_EPR)
            Done(start_time)

            if outputFieldAligned:
                prgMsg('Creating FieldAligned output file')

                # create the output data_dict
                data_dict = deepcopy(data_dict_mag)
                comps = ['Bx', 'By', 'Bz', 'Bmag']
                newComps = ['B_e', 'B_p', 'B_r', 'B_Mag']
                data_for_output_despin = np.array([[data_for_output[i][0], data_for_output[i][1], data_for_output[i][2], np.linalg.norm(data_for_output[i])] for i in range(len(data_for_output))])

                # --- Magnetic Components ---
                # get the attributes of the old components and replace them
                for i, key in enumerate(comps):
                    newAttrs = deepcopy(data_dict[key][1])
                    newAttrs['LABLAXIS'] = newComps[i]

                    # remove the old key
                    del data_dict[key]

                    # append the new key
                    data_dict = {**data_dict, **{newComps[i]: [data_for_output_despin[:, i], newAttrs]}}

                outputPath = f'{rocketFolderPath}{outputPath_modifier_despin}\{fliers[wflyer]}\\{fileoutName_despin}_FieldAligned.cdf'

                outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, 'RingCore')

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

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        RingCore_L1_to_L2_Despin(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            RingCore_L1_to_L2_Despin(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            RingCore_L1_to_L2_Despin(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)