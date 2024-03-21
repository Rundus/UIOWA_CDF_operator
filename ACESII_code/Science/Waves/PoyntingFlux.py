# --- PoyntingFlux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Determine the PoyntingFLux of the data using E-Field and B-Field Measurements.
# For the low flyer, it ONLY accepts despun data, high flyer has its own case


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from ACESII_code.myImports import *

start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False
printMagFiles = True
printElecFiles = True

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

modifier = ''
inputPath_modifier_elec = 'L3\deltaE'
wMagFile = 0
inputPath_modifier_mag = 'L3\deltaB' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
wEFIFile = 0
outputPath_modifier = 'science/PoyntingFlux' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# --- --- --- Which Data --- -- ---
useDelta_E_B = True # use the deltaB, deltaE data
targetVar = [[dt.datetime(2022,11,20,17,20,00,00),dt.datetime(2022,11,20,17,28,00,00)],'Epoch']
# --- --- --- PLOT --- --- ---
plotSPoynting = True
# --- --- --- OUTPUT --- --- ---
outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---

from ACESII_code.class_var_func import u0, IonMasses, InterpolateDataDict, dateTimetoTT2000

def PoyntingFlux(wRocket, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]

    inputFiles_elec = glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wflyer]}{modifier}\*Field_Aligned*')
    inputFiles_mag = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}{modifier}\*Field_Aligned*')

    # determine which coordinates system its in based on magnetic field data
    FileName = inputFiles_mag[wMagFile].replace('.cdf', '')
    wCoordinate = 'Field_Aligned'
    perpCoordinates = [0, 2]
    parCoordinates = [1]

    fileoutName = f'ACESII_{rocketID}_PoyntingFlux_deltaS_{wCoordinate}.cdf' if useDelta_E_B else f'ACESII_{rocketID}_PoyntingFlux_flight_{wCoordinate}.cdf'

    if justPrintFileNames:
        if printMagFiles:
            print('--- B-FIELD FILES ---')
            for i, file in enumerate(inputFiles_mag):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_mag[i], round(getsize(file) / (10 ** 6), 1)))
            print('\n')

        if printElecFiles:
            print('--- E-FIELD FILES ---')
            for i, file in enumerate(inputFiles_elec):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_elec[i], round(getsize(file) / (10 ** 6), 1)))
            print('\n')

        return
    else:
        print('\n')
        print(color.UNDERLINE + f'Calculating Poynting flux for ACESII {rocketID}' + color.END)

        # --- get the data from the mag file ---
        prgMsg(f'Loading data from mag Files')
        data_dict_mag = loadDictFromFile(inputFiles_mag[wMagFile],targetVar=targetVar)

        # component names for the magnetic field
        compNames_mag = [key for key, val in data_dict_mag.items() if key.lower() not in ['epoch','db_mag','alt', 'lat', 'long', 'alt_geom', 'lat_geom', 'long_geom']]

        # component names for the poynting flux
        compNamesS = [name.replace('B', 'S') for name in compNames_mag]

        # create vector variable and convert to tesla
        Done(start_time)

        if wRocket == 4:

            # collect the Magnitude of B from L1 spun data
            prgMsg('Getting Bmag')
            inputFileBmag = glob('C:\Data\ACESII\L1\high\*RingCore_rktFrm*')
            data_dict_Bmag = loadDictFromFile(inputFileBmag[0], {}, reduceData, targetTimes=targetTimes)
            badIndicies = [i for i, val in enumerate(data_dict_Bmag['Bmag'][0]) if not np.isfinite(val)] # quality assurance step to remove and NAN values
            data_dict_Bmag['Epoch'][0] = np.delete(data_dict_Bmag['Epoch'][0],badIndicies)
            data_dict_Bmag['Bmag'][0] = np.delete(data_dict_Bmag['Bmag'][0], badIndicies)
            Done(start_time)

            prgMsg('Getting Plasma Density')
            inputFileDensity = glob('C:\Data\ACESII\science\Langmuir\high\*fixed*')
            data_dict_density = loadDictFromFile(inputFileDensity[0], {}, reduceData, targetTimes=targetTimes)
            Done(start_time)

            #########################
            # --- DownSample Data ---
            #########################
            prgMsg('DownSampling Density Data')


            data_dict_BmagInterp = InterpolateDataDict(InputDataDict= data_dict_Bmag,
                                                       InputEpochArray= dateTimetoTT2000(data_dict_Bmag['Epoch'][0], False),
                                                       wKeys=['Bmag'],
                                                       targetEpochArray=dateTimetoTT2000(data_dict_mag['Epoch'][0], False))

            data_dict_densityInterp = InterpolateDataDict(InputDataDict=data_dict_density,
                                                       InputEpochArray=dateTimetoTT2000(data_dict_density['Epoch'][0], False),
                                                       wKeys=['ni', 'Epoch'],
                                                       targetEpochArray=dateTimetoTT2000(data_dict_mag['Epoch'][0], False))

            # --- prepare some variables ---
            plasmaDensity = (100**3)*data_dict_densityInterp['ni'][0]
            B_perp = np.array([
                [(1E-9)*data_dict_mag[compNames_mag[perpCoordinates[0]]][0][i],(1E-9)*data_dict_mag[compNames_mag[perpCoordinates[1]]][0][i]]
                for i in range(len(data_dict_mag['Epoch'][0]))
            ])

            Done(start_time)

            #################################
            # --- CALCULATE POYNTING FLUX ---
            #################################
            prgMsg('Calculating Poynting Flux using EigenFunction')

            # calculate Alfven Velocity
            AlfvenVelocity = np.array([((1E-9)*data_dict_BmagInterp['Bmag'][0][i])/np.sqrt(u0*plasmaDensity[i]*IonMasses[0]) for i in range(len(plasmaDensity))])

            # calculate TOTAL Poyning Flux
            S_tot = np.array([(AlfvenVelocity[i]/(2*u0))*np.dot(B_perp[i],B_perp[i]) for i in range(len(data_dict_mag['Epoch'][0]))])

            # construct the Poynting Flux vector
            S = []
            for i in range(len(S_tot)):
                temp = [0,0,0]
                temp[perpCoordinates[0]] = 0
                temp[perpCoordinates[1]] = 0
                temp[parCoordinates[0]] = S_tot[i]
                S.append(temp)

            S = np.array(S)

            Done(start_time)

            if plotSPoynting:
                prgMsg('Plotting HF Poynting Flux')
                Epoch = data_dict_mag['Epoch'][0]
                fig, ax = plt.subplots(3)
                fig.suptitle('Alfven EigenFunction')
                for i in range(3):
                    ax[i].plot(Epoch, S[:, i])
                    ax[i].set_ylabel(compNamesS[i] + '[$W/m^{2}$]')


                plt.show()
                Done(start_time)

        elif wRocket == 5:

            # --- get the data from the electric file ---
            prgMsg(f'Loading data from Electric Field Files')
            data_dict_elec = loadDictFromFile(inputFiles_elec[wEFIFile],targetVar=targetVar)
            data_dict_elec['Epoch'][0] = np.array([int(pycdf.lib.datetime_to_tt2000(tme)+ (0.1157*1E9)) for tme in data_dict_elec['Epoch'][0]])
            compNames_elec = [key for key, val in data_dict_elec.items() if key.lower() not in ['epoch', 'de_mag']]

            Done(start_time)

            ###########################################
            # --- INTERPOLATE EFI ONTO MAG TIMEBASE ---
            ###########################################

            prgMsg('Downsampling via Interpolate the EFI Data')
            data_dict_elec_interp = InterpolateDataDict(InputDataDict=data_dict_elec,
                                                        InputEpochArray=data_dict_elec['Epoch'][0],
                                                        wKeys=[],
                                                        targetEpochArray=data_dict_mag['Epoch'][0])

            # get the electric field and convert it to V/m, get the magnetic field and convert it to T
            B_Field = (1E-9) * np.array([np.array([data_dict_mag[compNames_mag[0]][0][i], data_dict_mag[compNames_mag[1]][0][i], data_dict_mag[compNames_mag[2]][0][i]]) for i in range(len(data_dict_mag['Epoch'][0]))])
            E_Field = (1/1000)*np.array([[data_dict_elec_interp[compNames_elec[0]][0][i], data_dict_elec_interp[compNames_elec[1]][0][i], data_dict_elec_interp[compNames_elec[2]][0][i]] for i in range(len(data_dict_mag['Epoch'][0]))])
            Done(start_time)

            #################################
            # --- CALCULATE POYNTING FLUX ---
            #################################
            prgMsg('Calculating Poynting Flux')

            S = np.array([(1/u0)*(np.cross(E_Field[i], B_Field[i])) for i in range(len(data_dict_mag['Epoch'][0]))])


            Done(start_time)
            if plotSPoynting:
                Epoch = data_dict_mag['Epoch'][0]
                fig, ax = plt.subplots(3)
                fig.suptitle('B_filtered')
                ax[0].plot(Epoch, S[:, 0])
                ax[0].set_ylabel(compNamesS[0])
                ax[1].plot(Epoch, S[:, 1])
                ax[1].set_ylabel(compNamesS[1])
                ax[2].plot(Epoch, S[:, 2])
                ax[2].set_ylabel(compNamesS[2])
                for i in range(3):
                    ax[i].set_ylim(-3E-5,3E-5)
                plt.show()


        # --- prepare data for output ---
        prgMsg('Preparing Data')
        data_for_output = np.array([[S[i][0], S[i][1], S[i][2], np.linalg.norm(S[i])] for i in range(len(S))])
        newComps = compNamesS + ['dS_tot'] if useDelta_E_B else compNamesS + ['S_tot']
        print(newComps)
        Done(start_time)

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('creating output file')

            data_dict_output = {}

            for i in range(len(newComps)):
                data = data_for_output[:, i]
                varAttrs = {'LABLAXIS': newComps[i], 'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None,
                            'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'W/m!A2!N',
                            'VALIDMIN': data.min(), 'VALIDMAX': data.max(),
                            'VAR_TYPE': 'data', 'SCALETYP': 'linear'}

                data_dict_output = {**data_dict_output, **{newComps[i]:[data, varAttrs]}}

            Epoch_output = deepcopy(data_dict_mag['Epoch'])

            Epoch_output[1]['VAR_TYPE'] = 'support_data'

            data_dict_output = {**data_dict_output, **{'Epoch': Epoch_output}}

            # add in the attitude data
            keys = ['Alt', 'Lat', 'Long', 'Alt_geom', 'Lat_geom', 'Long_geom','ILat','ILong']
            for key in keys:
                data_dict_output = {**data_dict_output, **{key:data_dict_mag[key]}}

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict_output,instrNam='PoyntingFlux')

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

if len(glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wflyer]}\*.cdf')) == 0 and wRocket ==5:
    print(color.RED + 'There are no electric field .cdf files in the specified directory' + color.END)
elif len(glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no B-field .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        PoyntingFlux(wRocket, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        PoyntingFlux(wRocket, rocketFolderPath, justPrintFileNames, wflyer)
