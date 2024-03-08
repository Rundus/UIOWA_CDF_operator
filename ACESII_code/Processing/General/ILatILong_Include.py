# --- ILatILong_Include.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: I need to create some ILatILong coordinates that will be
# appended to any datafile I want. The code finds the Epoch variable in the data
# and interpolates the new epoch onto the ILatILong variables, then appends it to the datafile
# the CHAOS model is used to get the direction of the magnetic field



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
justPrintFileNames = False
wRocket = 4
inputPath_modifier = 'L1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
wFiles = [0, 2,3,11]
refAlt = 150 # represents 150 km reference altitude that everything is tied to
# ---------------------------
generateILatILong = False # Calculates and Stores the ILat and ILong variables as a .cdf File
plotILatILong = False
updateCDFfile = True # takes a CDF file, interpolates ILat/ILong and then stores it
outputData = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import CHAOS,lat_to_meter,long_to_meter,meter_to_long, InterpolateDataDict


def ILatILong_Include(wRocket, rocketFolderPath, justPrintFileNames, wFile):

    # --- ACES II Flight/Integration Data ---
    wflyer = wRocket-4
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    ModelData = L2_TRICE_Quick(wflyer)


    # --- GET THE INPUT FILE ---
    inputFilesAttitude = glob(f'{rocketFolderPath}attitude\{fliers[wflyer]}\*.cdf')[0]
    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace(inputPath_modifier.lower()+'_', '') for ifile in input_names]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    if generateILatILong:
        # --- --- --- --- --- --- --- -
        # --- GET THE ATTITUDE DATA ---
        # --- --- --- --- --- --- --- -
        prgMsg(f'Loading Attitude Data')
        data_dict_attitude = loadDictFromFile(inputFilePath=inputFilesAttitude)
        Lat = data_dict_attitude['Lat'][0]
        Long = data_dict_attitude['Long'][0]
        Alt = data_dict_attitude['Alt'][0]/1000
        EpochAttitude = data_dict_attitude['Epoch'][0]
        Done(start_time)

        # --- --- --- --- --- --- -
        # --- CREATE ILAT ILONG ---
        # --- --- --- --- --- --- -
        prgMsg('Calculating ILatILong')
        B_CHAOS = CHAOS(lat=Lat,
                        long=Long,
                        alt=Alt,
                        times=EpochAttitude)

        b = B_CHAOS/np.linalg.norm(B_CHAOS) # normalize the magnetic Field


        Latkm = lat_to_meter*Lat
        Longkm = np.array([long_to_meter(Long[i],Lat[i]) for i in range(len(EpochAttitude))])
        ILat = np.zeros(shape=len(EpochAttitude))
        ILong = np.zeros(shape=len(EpochAttitude))

        for i in range(len(EpochAttitude)):
            t = (refAlt - Alt[i])/b[i][2]
            El = Longkm[i] + b[i][0]*t
            Nl = Latkm[i] + b[i][1]*t

            ILat[i] = Nl/lat_to_meter
            ILong[i] = meter_to_long(long_km=El,lat_km=Nl)

        if plotILatILong:
            fig, ax = plt.subplots()
            ax.plot(Long, Lat, color='blue')
            ax.plot(ILong,ILat, color='red')
            plt.show()
            Done(start_time)

        # --- UPDATE ATTITUDE FILE ---
        prgMsg('Updating Attitude Files')
        attitudeKeys = [key for key in data_dict_attitude.keys()]
        if 'ILat' and 'ILong' in attitudeKeys:
            data_dict_attitude['ILat'][0] = ILat
            data_dict_attitude['ILong'][0] = ILong
        else:
            data_dict_attitude = {**data_dict_attitude, **{'ILat':
                                             [ILat, {'LABLAXIS': f'ILat_{int(refAlt)}km',
                                                          'DEPEND_0': 'Epoch',
                                                          'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                          'UNITS': 'deg',
                                                          'VALIDMIN': ILat.min(), 'VALIDMAX': ILat.max(),
                                                          'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]}}

            data_dict_attitude = {**data_dict_attitude, **{'ILong':
                                                               [ILong, {'LABLAXIS': f'ILong_{int(refAlt)}km',
                                                                       'DEPEND_0': 'Epoch',
                                                                       'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                       'FORMAT': 'E12.2',
                                                                       'UNITS': 'deg',
                                                                       'VALIDMIN': ILong.min(), 'VALIDMAX': ILong.max(),
                                                                       'VAR_TYPE': 'support_data',
                                                                       'SCALETYP': 'linear'}]}}

        outputCDFdata(outputPath=inputFilesAttitude, data_dict=data_dict_attitude, instrNam='attitude')
        Done(start_time)


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if updateCDFfile:

        # --- --- --- --- --- --
        # --- LOAD THE FILES ---
        # --- --- --- --- --- --
        prgMsg(f'Loading data from {input_names_searchable[wFile]}')
        data_dict, globalAttrs = loadDictFromFile(inputFilePath=inputFiles[wFile], getGlobalAttrs=True)

        data_dict_attitude = loadDictFromFile(inputFilePath=inputFilesAttitude)

        try:
            dataEpoch = data_dict['Epoch'][0]
        except:
            dataEpoch = data_dict['fixed_Epoch'][0]

        ILat = data_dict_attitude['ILat'][0]
        ILong = data_dict_attitude['ILong'][0]
        Done(start_time)

        if len(dataEpoch) > len(data_dict_attitude['Epoch'][0]):
            # --- interpolate ILat/ILong onto dataset time base ---
            prgMsg('Interpolating ILat/ILong')
            data_dict_attitudeInterp = InterpolateDataDict(InputDataDict=data_dict_attitude,
                                                           InputEpochArray=data_dict_attitude['Epoch'][0],
                                                           wKeys= ['ILat', 'ILong'],
                                                           targetEpochArray=dataEpoch)
            newILat = np.array(data_dict_attitudeInterp["ILat"][0])
            newILong = np.array(data_dict_attitudeInterp["ILong"][0])

        elif len(dataEpoch) < len(data_dict_attitude['Epoch'][0]):
            prgMsg('DownSampling Ilat/ILong')
            # downsample ILatILong to fit data
            downSampleIndicies = np.array([np.abs(data_dict_attitude['Epoch'][0] - dataEpoch[i]).argmin() for i in range(len(dataEpoch))])
            newILat = np.array(ILat[downSampleIndicies])
            newILong = np.array(ILong[downSampleIndicies])
            Done(start_time)
        else:
            newILat = np.array(ILat)
            newILong = np.array(ILong)


        if outputData:

            prgMsg('Creating output file')

            data_dict = {**data_dict, **{'ILat': [newILat, {'LABLAXIS': f'ILat_{int(refAlt)}km',
                                                                       'DEPEND_0': 'Epoch',
                                                                       'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                       'FORMAT': 'E12.2',
                                                                       'UNITS': 'deg',
                                                                       'VALIDMIN': newILat.min(), 'VALIDMAX': newILat.max(),
                                                                       'VAR_TYPE': 'support_data',
                                                                       'SCALETYP': 'linear'}]}}

            data_dict = {**data_dict, **{'ILong': [newILong, {'LABLAXIS': f'ILong_{int(refAlt)}km',
                                                                        'DEPEND_0': 'Epoch',
                                                                        'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                        'FORMAT': 'E12.2',
                                                                        'UNITS': 'deg',
                                                                        'VALIDMIN': newILong.min(),
                                                                        'VALIDMAX': newILong.max(),
                                                                        'VAR_TYPE': 'support_data',
                                                                        'SCALETYP': 'linear'}]}}

            outputPath = inputFiles[wFile]

            outputCDFdata(outputPath, data_dict, globalAttrsMod=globalAttrs, instrNam=globalAttrs['Descriptor'])

            Done(start_time)





# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if generateILatILong or justPrintFileNames:
        ILatILong_Include(wRocket, rocketFolderPath, justPrintFileNames, 0)
    elif wFiles == []:
        for file in range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf'))):
            ILatILong_Include(wRocket, rocketFolderPath, justPrintFileNames, file)
    else:
        for file in wFiles:
            ILatILong_Include(wRocket, rocketFolderPath, justPrintFileNames, file)
