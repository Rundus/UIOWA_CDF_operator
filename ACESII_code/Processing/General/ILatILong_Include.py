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
wRocket = 5
# inputPath_modifier = 'L3\Langmuir' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier = 'L2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder

wFiles = [1,2,5,6]
refAlt = 150 # represents 150 km reference altitude that everything is tied to
# ---------------------------
generateILatILong = False # Calculates and Stores the ILat and ILong variables as a .cdf File
plotILatILong = False
footPrintDiffernece = False # reads in BOTH attitude files and updates both files with a "ILat/ILong difference" variable
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
        Done(start_time)
        if plotILatILong:
            fig, ax = plt.subplots()
            ax.plot(Long, Lat, color='blue')
            ax.plot(ILong,ILat, color='red')
            plt.show()
            Done(start_time)
        # --- --- --- --- --- --- ----
        # --- UPDATE ATTITUDE FILE ---
        # --- --- --- --- --- --- ----
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

    if footPrintDiffernece:
        # load the data
        inputFilesAttitude_high = glob(f'{rocketFolderPath}attitude\high\*.cdf')[0]
        inputFilesAttitude_low = glob(f'{rocketFolderPath}attitude\low\*.cdf')[0]

        data_dict_attitude_high,globalAttrsHigh = loadDictFromFile(inputFilePath=inputFilesAttitude_high, getGlobalAttrs=True)
        data_dict_attitude_low,globalAttrsLow = loadDictFromFile(inputFilePath=inputFilesAttitude_low, getGlobalAttrs=True)

        # determine what index the low flyer data starts w.r.t. the high flyer
        startIndex = np.abs(data_dict_attitude_high['Epoch'][0] - data_dict_attitude_low['Epoch'][0][0]).argmin()

        footPrint_Latkm_Difference = np.zeros(shape=len(data_dict_attitude_high['Epoch'][0]))
        footPrint_Longkm_Difference = np.zeros(shape=len(data_dict_attitude_high['Epoch'][0]))
        footPrint_total_Difference = np.zeros(shape=len(data_dict_attitude_high['Epoch'][0]))
        footPrint_time_Difference = np.zeros(shape=len(data_dict_attitude_high['Epoch'][0]))


        # --- Spatial Difference ---
        lat_km = [ [lat_to_meter*data_dict_attitude_high['ILat'][0][i] for i in range(len(data_dict_attitude_high['Epoch'][0]))],
                   [lat_to_meter*data_dict_attitude_low['ILat'][0][i] for i in range(len(data_dict_attitude_low['Epoch'][0]))]]

        long_km = [[long_to_meter(data_dict_attitude_high['ILong'][0][i], data_dict_attitude_high['ILat'][0][i]) for i in range(len(data_dict_attitude_high['Epoch'][0]))],
                  [long_to_meter(data_dict_attitude_low['ILong'][0][i], data_dict_attitude_low['ILat'][0][i]) for i in range(len(data_dict_attitude_low['Epoch'][0]))]]

        for i in range(len(data_dict_attitude_high['Epoch'][0])):

            if i < startIndex:
                latDif = lat_km[0][i] - lat_km[1][0]
                longDif = long_km[0][i] - long_km[1][0]
                totalDif = np.sqrt(latDif**2 + longDif**2)

            elif startIndex <= i <= startIndex + len(data_dict_attitude_low['Epoch'][0])-1:
                latDif = lat_km[0][i] - lat_km[1][i-startIndex]
                longDif = long_km[0][i] - long_km[1][i-startIndex]
                totalDif = np.sqrt(latDif**2 + longDif**2)

            elif i > startIndex + len(data_dict_attitude_low['Epoch'][0]):
                latDif = lat_km[0][i] - lat_km[1][-1]
                longDif = long_km[0][i] - long_km[1][-1]
                totalDif = np.sqrt(latDif**2 + longDif**2)

            footPrint_Latkm_Difference[i]=latDif
            footPrint_Longkm_Difference[i]=longDif
            footPrint_total_Difference[i]=totalDif

        # --- Temporal Difference ---
        for i in range(len(data_dict_attitude_high['ILat'][0])):

            HF_ILat_time = pycdf.lib.datetime_to_tt2000(data_dict_attitude_high['Epoch'][0][i])

            LF_sawSameILat_index = np.abs(data_dict_attitude_low['ILat'][0]-data_dict_attitude_high['ILat'][0][i]).argmin()

            timeBetweenSeeingSameILat = (HF_ILat_time - pycdf.lib.datetime_to_tt2000(data_dict_attitude_low['Epoch'][0][LF_sawSameILat_index]))/1E9
            footPrint_time_Difference[i] = timeBetweenSeeingSameILat



        data_dict_attitude_high = {**data_dict_attitude_high, **{'footPrint_ILat_km_Diff': [footPrint_Latkm_Difference, {'LABLAXIS': f'footPrint_ILat_km_Diff',
                                                    'DEPEND_0': 'Epoch',
                                                    'FILLVAL': rocketAttrs.epoch_fillVal,
                                                    'FORMAT': 'E12.2',
                                                    'UNITS': 'km',
                                                    'VALIDMIN': footPrint_Latkm_Difference.min(), 'VALIDMAX': footPrint_Latkm_Difference.max(),
                                                    'VAR_TYPE': 'support_data',
                                                    'SCALETYP': 'linear'}]}}

        data_dict_attitude_high = {**data_dict_attitude_high, **{'footPrint_ILong_km_Diff': [footPrint_Longkm_Difference, {'LABLAXIS': f'footPrint_ILong_km_Diff',
                                                                    'DEPEND_0': 'Epoch',
                                                                    'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                    'FORMAT': 'E12.2',
                                                                    'UNITS': 'km',
                                                                    'VALIDMIN': footPrint_Longkm_Difference.min(),
                                                                    'VALIDMAX': footPrint_Longkm_Difference.max(),
                                                                    'VAR_TYPE': 'support_data',
                                                                    'SCALETYP': 'linear'}]}}

        data_dict_attitude_high = {**data_dict_attitude_high, **{'footPrint_total_km_Diff': [footPrint_total_Difference, {'LABLAXIS': f'footPrint_total_km_Diff',
                                                                      'DEPEND_0': 'Epoch',
                                                                      'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                      'FORMAT': 'E12.2',
                                                                      'UNITS': 'km',
                                                                      'VALIDMIN': footPrint_total_Difference.min(),
                                                                      'VALIDMAX': footPrint_total_Difference.max(),
                                                                      'VAR_TYPE': 'support_data',
                                                                      'SCALETYP': 'linear'}]}}

        data_dict_attitude_high = {**data_dict_attitude_high, **{
            'footPrint_lattitude_time_Difference': [footPrint_time_Difference, {'LABLAXIS': f'footPrint_lattitude_time_Difference',
                                                                     'DEPEND_0': 'ILat',
                                                                     'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                     'FORMAT': 'E12.2',
                                                                     'UNITS': 'Seconds',
                                                                     'VALIDMIN': footPrint_time_Difference.min(),
                                                                     'VALIDMAX': footPrint_time_Difference.max(),
                                                                     'VAR_TYPE': 'support_data',
                                                                     'SCALETYP': 'linear'}]}}

        outputPath = inputFilesAttitude_high
        outputCDFdata(outputPath= outputPath, data_dict=data_dict_attitude_high, globalAttrsMod=globalAttrsHigh)

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

            data_dict = {**data_dict, **{'Alt': [newILat, {'LABLAXIS': f'ILat_{int(refAlt)}km',
                                                            'DEPEND_0': 'Epoch',
                                                            'FILLVAL': rocketAttrs.epoch_fillVal,
                                                            'FORMAT': 'E12.2',
                                                            'UNITS': 'deg',
                                                            'VALIDMIN': newILat.min(), 'VALIDMAX': newILat.max(),
                                                            'VAR_TYPE': 'support_data',
                                                            'SCALETYP': 'linear'}]}}

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
            try:
                outputCDFdata(outputPath, data_dict, globalAttrsMod=globalAttrs, instrNam=globalAttrs['Descriptor'])
            except:
                outputCDFdata(outputPath, data_dict, globalAttrsMod=globalAttrs)

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
