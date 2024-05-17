# --- PlotExtra.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Use the CHAOS model to general the Lat/Long/Alt coordinates for "n" B_geo lines and store them in a .cdf file

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
justPrintFileNames = False
wRocket = 4
wFile = 0
modifier = ''
inputPath_modifier = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\CHAOSBgeoProjections' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder



tTimes_high = [dt.datetime(2022,11,20,17,21,00,000000),
               dt.datetime(2022,11,20,17,23,00,000000),
               dt.datetime(2022,11,20,17,25,00,000000),
               dt.datetime(2022,11,20,17,27,00,000000),
               dt.datetime(2022,11,20,17,29,00,000000)]

tTimes_low = [dt.datetime(2022,11,20,17,21,00,000000),
               dt.datetime(2022,11,20,17,23,00,000000),
               dt.datetime(2022,11,20,17,25,00,000000),
               dt.datetime(2022,11,20,17,27,00,000000),
               dt.datetime(2022,11,20,17,29,00,000000)]
tTimes = [tTimes_high,tTimes_low]
refAlts = [0 + 100*i*1000 for i in range(10)]
# ---------------------------
outputData = False
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import CHAOS,long_to_meter,lat_to_meter,ENUtoECEF


def PlotExtra_conjugactconcept(wRocket, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    inputFiles, input_names, input_names_searchable = getInputFiles(rocketFolderPath=rocketFolderPath,wRocket=wRocket,inputPath_modifier=inputPath_modifier)
    fileoutName = f'ACESII_{rocketID}_CHAOS_Bgeo_LatLongAlt.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return


    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    prgMsg(f'Loading data from {inputPath_modifier} Files')
    data_dict_attitude = loadDictFromFile(inputFiles[wFile])
    Done(start_time)

    exampleAttrs = {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                    'FORMAT': 'E12.2', 'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear',
                    'LABLAXIS': None}

    targetTimes = tTimes[wRocket-4]

    # --- --- --- --- --- --- --
    # --- GET THE CHAOS MODE ---
    # --- --- --- --- --- --- --
    prgMsg('Calculating Bgeo Lines')

    B_geo_line_data = []

    for tme in targetTimes:

        # STEP 1: Get B_geo in ENU coordinates
        idx = np.abs(data_dict_attitude['Epoch'][0] - tme).argmin()
        rktLat = data_dict_attitude['Lat'][0][idx]
        rktLong = data_dict_attitude['Long'][0][idx]
        rktLAlt = data_dict_attitude['Alt'][0][idx]

        B_model = CHAOS(lat=np.array([rktLat]),
                        long=np.array([rktLong]),
                        alt=np.array([rktLAlt]),
                        times=np.array([tme]))  # CHAOS in ENU coordinates

        b = np.array(B_model[0]) / np.linalg.norm(B_model[0])


        # STEP 2: Get B_geo in ECEF coordinates
        b_ECEF = np.matmul(b, ENUtoECEF(rktLat, rktLong))

        # STEP 3: Get Position of Rocket in ECEF coordinates
        rkt0_ECEF = np.matmul(b, ENUtoECEF(rktLat, rktLong))



        # STEP 4: Get B_geo in ECEF coordinates

        # STEP 5: Get B_geo points of interest in ECEF coordinates

        # STEP 6: Convert points of interest to Geodedic coordinates and store them






        B_model_Lat = []
        B_model_Long = []
        B_model_Alt = []

        Latkm = lat_to_meter * rktLat
        Longkm = long_to_meter(rktLong, rktLat)

        for refAlt in refAlts:
            t = (refAlt - rktLAlt) / b[2]
            El = Longkm + b[0] * t
            Nl = Latkm + b[1] * t

            B_model_Lat.append(Nl / lat_to_meter)
            B_model_Long.append(meter_to_long(long_km=El, lat_km=Nl))
            B_model_Alt.append(refAlt)

        B_geo_line_data.append([B_model_Lat,B_model_Long,B_model_Alt])

    Done(start_time)



    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        prgMsg('Creating output file')

        data_dict_output={}
        for tme in range(len(B_geo_line_data)):
            data_dict_output = {**data_dict_output, **{f'Line {tme}': [np.array(B_geo_line_data[tme]),deepcopy(exampleAttrs)]}}



        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wRocket-4]}\\{fileoutName}'

        outputCDFdata(outputPath, data_dict_output)

        Done(start_time)





# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

rocketFolderPath = ACES_data_folder


if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    PlotExtra_conjugactconcept(wRocket, rocketFolderPath, justPrintFileNames)

