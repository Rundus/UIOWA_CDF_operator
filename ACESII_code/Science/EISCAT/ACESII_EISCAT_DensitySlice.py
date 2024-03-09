# --- ACESII_EISCAT_DensitySlice.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Uses the ACESII attitde data to determine an expected density profile for the
# Rockets

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
wRocket = 4
modifier = ''
inputPath_modifier = r'science\EISCAT\tromso\UHF' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\EISCAT_ACESII_Slice' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# ---------------------------
outputData = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none


def ACESII_EISCAT_DensitySlice(wRocket, rocketFolderPath, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    ModelData = L2_TRICE_Quick(wflyer)

    inputFiles_EISCAT = glob(rf'{rocketFolderPath}{inputPath_modifier}\*.cdf')[0]
    inputFiles_attitude = glob(f'{rocketFolderPath}attitude\{fliers[wflyer]}{modifier}\*.cdf')[0]

    fileoutName = f'ACESII_{rocketID}_EISCAT_Tromso_DensitySlice.cdf'


    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    prgMsg(f'Loading data EISCAT Files')
    data_dict_attitude = loadDictFromFile(inputFiles_attitude)
    data_dict_EISCAT = loadDictFromFile(inputFiles_EISCAT)
    Done(start_time)

    exampleAttrs = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
     'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear',
     'LABLAXIS': None}

    data_dict = {
        'ne': [ np.zeros(shape=(len(data_dict_attitude['Epoch'][0]))),deepcopy(exampleAttrs)],
        'Ti': [np.zeros(shape=(len(data_dict_attitude['Epoch'][0]))), deepcopy(exampleAttrs)],
        'T_ratio': [np.zeros(shape=(len(data_dict_attitude['Epoch'][0]))), deepcopy(exampleAttrs)],
        'Nu_in': [np.zeros(shape=(len(data_dict_attitude['Epoch'][0]))), deepcopy(exampleAttrs)],
        'Ion_Comp': [np.zeros(shape=(len(data_dict_attitude['Epoch'][0]))), deepcopy(exampleAttrs)], # The number of ions with molecular weight 28 to 32 AS A FRACTION OF TOTAL ELECTRON DENSITY (# of ions/Ne)
        'Op_Comp': [np.zeros(shape=(len(data_dict_attitude['Epoch'][0]))), deepcopy(exampleAttrs)], # composition [O+]/Ne
    }

    # loop over epoch of rocket
    prgMsg('Sampling Tromso EISCAT Data')
    for tme in range(len(data_dict_attitude['Epoch'][0])):

        # get the rockets altitude and epoch
        altRkt = data_dict_attitude['Alt'][0][tme]/1000
        epochValRkt = data_dict_attitude['Epoch'][0][tme]

        # get the relevant point in the EISCAT data
        closestAlt_Index = np.abs(data_dict_EISCAT['range'][0] - altRkt).argmin()
        closestTime_Index =np.abs(data_dict_EISCAT['Epoch'][0] - epochValRkt).argmin()
        closestAlt =  data_dict_EISCAT['range'][0][closestAlt_Index]
        print(altRkt,closestAlt)

        # check if we're 1km away from the cloests EISCAT point. If so, don't include this point
        if np.abs(altRkt - closestAlt) > 10:

            # set everything to zero if we're too far away
            for key in data_dict.keys():
                data_dict[key][0][tme] = 0

        else:
            data_dict['ne'][0][tme] = data_dict_EISCAT['ne'][0][closestTime_Index][closestAlt_Index]
            data_dict['Ti'][0][tme] = data_dict_EISCAT['ti'][0][closestTime_Index][closestAlt_Index]
            data_dict['T_ratio'][0][tme] = data_dict_EISCAT['tr'][0][closestTime_Index][closestAlt_Index]
            data_dict['Nu_in'][0][tme] = data_dict_EISCAT['co'][0][closestTime_Index][closestAlt_Index]
            data_dict['Ion_Comp'][0][tme] = data_dict_EISCAT['pm'][0][closestTime_Index][closestAlt_Index]
            data_dict['Op_Comp'][0][tme] = data_dict_EISCAT['po+'][0][closestTime_Index][closestAlt_Index]

    Done(start_time)
    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        prgMsg('Creating output file')

        # update the data_dict with other known variables
        data_dict = {**data_dict, **{
            'Epoch': deepcopy(data_dict_attitude['Epoch']),
            'Alt' : deepcopy(data_dict_attitude['Alt']),
            'Long': deepcopy(data_dict_attitude['Long']),
            'Lat': deepcopy(data_dict_attitude['Lat']),
            'ILat': deepcopy(data_dict_attitude['ILat']),
            'ILong': deepcopy(data_dict_attitude['ILong'])
        }}

        # update the attributs
        data_dict['ne'][1]['LABLAXIS'] ='Density'
        data_dict['ne'][1]['UNITS'] = 'cm!A-3!N'

        data_dict['Ti'][1]['LABLAXIS'] = 'Ion Temperature'
        data_dict['Ti'][1]['UNITS'] = 'eV'

        data_dict['T_ratio'][1]['LABLAXIS'] = 'Temperature Ratio'
        data_dict['T_ratio'][1]['UNITS'] = 'Te/Ti'

        data_dict['Nu_in'][1]['LABLAXIS'] = 'Ion-Neutral Collision Freq'
        data_dict['Nu_in'][1]['UNITS'] = 'Hz'

        data_dict['Ion_Comp'][1]['LABLAXIS'] = 'Ion Weight Ratio'
        data_dict['Ion_Comp'][1]['UNITS'] = '(Ions mol 28 to 32)/Ne'

        data_dict['Op_Comp'][1]['LABLAXIS'] = 'Oxygen Ratio'
        data_dict['Op_Comp'][1]['UNITS'] = '[O+]/Ne'

        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

        outputCDFdata(outputPath, data_dict, instrNam='EISCAT_Tromso')

        Done(start_time)





# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
if wRocket == 4:  # ACES II High
    wflyer = 0
elif wRocket == 5: # ACES II Low
    wflyer = 1

if len(glob(rf'{rocketFolderPath}\science\EISCAT\tromso\UHF\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    ACESII_EISCAT_DensitySlice(wRocket, rocketFolderPath,wflyer)
