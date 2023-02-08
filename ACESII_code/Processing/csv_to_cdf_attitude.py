# --- csv_to_cdf_attitude.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Turn the .cdf files of the ACESII attitude solution data into cdf files


# --- --- --- --- ---
import time
from ACESII_code.class_var_func import Done, setupPYCDF,prgMsg
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

wRocket = 1
# 0 --> ACESII High Flyer
# 1 --> ACESII Low Flyer


Directories = [r'D:\Data\ACESII\attitude\high\\',r'D:\Data\ACESII\attitude\low\\']
inputFilenames = [r'36359_Synchronized.csv',r'36364_Synchronized.csv']
outputFilenames = [r'36359_Attitude_Synchronized.cdf',r'36364_Attitude_Synchronized.cdf']
inputDirectory = Directories[wRocket]
inputFileName = inputFilenames[wRocket]
outputDirectory = inputDirectory
outputFileName = outputFilenames[wRocket]

wRow_names = 5 -1 # Which row contains the variable names
wRow_units = 6 - 1 # Which row contains the variable units
wRow_startData = 8 - 1 # Which row contains the first row of actual data


# --- special modifications to the variables already in the csv file ---
#  in a dictionary format [ {string name of var: {string name of attributeKey1: attributeval1, attributeKey2: attributeval2  }},{...},...  ]
# Example: specialMods = [time: {'Units':'ns','FORMAT':'I5'} , Alt: {'Units':'km','FORMAT':'I5'} ]
# Defaults to attributes = None if specialMods == [], except for VALIDMIN/VALIDMAX


specialMods = {
    'Alt':               {'DEPEND_0':'Epoch','LABLAXIS':'Alt'},
    'Epoch':             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -9223372036854775808,'FORMAT': 'I5','UNITS': 'ns','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','MONOTON':'INCREASE','TIME_BASE':'J2000','TIME_SCALE':'Terrestrial Time','REFERENCE_POSITION':'Rotating Earth Geoid','SCALETYP':'linear','LABLAXIS':'Epoch'}
}



# --- --- --- ---
# --- import ---
# --- --- --- ---
import numpy as np
from os import remove,path
from csv import reader
from ACESII_code.data_paths import fliers, ACES_data_folder,ACES_csv_trajectories,TRICE_data_folder,ACES_L0_files,TRICE_L0_files
from ACESII_code.class_var_func import color, L2_TRICE_Quick
from ACESII_code.missionAttributes import ACES_mission_dicts,TRICE_mission_dicts
from copy import deepcopy
setupPYCDF()
from spacepy import pycdf
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field

print(color.BOLD + color.CYAN + 'csv_to_cdf.py' + color.END + color.END)

def csv_to_cdf(inputFile,outputFile,missionDicts,specialMods,wRow_names,wRow_units,wRow_startData):

    # collect the entire csv data
    with open(inputFile) as csvFile:
        csvAllData = [row for row in reader(csvFile)]

    varNames = np.array(csvAllData[wRow_names])
    varNames = np.array([nam.replace(' ','') for nam in varNames])


    varUnits = np.array(csvAllData[wRow_units])
    csvData = np.array(csvAllData[wRow_startData:], dtype='float64').transpose()

    # --- Create Data Dict ---
    data_dict = {}
    rocketAttrs, missionDict, data_dict_temp = missionDicts
    exampleAttrs = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808,
                    'FORMAT': 'I5', 'UNITS': None, 'LABLAXIS':None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                    'SCALETYP': 'linear'}

    for index, var in enumerate(varNames):
        varAttrs = deepcopy(exampleAttrs)
        unit = str(varUnits[index]).replace('(','').replace(')','').replace(' ','')
        if unit == 'm':
            unit = 'meters'

        varAttrs['LABLAXIS'] = var
        varAttrs['UNITS'] = unit
        data_dict = {**data_dict,**{varNames[index] : [csvData[index],  varAttrs ]}}



    prgMsg(f'Converting csv Files for {inputFileName}')

    # --- Modify Variable Attributes in data dict to the  ---
    for key, val in specialMods.items():

        if key in data_dict:

            for keyAttr, valAttr in specialMods[key].items():

                if keyAttr in data_dict[key][1]: # modify an existing attribute
                    data_dict[key][1][keyAttr] = valAttr

                else: # append a new attribute
                    data_dict[key][1] = {**data_dict[key][1], **{keyAttr:valAttr}}

        else: # Create a whole new variable
            data_dict = {**data_dict, **{key: [[], deepcopy(val) ]} }


    # Create the Epoch Variable
    Launch_time = rocketAttrs.Launch_Times[wRocket]
    Epoch = np.array([ pycdf.lib.tt2000_to_datetime(int(Launch_time + data_dict['Time'][0][i]*(10**9))) for i in range(len(data_dict['Time'][0])) ])
    data_dict['Epoch'][0] = Epoch
    Done(start_time)

    # --- --- --- --- ---
    # --- Output Data ---
    # --- --- --- --- ---
    prgMsg('Writing Out Data')

    # --- delete output file if it already exists ---
    if path.exists(outputFile):
        remove(outputFile)
    pycdf.lib.set_backward(False)


    # --- open the output file ---
    with pycdf.CDF(outputFile, '') as cdfFile:
        cdfFile.readonly(False)

        # --- write out global attributes ---
        globalAttrsMod = rocketAttrs.globalAttributes[wRocket]
        ModelData = L2_TRICE_Quick(wRocket)
        inputGlobDic = ModelData.cdfFile.globalattsget()

        for key, val in inputGlobDic.items():
            if key == 'Descriptor':
                globalAttrsMod[key] = 'None'
            if key in globalAttrsMod:
                cdfFile.attrs[key] = globalAttrsMod[key]
            else:
                cdfFile.attrs[key] = val

        # --- WRITE OUT DATA ---
        for varKey, varVal in data_dict.items():


            if 'Epoch' in varKey:  # epoch data
                cdfFile.new(varKey, data=varVal[0], type=33)
            else:  # other data
                cdfFile.new(varKey, data=varVal[0])

            # --- Write out the attributes and variable info ---
            for attrKey, attrVal in data_dict[varKey][1].items():
                if attrKey == 'VALIDMIN':
                    cdfFile[varKey].attrs[attrKey] = varVal[0].min()
                elif attrKey == 'VALIDMAX':
                    cdfFile[varKey].attrs[attrKey] = varVal[0].max()
                elif attrVal != None:
                    cdfFile[varKey].attrs[attrKey] = attrVal

    Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
inputFile = inputDirectory+inputFileName
outputFile = outputDirectory + outputFileName
if path.exists(inputFile):
    if wRocket == 0 or wRocket == 1:#ACESII
        missionDicts = ACES_mission_dicts()
    elif wRocket == 2 or wRocket == 3: #TRICEII
        missionDicts = TRICE_mission_dicts()

    csv_to_cdf(inputFile,outputFile,missionDicts,specialMods,wRow_names,wRow_units,wRow_startData)

else:
    print(color.RED + f'{inputFile} does not exist' + color.END)
