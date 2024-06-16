# --- conversions.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to store all the functions I used to load CDFs

# IMPORTS
from os import environ,remove,path
from numpy import abs
from ACESII_code import data_paths
def setupPYCDF():
    environ['HOMEDRIVE'] = data_paths.HOMEDRIVE
    environ['HOMEPATH'] = data_paths.HOMEPATH
    environ["CDF_LIB"] = data_paths.CDF_LIB

setupPYCDF()
from spacepy import pycdf

# VARIABLES
coordinatesSets = [['_east','_north','_up'],['_x','_y','_z'],['_e','_p','_r']]
coordinatesNames = ['ENU','RktFrm','Field_Aligned']

def getInputFiles(rocketFolderPath,wRocket,inputPath_modifier,**kwargs):

    from glob import glob
    modifier = kwargs.get('modifier', '')

    from ACESII_code.data_paths import fliers

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{modifier}\{fliers[wRocket - 4]}\*.cdf')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{modifier}\{fliers[wRocket - 4]}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace(inputPath_modifier.lower() + '_', '') for ifile in input_names]

    return inputFiles,input_names,input_names_searchable

def getCoordinateKeys(data_dict):

    keys = [key for key in data_dict.keys()]

    coordcompNames = []
    coordSetName = []
    coordSet = []

    for i, set in enumerate(coordinatesSets):

        tempCoords = []

        for key in keys:

            for coordStr in set:
                if coordStr in key.lower():
                    tempCoords.append(key)

        if len(tempCoords) == 3:
            coordcompNames = tempCoords
            coordSetName = coordinatesNames[i]
            coordSet = coordinatesSets[i]

    return coordcompNames, coordSetName, coordSet



def loadDictFromFile(inputFilePath, **kwargs):
    input_data_dict = kwargs.get('input_data_dict', {})
    globalAttrs = {}
    targetVar = kwargs.get('targetVar', [])
    getGlobalAttrs = kwargs.get('getGlobalAttrs', False)
    reduceData = True if targetVar != [] else kwargs.get('reduceData', False)
    wKeys_Load = kwargs.get('wKeys_Load', [])
    wKeys_Reduce = kwargs.get('wKeys_Reduce', [])

    # load the data dict
    with pycdf.CDF(inputFilePath) as inputDataFile:
        for key, val in inputDataFile.attrs.items():
            globalAttrs = {**globalAttrs, **{str(key): str(val)}}

        for key, val in inputDataFile.items():
            input_data_dict = {**input_data_dict, **{key: [inputDataFile[key][...], {key: val for key, val in inputDataFile[key].attrs.items()}]}}

    # load only the data in wKeys_Load
    output_data_dict = {}
    wKeys_Load = [key for key in input_data_dict.keys()] if wKeys_Load == [] else wKeys_Load
    for key in wKeys_Load:
        output_data_dict = {**output_data_dict, **{key:input_data_dict[key]}}

    # reduce the data
    if reduceData:
        try:
            h = input_data_dict[targetVar[1]]
        except:
            raise Exception(f'no Target Variable found: {targetVar[1]}')

        lowerIndex,higherIndex = abs(output_data_dict[targetVar[1]][0] - targetVar[0][0]).argmin(),abs(output_data_dict[targetVar[1]][0] - targetVar[0][1]).argmin()

        # determine which keys to reduce
        wKeys_Reduce = [key for key in output_data_dict.keys()] if wKeys_Reduce == [] else wKeys_Reduce
        for key in wKeys_Reduce:
            output_data_dict[key][0] = output_data_dict[key][0][lowerIndex:higherIndex]


    if getGlobalAttrs:
        return output_data_dict,globalAttrs
    else:
        return output_data_dict

def outputCDFdata(outputPath, data_dict, **kwargs):

    ModelData = kwargs.get('ModelData', [])
    globalAttrsMod = kwargs.get('globalAttrsMod', {})
    instrNam = kwargs.get('instrNam', 'None')

    # if ModelData == []:
    #     ModelData = L2_ACES_Quick(0)
    if globalAttrsMod == {}:
        from ACESII_code.missionAttributes import ACES_mission_dicts
        rocketAttrs, b, c = ACES_mission_dicts()
        globalAttrsMod = rocketAttrs.globalAttributes[0]

    # --- delete output file if it already exists ---
    if path.exists(outputPath):
        remove(outputPath)

    # --- open the output file ---
    with pycdf.CDF(outputPath, '') as sciFile:
        sciFile.readonly(False)
        inputGlobDic = globalAttrsMod

        # # --- write out global attributes ---
        # try:
        #     inputGlobDic = ModelData.cdfFile.globalattsget()
        # except:
        #     inputGlobDic = ModelData



        for key, val in inputGlobDic.items():
            if key == 'Descriptor':
                globalAttrsMod[key] = instrNam
            if key in globalAttrsMod:
                sciFile.attrs[key] = globalAttrsMod[key]
            else:
                sciFile.attrs[key] = val

        # --- WRITE OUT DATA ---
        for varKey, varVal in data_dict.items():

            if varKey in ['Epoch', 'Epoch_monitors', 'Epoch_esa']:  # epoch data
                sciFile.new(varKey, data=varVal[0], type=33)
            elif 'Function' in varKey:
                sciFile.new(varKey, data=varVal[0], type=pycdf.const.CDF_REAL8)
            else:  # other data
                sciFile.new(varKey, data=varVal[0])

            # --- Write out the attributes and variable info ---
            for attrKey, attrVal in data_dict[varKey][1].items():
                if attrKey == 'VALIDMIN':
                    sciFile[varKey].attrs[attrKey] = varVal[0].min()
                elif attrKey == 'VALIDMAX':
                    sciFile[varKey].attrs[attrKey] = varVal[0].max()
                elif attrVal != None:
                    sciFile[varKey].attrs[attrKey] = attrVal