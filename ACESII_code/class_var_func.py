# --- class_var_func.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to store all the classes/variables/functions

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
# -------------------

import numpy as np,time
from ACESII_code import data_paths
from cdflib import cdfread
from os import environ,path,remove
def setupPYCDF():
    environ['HOMEDRIVE'] = data_paths.HOMEDRIVE
    environ['HOMEPATH'] = data_paths.HOMEPATH
    environ["CDF_LIB"] = data_paths.CDF_LIB

setupPYCDF()
from spacepy import pycdf

# -------------------
# ----- CLASSES -----
# -------------------
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

class Initialize_cdfFile:
    def __init__(self, cdfFile):
        self.Filepath = cdfFile
        self.cdfFile = cdfread.CDF(cdfFile)
        self.info = self.cdfFile.cdf_info()
        self.zVars = self.info['zVariables']

        # construct dataL1 dictionary
        Zvar_data = []
        for variable in self.zVars:
            try:
                Zvar_data.append(self.cdfFile.varget(variable))
            except:
                Zvar_data.append([])

        Data_dict = dict(zip(self.zVars, Zvar_data))

        # initialize all zVars in __init__
        for key, value in Data_dict.items():
            try:
                setattr(self, key, value)
            except:
                print(f'Problem with {key}:')
    def getVarDat(self, variable):
        return self.cdfFile.varget(variable)



class newCDFvar:
    def __init__(self, zVar_str, data, modParam, model_data_object):
        # --- Convert data to a numpy array ---
        if isinstance(data, list):
            self.var_data = np.array(data)
        elif data.dtype.type is np.str_: # exception: dataL1 is ndarray of strings
            self.var_data = list(data[0])
        else:  # dataL1 is ndarray of numbers
            self.var_data = data

        # --- Modify variable attrs/info ---
        self.zVar_str = zVar_str
        self.attrs_types = {'CATDESC': None, 'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FIELDNAM': None,
                            'FILLVAL': None, 'FORMAT': None, 'LABLAXIS': None, 'UNITS': None, 'VALIDMIN': None,
                            'VALIDMAX': None, 'VAR_TYPE': None, 'SCALETYP': None, 'VAR_NOTES': None}
        self.varInfo_types = {'Variable': None, 'Num': None, 'Var_Type': None, 'Data_Type': None, 'Data_Type_Description': None,
                              'Num_Elements': None, 'Num_Dims': None, 'Dim_Sizes': None, 'Sparse': None,'Last_Rec': None,
                              'Rec_Vary': None, 'Dim_Vary': None, 'Pad': None, 'Compress': None, 'Block_Factor': None}

        self.var_info = model_data_object.cdfFile.varinq(zVar_str)
        self.var_attrs = model_data_object.cdfFile.varattsget(zVar_str)



        for key, value in modParam.items(): # Overwrite the model dataL1 if I specify a particular attr/var param in the input of the constructor
            if value != None and key in self.attrs_types:
                self.var_attrs[key] = value
            elif value != None and key in self.varInfo_types:
                self.var_info[key] = value

    def writeToFile(self,outputFile):
        outputFile.write_var(self.var_info, var_attrs = self.var_attrs, var_data= self.var_data)

# ---------------------
# ----- FUNCTIONS -----
# ---------------------
def prgMsg(message):
    print(f'{message}: ',end='')

def Done(start_time):
    print(f'{color.GREEN}Done{color.END} at {color.YELLOW}{round(time.time() - start_time,1)} seconds {color.END}' )



def setupPYGMT():
    environ["GMT_LIBRARY_PATH"] = data_paths.CDF_LIB

def loadCDFdata(inputFiles,wFile):
    data_dict = {}
    with pycdf.CDF(inputFiles[wFile]) as dataFile:
        for key, val in dataFile.items():
            data_dict = {**data_dict, **{key: [dataFile[key][...], {key: val for key, val in dataFile[key].attrs.items()}]}}
    return data_dict
def outputCDFdata(outputPath, data_dict, ModelData,globalAttrsMod,instrNam):

    # --- delete output file if it already exists ---
    if path.exists(outputPath):
        remove(outputPath)

    # --- open the output file ---
    with pycdf.CDF(outputPath, '') as sciFile:
        sciFile.readonly(False)

        # --- write out global attributes ---
        inputGlobDic = ModelData.cdfFile.globalattsget()
        for key, val in inputGlobDic.items():
            if key == 'Descriptor':
                globalAttrsMod[key] = instrNam
            if key in globalAttrsMod:
                sciFile.attrs[key] = globalAttrsMod[key]
            else:
                sciFile.attrs[key] = val

        # --- WRITE OUT DATA ---
        for varKey, varVal in data_dict.items():
            if 'Epoch' in varKey:  # epoch data
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


# --- The Basic rotation matricies
def Rx(angle):
    angleRad = np.radians(angle)
    return np.array([[1,0,0],[0,np.cos(angleRad),-np.sin(angleRad)],[0,np.sin(angleRad),np.cos(angleRad)]])
def Ry(angle):
    angleRad = np.radians(angle)
    return np.array([[np.cos(angleRad),0,np.sin(angleRad)],[0,1,0],[-np.sin(angleRad),0,np.cos(angleRad)]])
def Rz(angle):
    angleRad = np.radians(angle)
    return np.array([[np.cos(angleRad),-np.sin(angleRad),0],[np.sin(angleRad),np.cos(angleRad),0],[0,0,1]])
# ---------------------
# ----- VARIABLES -----
# ---------------------
m_e = 9.11 * 10**(-31)
q0 = 1.602176565 * 10**(-19)
kB = 1.380649 * 10**(-23)
cm_to_m = 100
IonMasses = [1.67 * 10**(-27)] # proton


# --- Model Data ---
def tmCDF_ACES(flier):
    return [Initialize_cdfFile(data_paths.ACES_tmCDF_files[flier][i]) for i in range(len(data_paths.ACES_tmCDF_files))]
def tmCDF_TRICE(flier):
    return [Initialize_cdfFile(data_paths.TRICE_tmCDF_files[flier][i]) for i in range(len(data_paths.TRICE_tmCDF_files))]
def tmCDF_TRICE_Quick(flier):
    return Initialize_cdfFile(data_paths.TRICE_tmCDF_files[flier][0])
def tmCDF_ACES_Quick(flier):
    return Initialize_cdfFile(data_paths.ACES_tmCDF_files[flier][0])
def L0_TRICE(flier):
    return [Initialize_cdfFile(data_paths.TRICE_L0_files[flier][i]) for i in range(len(data_paths.TRICE_L0_files))]
def L0_ACES(flier):
    return [Initialize_cdfFile(data_paths.ACES_L0_files[flier][i]) for i in range(len(data_paths.ACES_L0_files))]
def L0_TRICE_Quick(flier):
    return Initialize_cdfFile(data_paths.TRICE_L0_files[flier][0])
def L0_ACES_Quick(flier):
    return Initialize_cdfFile(data_paths.ACES_L0_files[flier][0])
def L1_TRICE(flier):
    return [Initialize_cdfFile(data_paths.TRICE_L1_files[flier][i]) for i in range(len(data_paths.TRICE_L1_files))]
def L1_ACES(flier):
    return [Initialize_cdfFile(data_paths.ACES_L1_files[flier][i]) for i in range(len(data_paths.ACES_L1_files))]
def L1_TRICE_Quick(flier):
    return Initialize_cdfFile(data_paths.TRICE_L1_files[flier][0])
def L1_ACES_Quick(flier):
    return Initialize_cdfFile(data_paths.ACES_L1_files[flier][0])
def L2_TRICE_Quick(flier):
    return Initialize_cdfFile(data_paths.TRICE_L2_files[flier][0])
def L2_ACES_Quick(flier):
    return Initialize_cdfFile(data_paths.ACES_L2_files[flier][0])
