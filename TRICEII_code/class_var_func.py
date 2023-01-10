# --- class_var_func.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to store all the classes/variables/functions

import numpy as np,time
import data_paths
from cdflib import cdfread
from os import environ
from bisect import bisect_left

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
        if isinstance(data,list):
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

def setupPYCDF():
    environ['HOMEDRIVE'] = data_paths.HOMEDRIVE
    environ['HOMEPATH'] = data_paths.HOMEPATH
    environ["CDF_LIB"] = data_paths.CDF_LIB

def take_closest(myList, myNumber):
    # Assumes myList is sorted. Returns index to closest value to myNumber.

    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return pos
    if pos == len(myList):
        return pos

    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return pos
    else:
       return pos - 1

# ---------------------
# ----- VARIABLES -----
# ---------------------
t_noise_start1_low = 10089 # corresponds to TRICE_low 08:31:10.014
t_noise_end1_low = 15289 # corresponds to TRICE_low 08:35:30.016
t_noise_start2_low = 22789 # corresponds to TRICE_low 08:41:45.018
t_noise_end2_low = 23589 # corresponds to TRICE_low 08:42:25.18
t_noise_start1_high = 7610# corresponds to TRICE_high 08:28:59.99
t_noise_end1_high = 14090 # corresponds to TRICE_high 08:34:23.99
t_noise_start2_high = 20510 # corresponds to TRICE_high 08:39:44.994
t_noise_end2_high = 22410 # corresponds to TRICE_high 08:41:19.995
noise_region_before_event_in_counts = [[7610, 14090], [10089, 15289]]
noise_region_after_event_in_counts = [[20510, 22410], [22789, 23589]]

m_e = 9.11 * 10**(-31)
q = 1.602176565 * 10**(-19)
cm_to_m = 100

def tmCDF_TRICE(flier):
    return [Initialize_cdfFile(data_paths.TRICE_tmCDF_files[flier][i]) for i in range(len(data_paths.TRICE_tmCDF_files[flier]))]

def L0_TRICE(flier):
    return [Initialize_cdfFile(data_paths.TRICE_L0_files[flier][i]) for i in range(len(data_paths.TRICE_L0_files[flier]))]

def L1_TRICE(flier):
    return [Initialize_cdfFile(data_paths.TRICE_L1_files[flier][i]) for i in range(len(data_paths.TRICE_L1_files[flier]))]

def L2_TRICE(flier):
    return [Initialize_cdfFile(data_paths.TRICE_L2_files[flier][i]) for i in range(len(data_paths.TRICE_L2_files[flier]))]

def Mag_TRICE(flier):
    return [Initialize_cdfFile(data_paths.TRICE_mag_files[flier][i]) for i in range(len(data_paths.TRICE_mag_files[flier]))]

def csv_attitude_TRICE(flier):
    return [data_paths.TRICE_attitude_csv_files[flier][i] for i in range(len(data_paths.TRICE_attitude_csv_files[flier]))]

def cdf_attitude_TRICE(flier):
    return [data_paths.TRICE_attitude_cdf_files[flier][i] for i in range(len(data_paths.TRICE_attitude_cdf_files[flier]))]

