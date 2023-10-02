# # --- scratchpaper.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering
import matplotlib.pyplot as plt

# the path to the input file
inputFile = 'C:\Data\ACESII\L2\low\ACESII_36364_l2_flight_E_Field.cdf'

# YOU CAN PROBABLY REMOVE THE NEXT TWO LINES
from ACESII_code.class_var_func import setupPYCDF
setupPYCDF()

# imports
from spacepy import pycdf
import numpy as np
import autoplot as ap

# Connor's special pycdf function to format/load data
def loadDictFromFile(inputFilePath,data_dict):
    with pycdf.CDF(inputFilePath) as inputDataFile:
        for key, val in inputDataFile.items():
            data_dict = {**data_dict, **{key: [inputDataFile[key][...], {key: val for key, val in inputDataFile[key].attrs.items()}]}}
    return data_dict

# --- load the data into a file ---
# Note: FORMAT of the data when loaded with "loadDictFromFile" will be:
# data_dict[some_key] = [ [DATA], [Attributes of DATA]]

data_dict = loadDictFromFile(inputFile,{})

Epoch = [pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict['Epoch'][0]]
EpochDiff = [Epoch[i+1] - Epoch[i] for i in range(len(Epoch)-1)]
print(Epoch[0:10])
