# # --- scratchpaper.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering
import matplotlib.pyplot as plt

# the path to the input file
inputFile = 'C:\Data\ACESII\L2\high\ACESII_36359_l2_eepaa_fullCal.cdf'

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


# --- prepare the data for plotting ---
Epoch = [pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict['Epoch_esa'][0]] # tt2000
Epoch_dt = data_dict['Epoch_esa'][0] # datetime
Epoch_dt64 = [np.datetime64(tme) for tme in Epoch_dt] # datetime64
DiffE = data_dict['Differential_Energy_Flux'][0]

# slice the 2D array out of the data
# Note: The data is 3D with 1st D: time, 2nd D: Pitch Angle, 3rd D: Energy
Data = DiffE[:, 2, :] # <-- Should be 10deg pitch angle
dataForPlotting = np.transpose(Data)

print(f'The shape of the data is: {np.shape(Data)}')

# --- --- -- --- --- --- --- ---
# --- MATPLOTLIB EQUIVALENCE ---
# --- --- -- --- --- --- --- ---
fig, ax = plt.subplots()
X,Y = np.meshgrid(Epoch_dt64,data_dict['Energy'][0])
cmap = ax.pcolormesh(X,Y,dataForPlotting,cmap='turbo',vmin=1E6,vmax=1E9)
cbar = plt.colorbar(cmap)
plt.show()

# --- --- -- --- --- ---
# --- AUTOPLOT STUFF ---
# --- --- -- --- --- ---

ap.start()
chosenEpoch = Epoch_dt64
ap.plot(Epoch,data_dict['Energy'][0], dataForPlotting)
a = input('Press key to continue')

# NOTE code DOES work if you change the time variable to "Epoch"

