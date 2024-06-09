# --- EvanModel.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Use the Evans 1974 model to validate integration techniques AND calculate the
# secondary/backscatter responses for our data


# --- IMPORTS ---

import matplotlib.pyplot as plt
import numpy as np
from ACESII_code.myImports import *
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from ACESII_code.class_var_func import EpochTo_T0_Rocket, q0,m_e
from scipy.interpolate import griddata
start_time = time.time()





##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
prgMsg('Loading Data')

rocketAttrs, b, c = ACES_mission_dicts()

# EEPAA Distribution Data
inputEEPAA_diffFlux_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0]
data_dict_diffFlux = loadDictFromFile(inputEEPAA_diffFlux_high, targetVar=targetVar, wKeys_Reduce=['Differential_Energy_Flux','Differential_Number_Flux', 'Epoch'])

# Define the data
diffEFlux = deepcopy(data_dict_diffFlux['Differential_Energy_Flux'][0])
diffNFlux = deepcopy(data_dict_diffFlux['Differential_Number_Flux'][0])
Epoch = deepcopy(data_dict_diffFlux['Epoch'][0])
Energy = deepcopy(data_dict_diffFlux['Energy'][0])
Pitch = deepcopy(data_dict_diffFlux['Pitch_Angle'][0])
Done(start_time)