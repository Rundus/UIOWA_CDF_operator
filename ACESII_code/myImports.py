# --- myImports.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: There are some common imports that every file uses. In order to de-clutter my code
# I can place these imports here. Only the imports which EVERY file uses will go here.

#########################
# --- IMPORTS IMPORTS ---
#########################
import scipy
import itertools
import time
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime as dt

######################
# --- FROM IMPORTS ---
######################
from myspaceToolsLib.physicsVariables import IonMasses,u0,q0,m_e,cm_to_m,kB,ep0,Re
from myspaceToolsLib.CDF_load import getCoordinateKeys
from myspaceToolsLib.time import dateTimetoTT2000
from myspaceToolsLib.interpolate import InterpolateDataDict
from ACESII_code.class_var_func import Done,setupPYCDF

from tqdm import tqdm
from glob import glob
from os.path import getsize
from scipy.optimize import curve_fit
from copy import deepcopy
from myspaceToolsLib.conversions import Re
from myspaceToolsLib.CDF_load import loadDictFromFile,outputCDFdata,getCoordinateKeys

from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg
from ACESII_code.class_var_func import  L1_TRICE_Quick,L2_TRICE_Quick,L0_TRICE_Quick,L0_ACES_Quick,L1_ACES_Quick
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap

#####################
# --- SETUP PYCDF ---
#####################
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)