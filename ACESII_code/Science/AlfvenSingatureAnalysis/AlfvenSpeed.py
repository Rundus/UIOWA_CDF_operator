# --- ModelAlfvenSpeed.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Plot a model of the Alfven Speed vs altitude in order to analyze
# Alfven dispersion features


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import copy
import itertools
# --- --- --- --- ---

import time

import matplotlib.pyplot as plt

from ACESII_code.class_var_func import Done, setupPYCDF

start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# assume the kinetic term is just 2^-0.5
simplifyKinetic = False

# Altitude Range (in KILOMETERS) to plot alfven speed over
wkperp = (2*3.1415926535)/100 # This is ANGULAR wave number: k = 2 pi/lambda
AltLow = 1
AltHigh = 22000
xscaling= (1E3) # how to scale the xaxis of the plot
yscaling= 299792458 # how to scale the yaxis of the plot. Nominally: (1E4)*(1E3)
plotylabel = 'Kinetic Alfven Speed [$V_{A}/c$]' # nominally: [10,000 km/s]
plotxlabel = 'Altitude (1000 km)'
ylimits = (0,0.05)
xlimits = (22,0) # the x-axis is inverted remember

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import pyIGRF
import datetime as dt
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg,outputCDFdata,L1_TRICE_Quick,L2_TRICE_Quick,loadDictFromFile, q0,m_e,u0,lightSpeed,IonMasses,ep0,Re,cm_to_m
from glob import glob
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


# profile taken from Kletzing: Alfven Wave Electron Acceleration
def density(z): # returns density for altitude "z [km]" in cm^-3
    h = 0.06*Re # in km from E's surface
    n0 = 6E4
    n1 = 1.34E7
    z0 = 0.05*Re # in km from E's surface
    n = n0*np.exp(-1*(z-z0)/h) + n1*(z**(-1.55)) # calculated density (in cm^-3)
    return (cm_to_m**3)*n

def kineticTerm(kperp, z, simplify): # represents the denominator of the Alfven velocity term: 1/(1 + (kperp*c/omega_pe)^2)^1/2
    if simplify:
        y = 1/np.sqrt(2)
    else:
        plasmaFreq = np.sqrt((density(z) * q0 * q0) / (m_e * ep0))
        y = 1 / np.sqrt(1 + (kperp * lightSpeed / plasmaFreq) ** 2)
    return y

def AlfvenSpeed(z,lat,long,year,kperp,simplify):
    # -- Output order forpyIGRF.igrf_value ---
    # [0] Declination (+ E | - W)
    # [1] Inclination (+ D | - U)
    # [2] Horizontal Intensity
    # [3] North Comp (+ N | - S)
    # [4] East Comp (+ E | - W)
    # [5] Vertical Comp (+ D | - U)
    # [6] Total Field

    B = pyIGRF.igrf_value(lat, long, z, year)
    V_A = (B[6]*1E-9)/np.sqrt(u0 * density(z) * IonMasses[0])

    if simplify:
        V = V_A*kineticTerm(1, z, simplify)
    else:
        V = V_A*kineticTerm(kperp, z, simplify)
    return V


def main(AltLow, AltHigh):

    # altitude range
    N = 10000
    altitudeAxis = np.linspace(AltLow, AltHigh, N)

    # fig,ax = plt.subplots()
    # ax.plot(altitudeAxis,density(altitudeAxis))
    # plt.show()

    # determine the alfven speed
    year = 2022 + 323 / 365  # Corresponds to 11/20/2022
    LaunchLat = 69.294167
    LaunchLong = 16.020833
    alfvenAxis = np.sqrt(2)*np.array([AlfvenSpeed(z, LaunchLat, LaunchLong, year, wkperp, simplifyKinetic) for z in altitudeAxis])

    # normalize the data to whatever you wish
    xData = altitudeAxis/(xscaling) # to 1000 km
    yData = alfvenAxis/(yscaling) # to 10,000 km/s

    # --- --- --- ---
    # --- PLOTTING ---
    # --- --- --- ---
    fig, ax = plt.subplots()
    ax.set_ylabel(plotylabel)
    ax.set_xlabel(plotxlabel)
    fig.suptitle('Alfven Speed vs Altitude')
    ax.plot(xData, yData)
    ax.invert_xaxis()
    ax.set_ylim(*ylimits)
    ax.set_xlim(*xlimits)
    plt.show()







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
main(AltLow,AltHigh)