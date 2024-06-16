# --- class_var_func.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to store all the classes/variables/functions

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"


import matplotlib.pyplot as plt
# -------------------
import numpy as np,time

from cdflib import cdfread
from ACESII_code import data_paths
from os import environ
import datetime as dt
def setupPYCDF():
    environ['HOMEDRIVE'] = data_paths.HOMEDRIVE
    environ['HOMEPATH'] = data_paths.HOMEPATH
    environ["CDF_LIB"] = data_paths.CDF_LIB

setupPYCDF()
from spacepy import pycdf

#####################################################################
# ----- CLASSES -----
#####################################################################
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

#####################################################################
# ----- VARIABLES -----
#####################################################################
Re = 6357 # radius of earth in kilometer
m_e = 9.11 * 10**(-31)
q0 = 1.602176565 * 10**(-19)
kB = 1.380649 * 10**(-23)
cm_to_m = 100
erg_to_eV = 6.2415E11 # how many eVs are in 1 erg
IonMasses = [1.67E-27, 2.6567E-26, 1.6738E-27, 6.646477E-27, 5.3134E-26, 4.9826E-26, 2.3259E-26,   4.6518E-26]
ionNames =  ['proton',       'O+',       'H+',        'He+',      'O2+',      'NO+',       'N+',   'N2+']
ep0 = 8.8541878128E-12 # permittivity of free space
u0 = 4*np.pi*(10**(-7))
lightSpeed = 299792458



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



#######################
# ----- FUNCTIONS -----
#######################

def prgMsg(message):
    print(f'{message}: ',end='')

def Done(start_time):
    print(f'{color.GREEN}Done{color.END} at {color.YELLOW}{round(time.time() - start_time,1)} seconds {color.END}' )
def setupPYGMT():
    environ["GMT_LIBRARY_PATH"] = data_paths.CDF_LIB





# def MaxwellianParticles(Te, Zpos, N): # returns an array of particle velocities [[v_perp1, v_perp2, ... v_perpN],[v_par1, v_par2, ..., v_parN]] for N particles with distribution temperature T


def calcChiSquare(yData, fitData, yData_errors, fitData_Errors, nu):
    chiSum = []

    for i in range(len(yData)):

        if (yData_errors[i] + fitData_Errors[i]) != 0:

            chiSum.append( ((yData[i] - fitData[i])**2) / (yData_errors[i] + fitData_Errors[i])  )

    return (1/nu) * sum(chiSum)

# General plot for the ability to quickly plot three axes
def plot3Axis(Epoch, AxesData, ylabels):

    fig, ax = plt.subplots(3)

    for i in range(3):
        ax[i].plot(Epoch,AxesData[:, i])
        ax[i].set_ylabel(ylabels[i])

    ax[2].set_xlabel('Epoch')
    plt.show()









#############################
# ----- MODEL FUNCTIONS -----
#############################


# --- Dipole Magnetic Field ---
def Bdip_mag(Alt_km, Lat_deg):
    B0 = 3.12E-5

    try: # if input data is arrays
        test1 = len(Alt_km)
        test2 = len(Lat_deg)
        colat = [np.radians(90 - lat) for lat in Lat_deg]
        Bdip = [B0 * np.power(Re / (Re + alt), 3) * np.sqrt(1 + 3 * np.power(np.cos(clat), 2)) for alt, clat in zip(Alt_km, colat)]

    except: # if input data is single values
        colat = np.radians(90 - Lat_deg)
        Bdip = B0 * np.power(Re / (Re + Alt_km), 3) * np.sqrt(1 + 3 * np.power(np.cos(colat), 2))

    return Bdip


# --- Model Ionosphere ---
# Inputs: list containing altitudes of interest
# Outputs: list containing lists of plasma paramters at all the input altitudes with format:
# [InputAltitude,rho, Te (in K), Ti (in K), n(O2+), n(N)+), N(O+), n(e)]
# all number density are in cm^-3

def JonesRees1972_Ionosphere(inputAltitudes):

    # --- get the model data ---
    import pandas as pd
    modelFilePath = r'C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\supportCode\IonosphereModels\JonesRees_IonosphereValues.xlsx'
    pandas_dict = pd.read_excel(modelFilePath)
    VariableNams = [thing for thing in pandas_dict]
    modelData = [pandas_dict[key][1:] for key in VariableNams]


    # --- interpolate input altitude onto dataset ---
    interpData = []
    from scipy.interpolate import CubicSpline
    for varNam in VariableNams:
        if varNam.lower() not in 'height' and varNam != 'Variable':
            # --- cubic interpolation ---
            splCub = CubicSpline(pandas_dict['Height'][1:],pandas_dict[varNam][1:])

            # --- evaluate the interpolation at all the new Epoch points ---
            interpData.append(np.array([splCub(hVal) for hVal in inputAltitudes]))


    # calculate rho
    m_O2p= 5.3133E-26
    m_NOp= 4.9826E-26
    m_Op = 2.6566E-26
    rho = m_O2p*np.array(interpData[2]) + m_NOp*np.array(interpData[3]) + m_Op*np.array(interpData[4])

    finalData = [inputAltitudes, rho] + interpData

    return {'Height':[finalData[0], {'UNITS': 'km'}],
            'rho':   [finalData[1], {'UNITS': 'kg/cm^-3'}],
            'T_e':   [finalData[2], {'UNITS': 'Kelvin'}],
            'T_i':   [finalData[3], {'UNITS': 'Kelvin'}],
            'n_O2p':[finalData[4], {'UNITS': 'cm^-3'}],
            'n_NOp':[finalData[5], {'UNITS': 'cm^-3'}],
            'n_Op': [finalData[6], {'UNITS': 'cm^-3'}],
            'n_e':  [finalData[7], {'UNITS': 'cm^-3'}],
            }

def density(z): # returns density for altitude "z [km]" in m^-3
    h = 0.06*Re # in km from E's surface
    n0 = 6E4
    n1 = 1.34E7
    z0 = 0.05*Re # in km from E's surface
    n = n0*np.exp(-1*(z-z0)/h) + n1*(z**(-1.55)) # calculated density (in cm^-3)
    return (cm_to_m**3)*n



def CHAOS(lat, long, alt, times):

    # imports
    import datetime as dt
    from glob import glob
    from chaosmagpy import load_CHAOS_matfile
    from chaosmagpy.data_utils import mjd2000

    FILEPATH_CHAOS = glob(r'C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\supportCode\CHAOS/CHAOS-*.mat')[0]

    R_REF = 6371.2

    # give inputs
    theta = np.array([90 - lat[i] for i in range(len(lat))]) # colat in deg
    phi = np.array(long)
    radius = np.array(alt) + R_REF

    # convert datetime date to mjd2000
    if not isinstance(times[0], dt.date):
        raise Exception('Input times are not datetimes. Convert to python datetime')
    else:
        time = np.array([mjd2000(date.year, date.month, date.day, date.hour) for date in times])  # year, month, day

    # load the CHAOS model
    model = load_CHAOS_matfile(FILEPATH_CHAOS)

    # print('Computing core field.')
    B_core = model.synth_values_tdep(time, radius, theta, phi)

    # print('Computing crustal field up to degree 110.')
    B_crust = model.synth_values_static(radius, theta, phi, nmax=110)

    # complete internal contribution
    B_radius_int = B_core[0] + B_crust[0]
    B_theta_int = B_core[1] + B_crust[1]
    B_phi_int = B_core[2] + B_crust[2]

    # print('Computing field due to external sources, incl. induced field: GSM.')
    B_gsm = model.synth_values_gsm(time, radius, theta, phi, source='all')

    # print('Computing field due to external sources, incl. induced field: SM.')
    B_sm = model.synth_values_sm(time, radius, theta, phi, source='all')

    # complete external field contribution
    B_radius_ext = B_gsm[0] + B_sm[0]
    B_theta_ext = B_gsm[1] + B_sm[1]
    B_phi_ext = B_gsm[2] + B_sm[2]

    # complete forward computation
    B_radius = B_radius_int + B_radius_ext
    B_theta = B_theta_int + B_theta_ext
    B_phi = B_phi_int + B_phi_ext

    # output CHAOS_ENU
    B_ENU = np.array([[B_phi[i],-1*B_theta[i], B_radius[i]] for i in range(len(B_radius))])

    return B_ENU

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

    B = CHAOS(lat, long, z, year)
    V_A = (B[6]*1E-9)/np.sqrt(u0 * density(z) * IonMasses[0])

    if simplify:
        V = V_A*kineticTerm(1, z, simplify)
    else:
        V = V_A*kineticTerm(kperp, z, simplify)
    return V






