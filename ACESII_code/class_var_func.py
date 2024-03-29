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
from ACESII_code import data_paths
from cdflib import cdfread
import os
from os import environ,remove
from scipy.signal import butter, filtfilt
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
Re = 6357 # radius of earth in meter
m_e = 9.11 * 10**(-31)
q0 = 1.602176565 * 10**(-19)
kB = 1.380649 * 10**(-23)
cm_to_m = 100
erg_to_eV = 6.2415E11 # how many eVs are in 1 erg
IonMasses = [1.67E-27, 2.6567E-26, 1.6738E-27, 6.646477E-27, 5.3134E-26, 4.9826E-26, 2.3259E-26] # proton
ionNames =  ['proton',       'O+',       'H+',        'He+',      'O2+',      'NO+',       'N+']
ep0 = 8.8541878128E-12 # permittivity of free space
u0 = 4*np.pi*(10**(-7))
lightSpeed = 299792458
coordinatesSets = [['_east','_north','_up'],['_x','_y','_z'],['_e','_p','_r']]
coordinatesNames = ['ENU','RktFrm','Field_Aligned']


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



#####################################################################
# ----- FUNCTIONS -----
#####################################################################
lat_to_meter = 111.319488  # 1 deg latitude to kilometers on Earth

def butterworth(lowcutoff, highcutoff, fs, order, filtertype):
    if filtertype.lower() == 'bandpass':
        return butter(N = order, Wn = [lowcutoff, highcutoff], fs=fs, btype='bandpass')
    elif filtertype.lower() == 'lowpass':
        return butter(N=order, Wn = highcutoff, fs=fs, btype='lowpass')
    elif filtertype.lower() == 'highpass':
        return butter(N=order, Wn = lowcutoff, fs=fs, btype='highpass')
    elif filtertype.lower() == 'bandstop':
        return butter(N=order, Wn = [lowcutoff, highcutoff], fs=fs, btype='bandstop')
def butter_filter(data, lowcutoff, highcutoff, fs, order,filtertype):
    b, a = butterworth(lowcutoff, highcutoff, fs, order, filtertype)
    y = filtfilt(b, a, data)
    return y

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


    for i,set in enumerate(coordinatesSets):

        tempCoords = []

        for key in keys:

            for coordStr in set:
                if coordStr in key.lower():
                    tempCoords.append(key)

        if len(tempCoords) == 3:
            coordcompNames = tempCoords
            coordSetName = coordinatesNames[i]
            coordSet = coordinatesSets[i]

    return coordcompNames,coordSetName, coordSet

def long_to_meter(long, lat):
    return long*(lat_to_meter * np.cos(np.radians(lat)))

def meter_to_long(long_km, lat_km):
    latDeg = lat_km/lat_to_meter
    return (long_km/lat_to_meter) * 1/(np.cos(np.radians(latDeg)))

def calculateLong_to_meter(Lat): # determines the meters/long conversion for each input lattitude using earth's radius as a perfect sphere
    return (np.pi/180) * Re * np.cos(np.radians(Lat))

def prgMsg(message):
    print(f'{message}: ',end='')

def Done(start_time):
    print(f'{color.GREEN}Done{color.END} at {color.YELLOW}{round(time.time() - start_time,1)} seconds {color.END}' )
def setupPYGMT():
    environ["GMT_LIBRARY_PATH"] = data_paths.CDF_LIB

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

        lowerIndex,higherIndex = np.abs(output_data_dict[targetVar[1]][0] - targetVar[0][0]).argmin(),np.abs(output_data_dict[targetVar[1]][0] - targetVar[0][1]).argmin()

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
    if os.path.exists(outputPath):
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




# def MaxwellianParticles(Te, Zpos, N): # returns an array of particle velocities [[v_perp1, v_perp2, ... v_perpN],[v_par1, v_par2, ..., v_parN]] for N particles with distribution temperature T


def mSSA_components(data_dict_input, compNames, SSA_window_Size, mirrorData, **kwargs):

    # Functionality Statement:
    # [1] Input a data_dict
    # [2] Identify the components to mSSA
    # [3] If you need to mirror the data first, do so
    # [4] Calculate the mSSA components
    # [5] return a new data_dict with the Epoch variable also inside

    # get kwargs
    mirrorPercent = kwargs.get('mirrorPercent', 0.2)

    # --- Imports ---
    from ACESII_code.supportCode.Support_Libraries.pymssa import MSSA
    from pandas import DataFrame
    from copy import deepcopy

    # --- Ensure all the compNames keys are in the data_dict_input ---
    for key in compNames:
        if key not in data_dict_input:
            raise Exception('input Variable Key not in data dictionary')

    # --- Ensure all the Epoch variable is in the data_dict_input ---
    if 'Epoch' not in data_dict_input:
        raise Exception(r'No variable named "Epoch" in input dictionary')

    # --- Mirror the data if requested ---
    # in order to eliminate the effect of mSSA has on the edge of the dataset when stitching
    # things together, we concatinate 20% of the mirrored data to the sides of the dataset
    if mirrorData:
        for key in compNames:
            inputData = np.array(deepcopy(data_dict_input[key][0]))
            addThisMuchIndex = int(mirrorPercent*len(inputData))
            firstHalf = list(np.flip(np.array(inputData[0:addThisMuchIndex])))
            backHalf = list(np.flip(np.array(inputData[-addThisMuchIndex:])))
            mirroredData = np.array( firstHalf + list(inputData) + backHalf)
            data_dict_input[key][0] = mirroredData

    # --- --- --- --- ----
    # --- PERFORM mSSA ---
    # --- --- --- --- ----

    # create the MSSA object
    mssa = MSSA(n_components=None, window_size=SSA_window_Size, verbose=False)

    # convert data_dict input Data to pandas dataframe
    data = DataFrame({compNames[i]: data_dict_input[compNames[i]][0] for i in range(len(compNames))})

    # calculate the mSSA
    mssa.fit(data)

    # get the mSSA components
    components = mssa.components_

    # --- --- --- --- --- ---
    # --- output the data ---
    # --- --- --- --- --- ---

    # creat the output data_dict and populate it
    data_dict_output = {}

    for i in range(len(compNames)):
        dataToOutput = np.array(components[i, :, :])  # update the data for output to be the components
        attrs = deepcopy(data_dict_input[compNames[i]][1])
        attrs['LABLAXIS'] = compNames[i]
        attrs['VALIDMIN'] = dataToOutput.min()
        attrs['VALIDMAX'] = dataToOutput.max()
        data_dict_output = {**data_dict_output, **{compNames[i]: [dataToOutput, attrs]}}

    # add in the Epoch variable
    data_dict_output = {**data_dict_output, **{'Epoch': deepcopy(data_dict_input['Epoch'])}}

    return data_dict_output

def generateGrouping(SSA_window_Size, badCompIndicies, InvestigateIndicies, noiseLimit):

    if SSA_window_Size == 0: # if the file I'm working on isn't an SSAcomp file
        return []
    else: # make the proper grouping
        # first plot every component i.e. the original data
        group = [[i for i in range(SSA_window_Size * 3)]]

        # plot all the identified bad inidices
        grouping = [badCompIndicies]

        # investigate other components
        grouping += InvestigateIndicies
        # grouping += [[i] for i in range(9)]

        # the "noise" data
        grouping += [[i for i in range(noiseLimit, 3 * SSA_window_Size)]]

        # # get the good stuff
        goodStuff = []
        for i in range(noiseLimit):
            if i not in grouping[0]:
                goodStuff.append(i)

        grouping += [goodStuff]

        return group + grouping

def mSSA_grouping_Plotting(B_SSA, groupings, SSA_window_Size, compNames, wAxesSSA, wRocket,Epoch,mirrored):

    prgMsg('Grouping mSSA elements')

    # some preamble
    from numpy.fft import rfft, fftfreq
    spinFreq = 0.6442441031179716 if wRocket == 4 else 0.55

    # --- Plot the FFT and groupings for one wAxes axes ---
    fig, ax = plt.subplots(nrows=len(groupings), ncols=2)

    rktName = 'ACESII 36359\n' if wRocket == 4 else 'ACESII 36364\n'
    fig.suptitle(rktName+
                 compNames[wAxesSSA] +
                 f'\n Window Length: {SSA_window_Size}\n'+
                 f'{Epoch[0].strftime("%H:%M:%S")} to {Epoch[-1].strftime("%H:%M:%S")}')


    # loop over all the groups in grouping
    for i in range(len(groupings)):
        data = np.array(B_SSA[wAxesSSA])

        # combine the groupings
        plotThisData = np.array([0 for i in range(len(data[:, 0]))], dtype='float64')

        for j, componentIndex in enumerate(groupings[i]):
            plotThisData += np.array(data[:, componentIndex])

        # plot the component
        ax[i, 0].plot([i for i in range(len(plotThisData))], plotThisData)

        # if the data is mirrored, plot vertical bars on the LHS plots showing the cutoff
        if mirrored:
            # determine where the true data starts using the length of the Epoch Variable
            startPoint = int(0.2*len(Epoch))
            endPoint = len(B_SSA[0]) - int(0.2*len(Epoch))
            ax[i, 0].axvline(x=startPoint, color='black')
            ax[i, 0].axvline(x=endPoint, color='black')

        if i == 0:
            ax[i, 0].set_ylabel('Original')
            ax[i, 1].set_ylabel('Original')
        elif i == 1:
            ax[i, 0].set_ylabel('Harmonics')
            ax[i, 1].set_ylabel('Harmonics')
        elif i == len(groupings) - 2:
            ax[i, 0].set_ylabel('Noise')
            ax[i, 1].set_ylabel('Noise')
        elif i == len(groupings) - 1:
            ax[i, 0].set_ylabel('Physical Signal')
            ax[i, 1].set_xlabel('Frequency [Hz]')
            ax[i, 1].set_ylabel('FFT')
        else:
            ax[i, 0].set_ylabel('F{}'.format(i))
            ax[i, 1].set_ylabel(groupings[i])

        # calculate the FFT and plot it
        N, T = len(plotThisData), 1 / 128
        yf, xf = rfft(plotThisData), fftfreq(N, T)[:N // 2]
        FFT = 2.0 / N * np.abs(yf[0:N // 2])
        ax[i, 1].plot(xf, FFT)
        ax[i, 1].vlines([spinFreq * (i + 1) for i in range(50)], ymin=min(FFT), ymax=max(FFT), alpha=0.5,
                        color='red')
        ax[i, 1].set_xlim(-0.1, 10)
        ax[i, 1].set_ylim(0, max(FFT))

    plt.show()

# def mSSA_wCorMatrix_Plotting(B_SSA, grouping, SSA_window_Size,compNames):
#
#     # --- collect the group info for the LAST (noise) grouping ---
#
#     prgMsg('Plotting wCor Matrix (WARNING: THIS CAN TAKE A LONG TIME):')
#     from ACESII_code.supportCode.Support_Libraries.pymssa import MSSA
#     from pandas import DataFrame
#
#     mssa = MSSA(n_components=None, window_size=SSA_window_Size, verbose=False)
#
#     # convert data_dict input Data to pandas dataframe
#     data = DataFrame({compNames[i]: data_dict_input[compNames[i]][0] for i in range(len(compNames))})
#
#     # calculate the mSSA
#     mssa.fit(data)
#
#     # plot all three axes correlation matrix
#     mssa.fit(DataFrame(
#         {'Data1': [i for i in range(len(data_dict_SSAcomps['Epoch'][0]))],
#          'Data2': [i for i in range(len(data_dict_SSAcomps['Epoch'][0]))],
#          'Data3': [i for i in range(len(data_dict_SSAcomps['Epoch'][0]))]}
#     ))
#
#     for i in range(3):
#         mssa.components_[i, :, :] = data_dict_SSAcomps[compNames[i]][0]
#
#     for i in range(3):
#         # calculate correlation matrix
#         wcorr = np.abs(mssa.w_correlation(mssa.components_[i, :, :]))
#
#         # plot it
#         plt.title(compNames[i])
#         ax = plt.imshow(wcorr, cmap='turbo')
#         plt.xlabel(r"$\tilde{F}_i$")
#         plt.ylabel(r"$\tilde{F}_j$")
#         plt.colorbar(ax.colorbar, fraction=0.045)
#         ax.colorbar.set_label("$W_{i,j}$")
#         plt.clim(0, 1)
#         plt.show()



# --- The Basic rotation matricies ---
def Rx(angle):
    angleRad = np.radians(angle)
    return np.array([[1,0,0],
                     [0,np.cos(angleRad),-np.sin(angleRad)],
                     [0,np.sin(angleRad),np.cos(angleRad)]])
def Ry(angle):
    angleRad = np.radians(angle)
    return np.array([[np.cos(angleRad),0,np.sin(angleRad)],
                     [0,1,0],
                     [-np.sin(angleRad),0,np.cos(angleRad)]])
def Rz(angle):
    angleRad = np.radians(angle)
    return np.array([[np.cos(angleRad),-np.sin(angleRad),0],
                     [np.sin(angleRad),np.cos(angleRad),0],
                     [0,0,1]])


def R_roll(angle):
    angleRad = np.radians(angle)
    return np.array([[1,                0,                0],
                     [0, np.cos(angleRad), np.sin(angleRad)],
                     [0,-np.sin(angleRad), np.cos(angleRad)]])
def R_pitch(angle):
    angleRad = np.radians(angle)
    return np.array([[np.cos(angleRad), 0, -1*np.sin(angleRad)],
                     [0,                1, 0],
                     [np.sin(angleRad), 0, np.cos(angleRad)]])
def R_yaw(angle):
    angleRad = np.radians(angle)
    return np.array([[np.cos(angleRad),    np.sin(angleRad), 0],
                     [-1*np.sin(angleRad), np.cos(angleRad), 0],
                     [0,                   0,                1]])

def DCM(roll,pitch,yaw):
    return np.matmul(R_yaw(yaw),np.matmul(R_pitch(pitch),R_roll(roll)))

def ENUtoECEF(Lat,Long):
    angleLat = np.radians(Lat)
    angleLong = np.radians(Long)

    R = np.array([
        [-np.sin(angleLong), -np.cos(angleLong)*np.sin(angleLat), np.cos(angleLong)*np.cos(angleLat)],
        [ np.cos(angleLong), -np.sin(angleLong)*np.sin(angleLat), np.sin(angleLong)*np.cos(angleLat)],
        [0,                   np.cos(angleLat),                   np.sin(angleLat)]
    ])

    return R

def sphereToCartesian(r,theta,phi):
    thetaRad = np.radians(theta)
    phiRad = np.radians(phi)

    R = np.array([
        [np.sin(thetaRad)*np.cos(phiRad), np.cos(thetaRad)*np.cos(phiRad), -np.sin(phiRad)],
        [np.sin(thetaRad)*np.sin(phiRad), np.cos(thetaRad)*np.sin(phiRad),  np.cos(phiRad)],
        [               np.cos(thetaRad),                np.sin(thetaRad),               0]
    ])
    return R

def Rotation3D(yaw, pitch, roll):
    yawR = np.radians(yaw)
    pitchR = np.radians(pitch)
    rollR = np.radians(roll)
    return np.array([[np.cos(yawR)*np.cos(pitchR), np.cos(yawR)*np.sin(pitchR)*np.sin(rollR) - np.sin(yawR)*np.cos(rollR), np.cos(yawR)*np.sin(pitchR)*np.cos(rollR) + np.sin(yawR)*np.sin(rollR) ], [np.sin(yawR)*np.cos(pitchR), np.sin(yawR)*np.sin(pitchR)*np.sin(rollR) + np.cos(yawR)*np.cos(rollR), np.sin(yawR)*np.sin(pitchR)*np.cos(rollR) - np.cos(yawR)*np.sin(rollR)], [-1*np.sin(pitchR), np.cos(pitchR)*np.sin(rollR), np.cos(pitchR)*np.cos(rollR)]])

def RotationAboutAxes(theta, axX,axY,axZ):
    thetaR = np.radians(theta)
    return np.array([
    [np.cos(thetaR) + (axX*axX)*(1 - np.cos(thetaR)),    axX*axY*(1-np.cos(thetaR)) - axZ*np.sin(thetaR), axX*axZ*(1-np.cos(thetaR)) + axY*np.sin(thetaR)],
    [axY*axX*(1- np.cos(thetaR)) + axZ*np.sin(thetaR), np.cos(thetaR) + axY*axY*(1- np.cos(thetaR)),    axY*axZ*(1- np.cos(thetaR)) - axX*np.sin(thetaR)],
    [axZ*axX*(1- np.cos(thetaR)) - axY*np.sin(thetaR), axZ*axY*(1- np.cos(thetaR))+ axX*np.sin(thetaR), np.cos(thetaR)+axZ*axZ*(1- np.cos(thetaR))]
        ])

def GreatCircleDistance(lat1,lat2,long1,long2):
    return 2*Re*np.arcsin(np.sqrt( np.sin(np.radians( (lat2-lat1)/2  ))**2 + np.cos(np.radians(lat1))*np.cos(np.radians(lat2))*np.sin(np.radians( (long2-long1)/2  ))**2  ))

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

def dateTimetoTT2000(InputEpoch,inverse):

    if inverse: # tt2000 to datetime
        if isinstance(InputEpoch[0], dt.datetime):
            raise Exception(TypeError, "Input Epoch Array is datetime array!")
        else:
            return np.array([pycdf.lib.tt2000_to_datetime(tme) for tme in InputEpoch])
    else: # datetime to tt2000
        if isinstance(InputEpoch[0], (int, float, complex)):
            raise Exception(TypeError, "Input Epoch Array is TT2000 array!")
        else:
            return np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in InputEpoch])


def EpochTo_T0_Rocket(InputEpoch, T0):

    # Convert the T0
    if isinstance(T0, dt.datetime):  # Input Epoch is datetime array
        startPoint = pycdf.lib.datetime_to_tt2000(T0)
    elif isinstance(T0, (int, float, complex)):
        if T0 < 1E6: # if T0 is an integter but small-ish, then it's not a tt2000
            raise Exception('T0 is not a TT2000 Value')
        else: # If T0 is already a TT2000
            startPoint = T0

    # Convert the Epoch
    if isinstance(InputEpoch[0],dt.datetime): # Input Epoch is datetime array
        return np.array([pycdf.lib.datetime_to_tt2000(tme) - startPoint for tme in InputEpoch])/1E9
    elif isinstance(InputEpoch[0], (int, float, complex)): # Input Epoch is tt2000 array
        return (np.array(InputEpoch) - startPoint)/1E9


def InterpolateDataDict(InputDataDict,InputEpochArray,wKeys,targetEpochArray):

    # InputDataDict --> Contains a data_dict of the data which will be interpolated onto the new dataset
    # InputEpoch --> The epoch that InputDataDict uses
    # wKeys --> Keys of the variables in InputDataDict that we want to interpolate. If wKeys == [], do all the keys
    # targetEpoch --> Epoch that the data will be interpolated onto (MUST BE IN TT2000)

    from scipy.interpolate import CubicSpline
    import datetime as dt

    # get the keys to interpolate
    if wKeys == []:
        wKeys = [key for key, val in InputDataDict.items()]

    # Ensure the inputEpoch is in tt2000
    if isinstance(InputEpochArray[0], dt.datetime):
        InputEpochArray = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in InputEpochArray])
    if isinstance(targetEpochArray[0], dt.datetime):
        targetEpochArray = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in targetEpochArray])


    # --- Do the interpolation ---
    data_dict_interpolated = {}

    # interpolate over all the keys and store them in new dictonary
    for key in wKeys:
        if 'Epoch'.lower() not in key.lower():

            # --- cubic interpolation ---
            splCub = CubicSpline(InputEpochArray, InputDataDict[key][0])

            # --- evaluate the interpolation at all the new Epoch points ---
            newData = np.array([splCub(timeVal) for timeVal in targetEpochArray])

            # --- store the data in the interpolated data_dict ---
            data_dict_interpolated = {**data_dict_interpolated, **{key:[newData,InputDataDict[key][1]]}}

        else:
            newEpoch = np.array([pycdf.lib.tt2000_to_datetime(tme) for tme in targetEpochArray])
            data_dict_interpolated = {**data_dict_interpolated, **{key:[newEpoch,InputDataDict[key][1]]}}

    return data_dict_interpolated




#############################
# ----- MODEL FUNCTIONS -----
#############################

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


    # --- Interpolate input altitude onto dataset ---
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






