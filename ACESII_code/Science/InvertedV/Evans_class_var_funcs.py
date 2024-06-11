import numpy as np
from numpy import power, cos, radians, exp, sin,pi,abs,arctan2
from ACESII_code.class_var_func import m_e,q0
import datetime as dt


###################
# --- --- --- --- -
# --- VARIABLES ---
# --- --- --- --- -
###################
datasetReduction_TargetTime = [dt.datetime(2022,11, 20, 17, 24, 50, 000000), dt.datetime(2022,11,20,17,25,15,000000)]

###################
# --- --- --- --- -
# --- FUNCTIONS ---
# --- --- --- --- -
###################

def loadDiffNFluxData():
    from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
    from glob import glob
    from ACESII_code.class_var_func import loadDictFromFile
    from copy import deepcopy
    targetVar = [datasetReduction_TargetTime, 'Epoch']

    rocketAttrs, b, c = ACES_mission_dicts()

    # EEPAA Distribution Data
    inputEEPAA_diffFlux_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0]
    data_dict_diffFlux = loadDictFromFile(inputEEPAA_diffFlux_high, targetVar=targetVar, wKeys_Reduce=['Differential_Energy_Flux', 'Differential_Number_Flux', 'Epoch'])

    # Define the data
    diffNFlux = deepcopy(data_dict_diffFlux['Differential_Number_Flux'][0])
    Epoch = deepcopy(data_dict_diffFlux['Epoch'][0])
    Energy = deepcopy(data_dict_diffFlux['Energy'][0])
    Pitch = deepcopy(data_dict_diffFlux['Pitch_Angle'][0])

    return diffNFlux,Epoch,Energy,Pitch

# def diffNFlux_for_mappedMaxwellian(x, n, T, beta, V, alpha):
#     Vpara_sqrd = (2 * x * power(cos(radians(alpha)), 2) / m_e) - 2 * V / m_e + (1 - 1 / beta) * (2 * x / m_e) * (power(sin(radians(alpha)), 2))
#     Vperp_sqrd = ((2 * x) / (beta * m_e)) * power(sin(radians(alpha)), 2)
#
#     return (2 * x) * ((q0 / m_e) ** 2) * (1E2 * n) * power(m_e / (2 * pi * q0 * T), 3 / 2) * exp((-m_e / (2 * T)) * (Vpara_sqrd + Vperp_sqrd))


def dist_Maxwellian(Vperp,Vpara,n,T):
    # Input: density [cm^-3], Temperature [eV], Velocities [m/s]
    # output: the distribution function in SI units [s^3 m^-6]
    Emag = (0.5*m_e*(Vperp**2 + Vpara**2))/q0
    return (1E6 * n) * power(m_e / (2 * pi * q0 * T), 3 / 2) * exp(-1*Emag/T)
    # return (n) * power(m_e / (2 * pi * q0 * T), 3 / 2) * exp(-1*Emag/T)


def calc_diffNFlux(Vperp,Vpara,dist):

    # Input: Velocities [m/s], distribution function [s^3m^-6]
    # output: diffNFlux [cm^-2 s^-1 eV^-1 str^-1]
    Emag = 0.5 * m_e * (Vperp**2 + Vpara**2)/q0
    return (2 * Emag) * power((100*q0 / m_e),2) * dist
    # return (2 * Emag) * power((1E-2*q0 / m_e),2) * dist


def velocitySpace_to_PitchEnergySpace(EnergyBins, PitchBins, VperpGrid, VparaGrid,ZGrid):

    # description:
    # INPUT: Two velocity space grids  + Z-value grid
    # OUTPUT: Energy and Pitch grids + new Z-value grid

    # determine the type of input data
    VperpValues = VperpGrid.flatten()
    VparaValues = VparaGrid.flatten()
    ZgridValues = ZGrid.flatten()

    ZGrid_New = [[[] for engy in range(len(EnergyBins))] for ptch in range(len(PitchBins))]
    calcEnergies = [0.5 * m_e * (perp ** 2 + par ** 2) / q0 for perp, par in zip(VperpValues,VparaValues)]
    calcPitch = [(180 / pi) * arctan2(perp, par) for perp, par in zip(VperpValues,VparaValues)]

    # assign the values to ZGrid_new
    for i in range(len(ZgridValues)):
        engyIdx = abs(EnergyBins - calcEnergies[i]).argmin()
        ptchIdx = abs(PitchBins - calcPitch[i]).argmin()
        ZGrid_New[ptchIdx][engyIdx].append(ZgridValues[i])

    # flatten the values in the diffnFlux new array
    for ptch in range(len(PitchBins)):
        for engy in range(len(EnergyBins)):
            try:
                ZGrid_New[ptch][engy] = sum(ZGrid_New[ptch][engy]) / len(ZGrid_New[ptch][engy])
                # diffNFlux_model[tme][ptch][engy] = sum(diffNFlux_model[tme][ptch][engy])
            except:
                ZGrid_New[ptch][engy] = sum(ZGrid_New[ptch][engy])

    EnergyGrid,PitchGrid = np.meshgrid(EnergyBins,PitchBins)
    return np.array(ZGrid_New), EnergyGrid, PitchGrid