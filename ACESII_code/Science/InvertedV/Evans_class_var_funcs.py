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

def diffNFlux_for_mappedMaxwellian(x, n, T, beta, V, alpha):
    Vpara_sqrd = (2 * x * power(cos(radians(alpha)), 2) / m_e) - 2 * V / m_e + (1 - 1 / beta) * (2 * x / m_e) * (power(sin(radians(alpha)), 2))
    Vperp_sqrd = ((2 * x) / (beta * m_e)) * power(sin(radians(alpha)), 2)

    return (2 * x) * ((q0 / m_e) ** 2) * (1E2 * n) * power(m_e / (2 * pi * q0 * T), 3 / 2) * exp((-m_e / (2 * T)) * (Vpara_sqrd + Vperp_sqrd))


def dist_Maxwellian(Vperp,Vpara,n,T):
    # Input: density [cm^-3], Temperature [eV], Velocities [m/s]
    # output: the distribution function in SI units [s^3 m^-6]
    Emag = (0.5*m_e*(Vperp**2 + Vpara**2))/q0
    return (1E6 * n) * power(m_e / (2 * pi * q0 * T), 3 / 2) * exp(-1*Emag/T)


def calc_diffNFlux(Vperp,Vpara,dist):
    # Input: Velocities [m/s], distribution function [s^3m^-6]
    # output: diffNFlux [cm^-2 s^-1 eV^-1 str^-1]
    Emag = 0.5 * m_e * (Vperp**2 + Vpara**2)/q0
    return (2 * Emag) * power(q0 /( 100*m_e),2) * dist


def calc_DistributionMapping(Vperp_gridVals,Vpara_gridVals,model_T, model_n, model_V0, beta,modifyInitialBeam,beamPitchThreshold,beamEnergyThreshold):

    # --- Define a grid a velocities (static) ---
    VperpGrid, VparaGrid = np.meshgrid(Vperp_gridVals, Vpara_gridVals)
    distGrid = dist_Maxwellian(VperpGrid, VparaGrid, n=model_n, T=model_T)

    if modifyInitialBeam:
        for i in range(len(VperpGrid)):
            for j in range(len(VperpGrid[0])):

                pitchVal = np.degrees(np.arctan2(VperpGrid[i][j] ,VparaGrid[i][j]))
                EnergyVal = 0.5*m_e*(VperpGrid[i][j] **2 + VparaGrid[i][j]**2) / q0

                if np.abs(pitchVal) >= beamPitchThreshold:
                    distGrid[i][j] = 0
                if EnergyVal >= beamEnergyThreshold:
                    distGrid[i][j] = 0

    diffNFluxGrid = calc_diffNFlux(VperpGrid, VparaGrid, distGrid)

    # --- Determine the Accelerated Velocities ---
    Vperp_gridVals_Accel = Vperp_gridVals
    Vpar_gridVals_Accel = np.array([np.sqrt(val ** 2 + 2 * model_V0 * q0 / m_e) for val in Vpara_gridVals])
    VperpGrid_Accel, VparaGrid_Accel = np.meshgrid(Vperp_gridVals_Accel, Vpar_gridVals_Accel)
    diffNFluxGrid_Accel = calc_diffNFlux(VperpGrid_Accel, VparaGrid_Accel, distGrid)

    # --- Determine the new velocities at different beta ---
    VperpArray_magsph = VperpGrid_Accel.flatten()
    VparaArray_magsph = VparaGrid_Accel.flatten()
    Vperp_gridVals_iono = np.array([np.sqrt(beta) * val for val in VperpArray_magsph])
    Vpara_gridVals_iono_sqrd = np.array([Vpar_magsph ** 2 + (1 - beta) * (Vper_magsph ** 2) for Vper_magsph, Vpar_magsph in zip(VperpArray_magsph, VparaArray_magsph)])
    Vpara_gridVals_iono = np.array([np.sqrt(val) if val >= 0 else -1 * np.sqrt(np.abs(val)) for val in Vpara_gridVals_iono_sqrd])


    # make the grids
    VperpGrid_iono = Vperp_gridVals_iono.reshape(len(Vperp_gridVals), len(Vpara_gridVals))
    VparaGrid_iono = Vpara_gridVals_iono.reshape(len(Vperp_gridVals), len(Vpara_gridVals))
    diffNFluxGrid_iono = calc_diffNFlux(VperpGrid_iono, VparaGrid_iono, distGrid)

    return distGrid,VperpGrid, VparaGrid, diffNFluxGrid, VperpGrid_Accel, VparaGrid_Accel, diffNFluxGrid_Accel, VperpGrid_iono, VparaGrid_iono, diffNFluxGrid_iono



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