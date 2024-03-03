# --- imports ---
import numpy as np
from numpy import linspace
from numpy import array

# --- USEFUL CONSTANTS ---
m_to_km = 1E3
R_REF = 6371.2 * m_to_km  # in meters

######################
# ---GENERAL SETUP ---
######################
class GenToggles:
    simLen = 200 # how many delta T steps to simulate
    deltaT = 0.01 # in seconds
    simAltLow = 200*m_to_km # low altitude (in meters)
    simAltHigh = 20000*m_to_km # high altitude (in meters)
    obsHeight = 400*m_to_km # height of observation (in meters)
    alt_Rez = 2000
    simAlt = linspace(simAltLow, simAltHigh, alt_Rez)  # in METERS
    simTime = linspace(0, deltaT*simLen, simLen+1)  # in seconds
    simColors = ['tab:purple', 'tab:orange', 'tab:red', 'tab:blue', 'tab:green', 'tab:brown', 'tab:pink'] # the color choices available for the simulation to use
    simFolderPath = r'C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\Science\AlfvenSingatureAnalysis\Simulations\TestParticle'
    simOutputPath = r'C:\Data\ACESII\science\simulations\TestParticle'


################################
# --- PARTICLE DISTRIBUTIONS ---
################################
class ptclToggles:
    seedChoice = 10 # some value to define the randomness seed
    ptclTemperature = 10 # distribution temperature in eV
    Z0_ptcl_ranges = array([0.6]) * m_to_km * 6371 # initial altitude of particles (in meters)
    N_ptcls = 400  # number of particles. Example: The real data at s3 has 10598 particles
    ptcl_mass = 9.11 * 10 ** (-31)  # mass of the particle
    ptcl_charge = 1.602176565 * 10 ** (-19)  # charge of the particle
    simEnergyRanges = [[0.01, 1], [1, 5], [5, 10], [10, 15], [15, 30], [30, 60]]  # The range of energies for each color (the rightmost linspace point is ignored)
    totalNumberOfParticles = N_ptcls * len(Z0_ptcl_ranges)

###########################
# --- GEOMAGNETIC FIELD ---
###########################
class BgeoToggles:
    Lshell = 9

########################
# --- ELECTRIC FIELD ---
########################
class EToggles:
    Z0_wave = (R_REF) # initial altitude of the wave (in meters)
    waveFreq_Hz = 4 # in Hz
    waveFreq_rad = 2*np.pi*waveFreq_Hz
    waveFraction = 2 # What fraction of the initial bipolar wave we want to keep. e.g. 2 --> Half the wave, 3 --> 1/3 of wave etc

    # lambda_perp
    tau0 = 0.3 # risetime of pulse
    lambdaPerp0 = 2 * m_to_km  # in meters
    kperp0 = 2*np.pi / lambdaPerp0 # this is Kperp AT THE OBSERVATION POINT
    Eperp0 = 5 # V/m
    lambdaPerp_Rez = 11 # resolution of the x-direction (MUST BE ODD)
    static_Kperp = False
    flipEField = True

