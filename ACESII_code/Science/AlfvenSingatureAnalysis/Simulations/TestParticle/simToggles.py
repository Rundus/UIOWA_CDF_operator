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
    simLen = 330 # how many delta T steps to simulate
    deltaT = 0.01 # in seconds
    simAltLow = 1*m_to_km # low altitude (in meters)
    simAltHigh = 20000*m_to_km # high altitude (in meters)
    obsHeight = 400*m_to_km # height of observation (in meters)
    simAlt = linspace(simAltLow, simAltHigh, 2000)  # in METERS
    simTime = linspace(0, deltaT*simLen, simLen+1)  # in seconds
    simColors = ['tab:purple', 'tab:orange', 'tab:red', 'tab:blue', 'tab:green', 'tab:brown', 'tab:pink'] # the color choices available for the simulation to use


################################
# --- PARTICLE DISTRIBUTIONS ---
################################
class ptclToggles:
    seedChoice = 10 # some value to define the randomness seed
    ptclTemperature = 40 # distribution temperature in eV
    Z0_ptcl_ranges = array([1]) * m_to_km * 6371 # initial altitude of particles (in meters)
    N_ptcls = 10  # number of particles. Example: The real data at s3 has 10598 particles
    ptcl_mass = 9.11 * 10 ** (-31)  # mass of the particle
    ptcl_charge = 1.602176565 * 10 ** (-19)  # charge of the particle
    simEnergyRanges = [[0.01, 1], [1, 5], [5, 10], [10, 15], [15, 30], [30, 60]]  # The range of energies for each color (the rightmost linspace point is ignored)

###########################
# --- GEOMAGNETIC FIELD ---
###########################
class BgeoToggles:
    Lshell = 9

########################
# --- ELECTRIC FIELD ---
########################
class EToggles:
    Z0_wave = (1.3*6371+3000)*m_to_km # initial altitude of the wave (in meters)
    waveFreq_Hz = 4 # in Hz
    waveFreq_rad = 2*np.pi*waveFreq_Hz
    lambdaPerp0 = 2 * m_to_km  # in meters
    peakPotential = 400 # volts
    static_Kperp = True
    flipEField = True

    # --- Static Toggles ---
    Amplitude_static = 0.18/1000  # V/m
    lamda_static = (7500 - 4500)*m_to_km  # in meters
    WaveSpeed_Static = 1000*m_to_km  # in meters/s
    duty = 40 # duty cycle percent (Max is 50)
