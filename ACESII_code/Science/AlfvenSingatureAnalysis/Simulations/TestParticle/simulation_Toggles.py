# --- imports ---
from numpy import linspace

# --- USEFUL CONSTANTS ---
m_to_km = 1E3
R_REF = 6371.2 * m_to_km  # in meters

############################
# --- simulation toggles ---
############################
simLen = 330 # how many delta T steps to simulate
deltaT = 0.01 # in seconds
simAltLow = 300*m_to_km # low altitude in meters
simAltHigh = 12000*m_to_km # high altitude in meters
obsHeight = 400*m_to_km # height of observation (in meters)

simAlt = linspace(simAltLow, simAltHigh, 2000)  # in METERS
simTime = [deltaT * i for i in range(simLen)]