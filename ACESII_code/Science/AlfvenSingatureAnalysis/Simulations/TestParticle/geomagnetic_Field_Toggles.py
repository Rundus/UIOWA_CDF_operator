# --- imports ---
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.simulation_Toggles import simAlt
from numpy import array,degrees,arccos
from numpy.linalg import norm
from datetime import datetime
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from ACESII_code.class_var_func import CHAOS
from simulation_Toggles import m_to_km, R_REF


########################################
# --- GENERATE THE B-FIELD & TOGGLES ---
########################################

# --- B-Field Toggles ---
Lshell = 9

def generateBField():
    geomagAlts = [((alt + R_REF) / R_REF) for alt in simAlt]
    geomagLats = array([degrees(arccos(radi / Lshell)) for radi in geomagAlts])
    geomagLongs = array([111.83 for i in range(len(simAlt))])
    times = [datetime(2022, 11, 20, 17, 20, 00, 000) for i in range(len(simAlt))]
    Pos = array([geomagAlts, geomagLats, geomagLongs]).transpose()
    ISOtime = [times[i].isoformat() for i in range(len(times))]
    cvals_MAG = coord.Coords(Pos, 'MAG', 'sph')
    cvals_MAG.ticks = Ticktock(ISOtime, 'ISO')
    cvals_GDZ = cvals_MAG.convert('GEO', 'sph')
    Lat_geo = cvals_GDZ.lati

    # Get the Chaos model
    B = CHAOS(Lat_geo, [15.25 for i in range(len(simAlt))], array(simAlt) / m_to_km, times)
    Bmag = (1E-9) * array([norm(Bvec) for Bvec in B])
    Bgrad = [(Bmag[i + 1] - Bmag[i]) / (simAlt[i + 1] - simAlt[i]) for i in range(len(Bmag) - 1)] + [0]

    return Bmag, Bgrad