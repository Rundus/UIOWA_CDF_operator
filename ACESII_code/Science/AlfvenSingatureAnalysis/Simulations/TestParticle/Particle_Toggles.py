# --- imports ---
from numpy import array,linspace
from simulation_Toggles import m_to_km



##########################
# --- particle toggles ---
##########################

useRandomEnergies = False # will use random energies between the limits set above
customPtches = False
useRandomPtches = False

# initial altitude of particles (in meters)
# Z0_ptcl_ranges = array(
#             [0.7*6371,
#             0.65*6371,
#             0.6*6371,
#             0.55*6371]
#             )*m_to_km
Z0_ptcl_ranges = array([1*6371])*m_to_km

# NOTE: The real data at s3 has 10598 particles
N = 2*64 # number of particles IN EACH ENERGY
# N = 16 # number of particles IN EACH ENERGY
N_total = int(len(Z0_ptcl_ranges)*N)
ptcl_mass = 9.11 * 10**(-31) # mass of the particle
ptcl_charge = 1.602176565 * 10**(-19) # charge of the particle

# The range of energies for each color (the rightmost linspace point is ignored)
simEnergyRanges = [[0.01, 1],
                   [1, 5],
                   [5, 10],
                   [10, 15],
                   [15, 30],
                   [30, 60]
                   ]

# the color choice for each set of particles
simColors = [['tab:purple'],
             ['tab:orange'],
             ['tab:red'],
             ['tab:blue'],
             ['tab:green'],
             ['tab:brown'],
             ['tab:pink']]

# the number of Energies to linspace the energy ranges (MUST BE SAME LENGTH AS simEnergyRanges)
simPtclsPerEngyRange = [20, 20, 16, 14, 8, 6, 6]

# CUSTOM PITCH ANGLES
# which pitch angles should the initial particles be at

customPitchRanges = [[-30,  30, int(3*N/8)],
                     [ 30, 150, int(1*N/8)],
                     [150, 210, int(3*N/8)],
                     [210, 330, int(1*N/8)]]

simPtcls_Initial_Pitch_Range = [item for sublist in customPitchRanges for item in linspace(*sublist)]
