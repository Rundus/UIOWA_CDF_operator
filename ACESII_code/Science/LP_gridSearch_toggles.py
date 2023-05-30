# --- LP_gridSearch_toggles.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to put the toggles for the LP GridSearch technique

import numpy as np

#############################
# --- GRID SEARCH TOGGLES ---
#############################
Ti = 0.3 # value of ion temperature, in eV

Vsp_range = np.linspace(0.4, 1.6, 30)
voltageStartPoint = 0.2  # Only take data > than this voltage
rounding = 9  # rounding on the fit parameters

N = 100 # number of points used to perform grid search
a0_range = np.linspace(0, 2.11E7, N)
a1_range = np.linspace(0, 2.11E7, N)
a2_range = np.linspace(0, 2.11E7, N)

transGuess = [1E15, 1, 0.2]
transBounds = ([1E10, 0, 0], [np.Inf, 2.5, 100]) # Bounds for: a0 = nT_e, a1 = Vsp, a2 = Te
satGuess = [1E15, 1, 0.2]
satBounds = ([1E10, 0, 0], [np.Inf, 2.5, 100]) # Bounds for: a0 = nT_e, a1 = Vsp, a2 = Te

# Transition region fitting params
transFitParameters = dict(maxfev=1000000, bounds=transBounds, p0=transGuess)

# saturation region Fitting params
satFitParameters = dict(maxfev=1000000, bounds=satBounds, p0=satGuess)
