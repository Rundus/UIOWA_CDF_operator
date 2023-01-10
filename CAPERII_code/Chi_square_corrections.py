# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- CAPER Chi Square correction factors --------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# --------------------
# IDEA
# --------------------

# Goal : There may be pads that are a bit too high or two low on the TRICE data set and must be corrected via chi-square minimization to fix this

# Method:
# [1] We first look for a pair of detector pads that overlay in their actual pitch angle, namely pitch angle with fall within 5deg of one paid (e.g. pad 40deg is at 47deg true, and pad 50deg at 49 deg true)
#       -> The logic here is that two pads looking at approximately the same true pitch angle should have about the same flux value, if they don't its due to the pad calibration not the physics
# [2] For any one time slice, sweep through in energy the two pads of interest and store necesary data (counts, true pitch angle, and index identifiers).
# [4] The data is now paired; one count value & pitch for one pad also one count value and pitch for the second pad. One of these pads is the principal look angles, the other is the uncal pad.
# [5] Continue process [4] through all time slices and collect all the pairs between the two pads
# [6] You will now preform the chi-square minimization technique given by \chi ~ 1/(N-1) \sum \frac{data_{principal} - \alpha * data_{uncal}   }{ error_{principal} + error_{uncal}}
# [7] Create a method to sweep through possible values for the alpha fitting parameter such that \chi becomes ~1
# [8] Repeat steps [1] - [7] for the next pad pairs
# [9] Likely I will move from low pads to high pads, do the results become teh same if I perform them in reverse?



import numpy as np, itertools,time,math
from cdflib import cdfread

from Variables import counts_2,counts_1,pitch_angle_actual_1

print('importing variables: ', end='')
start_time = time.time()

print('Done')








# -----------------------------------------------------
# TEST CASE - Pad 90deg (10) principal, Pad 80deg uncal
# -----------------------------------------------------

prin = 10
uncal = prin - 1







# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



print("--- %s seconds for Sorting_DiffFlux---" % (time.time() - start_time) ,'\n')