# # --- scratchpaper.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering

import matplotlib.pyplot as plt
import numpy as np

a = np.array([1,-1,2,3])

print(a>0)

if all(a>0):
    print('ALl are true')
else:
    print('POTATO')

