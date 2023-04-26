# # --- scratchPaper.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering
import numpy as np
from ACESII_code.class_var_func import Ry
n = np.array([1,0,0])

print(np.matmul(Ry(90),n))