# Example of pcolormesh + colorbar
import itertools

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *

a = np.array([[1,2,3],[5,6,7]])

rowIdx,colIdx = np.where(a == a.max())

print(a[rowIdx[0]][colIdx[0]])