# Example of pcolormesh + colorbar
import itertools

import numpy as np

from ACESII_code.myImports import *

Xdat = [1,2,3,4]
Ydat = [5,6,7]

print(np.histogram_bin_edges(Xdat,bins=5))
