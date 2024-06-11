import numpy as np

from ACESII_code.myImports import *

a,b = np.meshgrid([1,2,3],[5,6,7])
print(a.flatten())
print(b.flatten())