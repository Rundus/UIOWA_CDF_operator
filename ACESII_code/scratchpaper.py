import numpy as np

from ACESII_code.myImports import *

a = np.array([[5,2,3],[0,0,3],[3,5,3]]).T

print(np.where(np.abs(a) < 2))
finder = np.where(np.abs(a) < 4)
a[finder] = 0
print(a)

