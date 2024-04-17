# Example of pcolormesh + colorbar
import itertools

from ACESII_code.myImports import *

Xdat = [1,2,3,4]
Ydat = [5,6,7]

Z = np.array([[1,2,3,4],
     [2,3,4,5],
     [3,4,5,6]])


tempSet = list(set(list(itertools.chain(*[list(set(ar)) for ar in Z]))))

print(tempSet)


