import numpy as np


a = np.array([[1,2,3],[4,5,6],[7,8,9]])
b  = a + 10
b[np.where(10>a >5)] = 0
print(b)