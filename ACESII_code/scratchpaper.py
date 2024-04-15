from ACESII_code.myImports import *

import numpy as np
from scipy import signal
from scipy.fft import fftshift
import matplotlib.pyplot as plt

a = np.array([1.001,2,3])
b = np.array([-2.77E7,-9.71E6,-8.29E6])

def fitFunc(x,A,B,C):
    return A*np.log(x-B) + C

guess= [1E6,1,-1E7]

params, cov = curve_fit(fitFunc,a,b,p0=guess,maxfev=1000000)

def invFunc(x,A,B,C):
    return np.exp((x-C)/A) + B

mySlopes = np.array([-9.156E6,-5.417E6,-2.919E6,-9.86E6])
print(invFunc(mySlopes,*params))

xDataFit = np.linspace(1,10,10000)
yDataFit = fitFunc(xDataFit,*params)
fig, ax = plt.subplots()
ax.plot(a,b,color='blue')
ax.plot(xDataFit,yDataFit,color='red')
plt.show()