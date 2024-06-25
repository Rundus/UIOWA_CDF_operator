import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *

EnergyBins = [13678.4, 11719.21, 10040.64, 8602.5, 7370.34, 6314.67, 5410.2, 4635.29,
                 3971.37, 3402.54, 2915.18, 2497.64, 2139.89, 1833.39, 1570.79, 1345.8,
                 1153.04, 987.89, 846.39, 725.16, 621.29, 532.3, 456.06, 390.74,
                 334.77, 286.82, 245.74, 210.54, 180.39, 154.55, 132.41, 113.45,
                 97.2, 83.28, 71.35, 61.13, 52.37, 44.87, 38.44, 32.94,
                 28.22, 24.18, 20.71, 17.75, 15.21, 13.03, 11.16, 9.56, 8.19]

energyResolution = 0.18
wBin = 20
sigma = energyResolution*EnergyBins[wBin] / (2*np.sqrt(2*np.log(2)))
mean = EnergyBins[wBin]
def normalDistribution(x,mean,sigma):
    return np.exp(- np.power((x-mean),2)/(2*np.power(sigma,2)))

print(EnergyBins[wBin],sigma)
fig, ax = plt.subplots()
xData = np.linspace(-1E3,1E3,1000)
yData = normalDistribution(xData,mean,sigma)
yData_Norm = yData/np.max(yData)
ax.plot(xData,yData)
plt.show()