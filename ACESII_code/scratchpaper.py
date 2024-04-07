from ACESII_code.myImports import *

import numpy as np
from scipy import signal
from scipy.fft import fftshift
import matplotlib.pyplot as plt

phaseDiffBins = [-187.5 + 15 * i for i in range(26)]
plotBins = [-180 + 15 * i for i in range(25)]

print(phaseDiffBins)
print(plotBins)