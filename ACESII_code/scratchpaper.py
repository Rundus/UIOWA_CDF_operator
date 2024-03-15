
from ACESII_code.class_var_func import u0
import numpy as np
import matplotlib.pyplot as plt


# Plot reflection function
fig, ax = plt.subplots()
freq = 2*np.pi*np.linspace(0,64,100)
VA_quiet= 963411.0134456921
VA_mid = 620542.3836004707
def IAR(SimgaP,freq, Alt):
    VA = 1228498.9582697104
    SigmaA = 1 /(u0*VA)
    R = (SigmaA - SimgaP)/(SigmaA + SimgaP)
    print(R)
    return VA * (1 - R*np.cos(-2*freq*Alt/VA)) / (1+R*np.cos(-2*freq*Alt/VA))

IAR_Ratio = [IAR(1,val,60000) for val in freq]
ax.plot(freq,IAR_Ratio, color='orange')
ax.set_yscale('log')
ax.set_ylim(1E5,1E7)
ax.set_xlim(0,10)
plt.show()