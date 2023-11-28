# # --- scratchpaper.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering



import matplotlib.pyplot as plt
import numpy as np


a = [66 + i for i in range(8)]
b = [100 + i*50 for i in range(8)]

fig, ax = plt.subplots()

slope = -1*(111/np.sin(np.radians(90 - 78.13))) # corresponds to line with -78.13deg inclination
for i in range(21):
    ax.axline(xy1=(66+i*0.5, 0), slope=slope,color='black', linewidth=2, linestyle='-.', alpha=0.2)

plt.plot(a,b,color='red')
plt.ylim(0, 500)
plt.xlim(66.5, 73.5)
plt.show()