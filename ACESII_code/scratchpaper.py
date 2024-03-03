
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0,2*np.pi,40)
y = np.sin(x)

negativeIndicies = np.where(y<0)[0]
positiveIndicies = np.where(y>0)[0]


fig, ax = plt.subplots()
ax.plot(x,y)
ax.fill_between(x[positiveIndicies], y[positiveIndicies],color='red',alpha=0.4)
ax.fill_between(x[negativeIndicies], y[negativeIndicies],color='blue',alpha=0.4)
plt.gca()
plt.show()