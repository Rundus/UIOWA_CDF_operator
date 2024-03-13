
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

x,y,c = zip(*np.random.rand(30,3)*4-2)

norm=plt.Normalize(0,1)
colors = ['m',"red",'orange','yellow','green','blue','black']
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors[::-1])

plt.scatter(x,y,c=c, cmap=cmap, norm=norm)
plt.colorbar()
plt.show()