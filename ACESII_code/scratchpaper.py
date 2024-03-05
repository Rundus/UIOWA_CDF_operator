
import matplotlib.pyplot as plt
import numpy as np

a = np.array([ [[1,2,3],[1,2,3]],[[1,2,3],[1,2,3]],[[1,2,3],[1,2,3]], [[1,2,3],[1,2,3]],[[1,2,3],[1,2,3]],[[1,2,3],[1,2,3]]])
b = np.array([ [[4,5,6],[4,5,6]],[[4,5,6],[4,5,6]],[[4,5,6],[4,5,6]],[[4,5,6],[4,5,6]],[[4,5,6],[4,5,6]],[[4,5,6],[4,5,6]]])

print(np.append(a,b,axis=2))


# x = np.linspace(0,2*np.pi,40)
# y = np.sin(x)
#
# negativeIndicies = np.where(y<0)[0]
# positiveIndicies = np.where(y>0)[0]
#
#
# fig, ax = plt.subplots()
# ax.plot(x,y)
# ax.fill_between(x[positiveIndicies], y[positiveIndicies],color='red',alpha=0.4)
# ax.fill_between(x[negativeIndicies], y[negativeIndicies],color='blue',alpha=0.4)
# plt.gca()
# plt.show()