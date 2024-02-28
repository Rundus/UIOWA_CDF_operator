import numpy as np
simLen = 331 # how many delta T steps to simulate
deltaT = 0.01 # in seconds
simTime =np.linspace(0,simLen*deltaT,simLen+1)
print(len(simTime))
print(simTime)