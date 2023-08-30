# # --- scratchPaper.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering

from ACESII_code.class_var_func import loadDictFromFile, DCM
import matplotlib.pyplot as plt
import numpy as np
import itertools



# Load the data
data_dict_IGRF = loadDictFromFile(r'C:\Data\ACESII\science\IGRF_interpolated\high\ACESII_36359_IGRF_interpolated_model.cdf',{})
data_dict_Brkt = loadDictFromFile(r'C:\Data\ACESII\L2\high\ACESII_36359_l2_RingCore_Despun.cdf',{})


# construct the raw data vectors
B_IGRF = np.array([[data_dict_IGRF['IGRF_X'][0][i],data_dict_IGRF['IGRF_Y'][0][i],data_dict_IGRF['IGRF_Z'][0][i]] for i in range(len(data_dict_IGRF['Epoch'][0]))])
B_rkt = np.array([[data_dict_Brkt['B_x'][0][i],data_dict_Brkt['B_y'][0][i],data_dict_Brkt['B_z'][0][i]] for i in range(len(data_dict_Brkt['Epoch'][0]))])
Epoch_rkt = np.array(data_dict_Brkt['Epoch'][0])
Epoch_IGRF = np.array(data_dict_IGRF['Epoch'][0])



fig, ax = plt.subplots(3)

# Initialize data
Bx_rkt, = ax[0].plot(Epoch_rkt, B_rkt[:, 0])
Bx_IGRF, = ax[0].plot(Epoch_IGRF, B_IGRF[:, 0])
By_rkt, = ax[1].plot(Epoch_rkt, B_rkt[:, 1])
By_IGRF, = ax[1].plot(Epoch_IGRF, B_IGRF[:, 1])
Bz_rkt, = ax[2].plot(Epoch_rkt, B_rkt[:, 2])
Bz_IGRF, = ax[2].plot(Epoch_IGRF, B_IGRF[:, 2])

# set the labels
ax[0].set_ylabel('$B_{X}$')
ax[1].set_ylabel('$B_{Y}$')
ax[2].set_ylabel('$B_{Z}$')
ax[2].set_xlabel('Epoch')

# set title
roll = 0
pitch = 0
yaw = 0
fig.suptitle(f'Roll: {roll}, Pitch: {pitch}, Yaw: {yaw}')

# draw and show it
fig.canvas.draw()
plt.show(block=False)

# loop to update the data
angleUpperLims = [5,5,5] # roll, pitch,yaw
N = 10
angleRanges = [np.linspace(0,angleUpperLims[0],N),np.linspace(0,angleUpperLims[1],N),np.linspace(0,angleUpperLims[2],N) ]

for i,j,k in itertools.product(*angleRanges):
    try:
        roll = 0 + i
        pitch = 0 + j
        yaw = 0 + k

        # change the data
        B_rkt_rot = np.array([ np.matmul(DCM(roll,pitch,yaw),B_rkt[k]) for k in range(len(Epoch_rkt))])

        # set the new data
        Bx_rkt.set_ydata(B_rkt_rot[:, 0])
        By_rkt.set_ydata(B_rkt_rot[:, 1])
        Bz_rkt.set_ydata(B_rkt_rot[:, 2])

        # update axes?
        # for i in range(3):
        # ax.relim()
        # ax.autoscale_view(True, True, True)

        # update title
        fig.suptitle(f'Roll: {roll}, Pitch: {pitch}, Yaw: {yaw}')

        fig.canvas.draw()
        plt.pause(0.005)

    except KeyboardInterrupt:
        plt.close('all')
        break