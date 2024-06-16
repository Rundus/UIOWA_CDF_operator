# --- plot.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to store some plotting functions I use

# Imports
import matplotlib.pyplot as plt


def plot3Axis(Epoch, AxesData, ylabels):

    fig, ax = plt.subplots(3)

    for i in range(3):
        ax[i].plot(Epoch,AxesData[:, i])
        ax[i].set_ylabel(ylabels[i])

    ax[2].set_xlabel('Epoch')
    plt.show()