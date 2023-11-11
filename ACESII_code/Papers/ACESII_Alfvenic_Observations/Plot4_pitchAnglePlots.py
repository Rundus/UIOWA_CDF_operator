# --- Plots4_pitchAnglePlots.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Recreation of DiffEFlux pitch angle Plots, focusing on
# a few particle signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from itertools import product

print(color.UNDERLINE + f'Plot4_pitchAnglePlots' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
reduceData =  False
targetTimes = [dt.datetime(2022, 11, 20, 17, 24, 55, 500000), dt.datetime(2022, 11, 20, 17, 25, 9, 500000)]
wEngyLim = 20 # Energy index used to determine limits of x-y plot when converted into V_perp,V_parallel
cbarLow,cbarHigh = [1E6,5E8]

# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---

# Trajectory data (for geomagnetic lat/long info)
inputTracj_low = glob(r'C:\Data\ACESII\trajectories\low\*ILat_ILong*')
data_dict_traj_low = loadDictFromFile(inputFilePath=inputTracj_low[0],input_data_dict={},reduceData=reduceData,targetTimes=targetTimes,wKeys=[])
inputTracj_high = glob(r'C:\Data\ACESII\trajectories\high\*ILat_ILong*')
data_dict_traj_high = loadDictFromFile(inputFilePath=inputTracj_high[0],input_data_dict={},reduceData=reduceData,targetTimes=targetTimes,wKeys=[])

# EEPAA Particle Data
inputEEPAA_low = glob('C:\Data\ACESII\L2\low\*eepaa_fullCal*')
data_dict_eepaa_low = loadDictFromFile(inputFilePath=inputEEPAA_low[0],input_data_dict={},reduceData=reduceData,targetTimes=targetTimes,wKeys=[])
inputEEPAA_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')
data_dict_eepaa_high = loadDictFromFile(inputFilePath=inputEEPAA_high[0],input_data_dict={},reduceData=reduceData,targetTimes=targetTimes,wKeys=[])

prgMsg('Introducing fillvals for log plot')
rocketAttrs,b,c = ACES_mission_dicts()
ranges = [range(len(data_dict_eepaa_high['Epoch'][0])),range(len(data_dict_eepaa_high['Pitch_Angle'][0])), range(len(data_dict_eepaa_high['Energy'][0]))]
for tme, ptch, engy in product(*ranges):
    if data_dict_eepaa_high['Differential_Energy_Flux'][0][tme][ptch][engy] == 0:
        data_dict_eepaa_high['Differential_Energy_Flux'][0][tme][ptch][engy] = 1
    elif data_dict_eepaa_high['Differential_Energy_Flux'][0][tme][ptch][engy] == rocketAttrs.epoch_fillVal:
        data_dict_eepaa_high['Differential_Energy_Flux'][0][tme][ptch][engy] = 1



# Calculate Vperp and Vparallel
Energy = data_dict_eepaa_high['Energy'][0]
Pitch = data_dict_eepaa_high['Pitch_Angle'][0]
EnergyLimit = Energy[wEngyLim]
diffEFlux = deepcopy(data_dict_eepaa_high['Differential_Energy_Flux'][0][6000])
Vperp = deepcopy(diffEFlux)
Vpara = deepcopy(diffEFlux)

for ptch in range(len(diffEFlux)):
    for engy in range(len(diffEFlux[0])):
        Emag = np.sqrt(2*q0*Energy[engy]/(m_e))
        Vperp[ptch][engy] = np.sin(np.radians(Pitch[ptch]))*Emag
        Vpara[ptch][engy] = np.cos(np.radians(Pitch[ptch]))*Emag

Vpara, Vperp = np.array(Vpara)/1000,np.array(Vperp)/1000
VelLimit = np.sqrt(2 *EnergyLimit*q0/(m_e))/1000


# get the time titles of the plots
timesIndicies = [[5963, 5969, 5974],
                 [5986, 5990, 5994], # scale x: -500 to 500, scale y: -100 to 500
                 [6014,6017,6023]]

Done(start_time)

for j in range(3):
    fig, ax = plt.subplots()
    colors=['red','blue','green']

    for i in range(3):
        AvFLuxNorm = np.array([np.average(data_dict_eepaa_high['Differential_Energy_Flux'][0][timesIndicies[j][i]][ptch]) for ptch in range(len(Pitch))])
        AvFLuxNorm=AvFLuxNorm[1:20]
        xData = np.array(Pitch[1:20])

        badIndicies = []
        for k in range(len(xData)):
            if AvFLuxNorm[k] < 1000:
                badIndicies.append(k)

        AvFLuxNorm = np.delete(AvFLuxNorm,badIndicies)/AvFLuxNorm.max()
        xData = np.delete(xData,badIndicies)

        ax.set_title(data_dict_eepaa_high['Epoch'][0][timesIndicies[j][i]])
        ax.plot(xData, AvFLuxNorm,color=colors[i])
        ax.set_xticks(Pitch[1:20])

    plt.show()




############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

prgMsg('Beginning Plot')


fig = plt.figure()
figure_height = 15
figure_width = 10
fig.set_figwidth(figure_width)
fig.set_figheight(figure_height)

# title
# fig.suptitle('ACESII')

##### Define the primary Gridspec #####
topPlotHeight = 0.15
gs0 = gridspec.GridSpec(nrows=4, ncols=4, figure=fig, width_ratios=[0.325, 0.325, 0.325, 0.025],height_ratios=[topPlotHeight, (1-topPlotHeight)/3, (1-topPlotHeight)/3, (1-topPlotHeight)/3], wspace=0.25, hspace=0.3) # splits figure between the plots and the colorbar at the very bottom


# --- ALFVENIC SIGNAUTRE PLOT ---
axAlfSig = fig.add_subplot(gs0[0, 0:3])
axAlfSig.set_title('ACESII 36359')
IndexLow,IndexHigh = np.abs(data_dict_eepaa_high['Epoch'][0] - targetTimes[0]).argmin(),np.abs(data_dict_eepaa_high['Epoch'][0] - targetTimes[1]).argmin()
Epoch = data_dict_eepaa_high['Epoch'][0][IndexLow:IndexHigh]
diffEFlux = np.array(data_dict_eepaa_high['Differential_Energy_Flux'][0][IndexLow:IndexHigh])
diffEFlux = np.transpose(diffEFlux[:,2,:])
alfSigPlot = axAlfSig.pcolormesh(Epoch, Energy, diffEFlux, cmap='turbo', shading='nearest', norm=matplotlib.colors.LogNorm(vmin=cbarLow,vmax=cbarHigh))
Eflip = Energy[::-1]
axAlfSig.set_ylim(Energy[-1],Eflip[22])
axAlfSig.set_ylabel('Energy [eV]')
# axAlfSig.set_yscale('log')




# --- Pitch Angle PLOTS ---
times = [ [dt.datetime(2022, 11, 20, 17, 24, 58, 212000), dt.datetime(2022, 11, 20, 17, 24, 58, 462000), dt.datetime(2022, 11, 20, 17, 24, 58, 712000)], # s2
          [dt.datetime(2022, 11, 20, 17, 24, 59, 312000), dt.datetime(2022, 11, 20, 17, 24, 59, 512000), dt.datetime(2022, 11, 20, 17, 24, 59, 712000)], # s3
          [dt.datetime(2022, 11, 20, 17, 25, 00, 712000), dt.datetime(2022, 11, 20, 17, 25,  0, 860000), dt.datetime(2022, 11, 20, 17, 25,  1, 162000)], # s5 + Inverted V
          ]

timesIndicies = [[5963, 5969, 5974],
                 [5986, 5990, 5994], # scale x: -500 to 500, scale y: -100 to 500
                 [6014,6017,6023]]




for j in range(len(times)): # row
    for i in range(len(times[0])): # column
        axS = fig.add_subplot(gs0[j+1,i])
        axS.set_title( 't='+ times[j][i].strftime('%H:%M:%S.%f'),fontsize=14, weight='bold')

        diffEFlux = data_dict_eepaa_high['Differential_Energy_Flux'][0][np.abs(data_dict_eepaa_high['Epoch'][0] - times[j][i]).argmin()]

        axS.pcolormesh(Vperp, Vpara, diffEFlux, cmap='turbo', shading='nearest',norm=matplotlib.colors.LogNorm(vmin=cbarLow,vmax=cbarHigh))
        axS.set_xlim(-5000, VelLimit)
        axS.set_ylim(-VelLimit, VelLimit)
        if i == 0:
            axS.set_ylabel('V$_{\parallel}$ [km/s]',fontsize=14)
        else:
            axS.set_yticklabels([])

        if j ==2:
            axS.set_xlabel('V$_{\perp}$ [km/s]',fontsize=14)
        else:
            axS.set_xticklabels([])

        # add in the average flux plot




        # axS.set_aspect('equal')
        axS.invert_yaxis()

Done(start_time)

# --- Colorbar Column - LF ---
axCbar = fig.add_subplot(gs0[0:4, 3])
cbar = fig.colorbar(mappable=alfSigPlot, cax=axCbar, norm='log')
cbar.set_label('Differential Energy Flux [eV/cm$^{2}$-s-sr-eV]',fontsize=16)
for l in cbar.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(16)




# --- show plot ---
plt.tight_layout()
# plt.show()

plt.savefig(r'C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\Papers\ACESII_Alfvenic_Observations\Plots\\Plot4_pitchAngle.png')