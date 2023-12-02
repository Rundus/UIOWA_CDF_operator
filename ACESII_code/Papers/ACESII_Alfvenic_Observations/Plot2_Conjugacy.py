# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing field-aligned
# particle data along with electric and magnetic signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

print(color.UNDERLINE + f'Plot2_Conjugacy' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# targetInvariantlats = [70.15, 70.5]
targetInvariantlats = [69.8, 70.8]
# cbarLimits = [1E7,1E9]
cbarLimits = [4E5,5E8]
wPitch_engy_val = 38 # the chosen energy (0 - 41) for thePitch Angle Plot

# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---
prgMsg('Loading Data')

# delta B
inputMagFiles_high = glob('C:\Data\ACESII\L3\deltaB\high\*Field_Aligned*')
data_dict_mag_high = loadDictFromFile(inputFilePath=inputMagFiles_high[0],input_data_dict={},reduceData=False,targetTimes=[],wKeys=[])
inputMagFiles_low = glob('C:\Data\ACESII\L3\deltaB\low\*Field_Aligned*')
data_dict_mag_low = loadDictFromFile(inputFilePath=inputMagFiles_low[0],input_data_dict={},reduceData=False,targetTimes=[],wKeys=[])

# delta E
inputEFIFiles_low = glob('C:\Data\ACESII\L3\deltaE\low\*Field_Aligned*')
data_dict_elec_low = loadDictFromFile(inputFilePath=inputEFIFiles_low[0],input_data_dict={},reduceData=False,targetTimes=[],wKeys=[])

# EEPAA Particle Data
inputEEPAA_low = glob('C:\Data\ACESII\L2\low\*eepaa_fullCal*')
data_dict_eepaa_low = loadDictFromFile(inputFilePath=inputEEPAA_low[0],input_data_dict={},reduceData=False,targetTimes=[],wKeys=[])
inputEEPAA_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')
data_dict_eepaa_high = loadDictFromFile(inputFilePath=inputEEPAA_high[0],input_data_dict={},reduceData=False,targetTimes=[],wKeys=[])


# Trajectory data (for geomagnetic lat/long info)
inputTracj_low = glob(r'C:\Data\ACESII\trajectories\low\*ILat_ILong*')
data_dict_traj_low = loadDictFromFile(inputFilePath=inputTracj_low[0],input_data_dict={},reduceData=False,targetTimes=[],wKeys=[])
inputTracj_high = glob(r'C:\Data\ACESII\trajectories\high\*ILat_ILong*')
data_dict_traj_high = loadDictFromFile(inputFilePath=inputTracj_high[0],input_data_dict={},reduceData=False,targetTimes=[],wKeys=[])

# Poynting Data
inputPoynting_low = glob(r'C:\Data\ACESII\science\PoyntingFlux\low\*Field_Aligned*')
data_dict_poynting_low = loadDictFromFile(inputFilePath=inputPoynting_low[0],input_data_dict={},reduceData=False,targetTimes=[],wKeys=[])
inputPoynting_high = glob(r'C:\Data\ACESII\science\PoyntingFlux\high\*Field_Aligned*')
data_dict_poynting_high = loadDictFromFile(inputFilePath=inputPoynting_high[0],input_data_dict={},reduceData=False,targetTimes=[],wKeys=[])
Done(start_time)

# --- --- --- --- --- --- --- --- --- -
# --- Calculate Invariant Lattitude ---
# --- --- --- --- --- --- --- --- --- -

# determine invariant lattitude for the particle
invariantLat_ptcle_low = np.array([np.degrees(np.arccos( np.cos(np.radians(data_dict_traj_low['geomagLat'][0][i]))/np.sqrt(data_dict_traj_low['geomagAlt'][0][i]) )) for i in range(len(data_dict_traj_low['Epoch_esa'][0]))])
invariantLat_ptcle_high = np.array([np.degrees(np.arccos( np.cos(np.radians(data_dict_traj_high['geomagLat'][0][i]))/np.sqrt(data_dict_traj_high['geomagAlt'][0][i]) )) for i in range(len(data_dict_traj_high['Epoch_esa'][0]))])

# determine invariant lattitude for the Wave Data
invariantLat_wave_low = np.array([np.degrees(np.arccos( np.cos(np.radians(data_dict_poynting_low['Lat_geom'][0][i]))/np.sqrt(data_dict_poynting_low['Alt_geom'][0][i]) )) for i in range(len(data_dict_poynting_low['Epoch'][0]))])
invariantLat_wave_high = np.array([np.degrees(np.arccos( np.cos(np.radians(data_dict_poynting_high['Lat_geom'][0][i]))/np.sqrt(data_dict_poynting_high['Alt_geom'][0][i]) ))for i in range(len(data_dict_poynting_high['Epoch'][0]))])
Done(start_time)


# --- --- --- --- --- --- --- --- --- --- ---
# --- REDUCE DATA BASED ON Invariant LAT ---
# --- --- --- --- --- --- --- --- --- --- ---
prgMsg('Reducing Data')

def reduceDictonary(inputDict, InvariantLat):
    lowCut,highCut = np.abs(InvariantLat - targetInvariantlats[0]).argmin(),np.abs(InvariantLat - targetInvariantlats[1]).argmin()
    for key,val in inputDict.items():
        inputDict[key][0] = inputDict[key][0][lowCut:highCut]
    return inputDict

data_dict_mag_low = reduceDictonary(data_dict_mag_low, invariantLat_wave_low)
data_dict_mag_high = reduceDictonary(data_dict_mag_high, invariantLat_wave_high)
data_dict_elec_low = reduceDictonary(data_dict_elec_low, invariantLat_wave_low)
data_dict_poynting_low = reduceDictonary(data_dict_poynting_low, invariantLat_wave_low)
data_dict_poynting_high = reduceDictonary(data_dict_poynting_high, invariantLat_wave_high)
data_dict_traj_low = reduceDictonary(data_dict_traj_low, invariantLat_ptcle_low)
data_dict_traj_high = reduceDictonary(data_dict_traj_high, invariantLat_ptcle_high)

# special case for EEPAA data
lowCut, highCut = np.abs(invariantLat_ptcle_high - targetInvariantlats[0]).argmin(), np.abs(invariantLat_ptcle_high - targetInvariantlats[1]).argmin()
for key in ['Differential_Energy_Flux', 'Epoch']:
    data_dict_eepaa_high[key][0] = data_dict_eepaa_high[key][0][lowCut:highCut]

lowCut, highCut = np.abs(invariantLat_ptcle_low - targetInvariantlats[0]).argmin(), np.abs(invariantLat_ptcle_low - targetInvariantlats[1]).argmin()
for key in ['Differential_Energy_Flux', 'Epoch']:
    data_dict_eepaa_low[key][0] = data_dict_eepaa_low[key][0][lowCut:highCut]

# define some variables
Energy = data_dict_eepaa_high['Energy'][0]
Pitch = data_dict_eepaa_high['Pitch_Angle'][0]

# reduce the invariant lattitudes themselves
lowCut, highCut = np.abs(invariantLat_ptcle_low - targetInvariantlats[0]).argmin(), np.abs(invariantLat_ptcle_low - targetInvariantlats[1]).argmin()
invariantLat_ptcle_low =invariantLat_ptcle_low[lowCut:highCut]
lowCut, highCut = np.abs(invariantLat_ptcle_high - targetInvariantlats[0]).argmin(), np.abs(invariantLat_ptcle_high - targetInvariantlats[1]).argmin()
invariantLat_ptcle_high =invariantLat_ptcle_high[lowCut:highCut]
lowCut, highCut = np.abs(invariantLat_wave_low - targetInvariantlats[0]).argmin(), np.abs(invariantLat_wave_low - targetInvariantlats[1]).argmin()
invariantLat_wave_low =invariantLat_wave_low[lowCut:highCut]
lowCut, highCut = np.abs(invariantLat_wave_high - targetInvariantlats[0]).argmin(), np.abs(invariantLat_wave_high - targetInvariantlats[1]).argmin()
invariantLat_wave_high =invariantLat_wave_high[lowCut:highCut]


# --- --- --- --- --- --- --
# --- GET THE XTICK INFO ---
# --- --- --- --- --- --- --

# description: the poynting flux data is the last databox. Its' ticks
# define the x-axis labels for all the data. We need to collect its ticks for display

N = 6 # number of ticks
datLength = len(data_dict_poynting_low['S_p'][0])
wIndices_low = [i for i in range(0, datLength, int(len(data_dict_poynting_low['S_p'][0])/N))]
UTTicks_low = [data_dict_poynting_low['Epoch'][0][i].strftime('%H:%M:%S') for i in wIndices_low]

datLength = len(data_dict_poynting_high['S_p'][0])
wIndices_high = [i for i in range(0, datLength, int(len(data_dict_poynting_high['S_p'][0])/N))]
UTTicks_high = [data_dict_poynting_high['Epoch'][0][i].strftime('%H:%M:%S') for i in wIndices_high]




# E_perp_low = np.array([ np. for i in range(len(data_dict_elec_low['Epoch'][0]))])



############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

fig = plt.figure()
figure_height = 10
figure_width = 24
fig.set_figwidth(figure_width)
fig.set_figheight(figure_height)

# title
# fig.suptitle('ACESII')

##### Define the primary Gridspec #####
gs0 = gridspec.GridSpec(nrows=4, ncols=4, figure=fig, width_ratios=[0.49, 0.01, 0.49, 0.01], hspace=0.1,wspace=0.2) # splits figure between the plots and the colorbar at the very bottom


#%%%%%%%%%% HF %%%%%%%%%%%
# --- First Data Column ---

# HF eepaa data
axEEPAA_HF = fig.add_subplot(gs0[0, 0])
axEEPAA_HF.get_xaxis().set_visible(False)
xData, yData = invariantLat_ptcle_high, Energy
zData = np.transpose(data_dict_eepaa_high['Differential_Energy_Flux'][0][:, 2, :])
cmap_EEPAA_HF = axEEPAA_HF.pcolormesh(xData, yData, zData, vmin=cbarLimits[0], vmax=cbarLimits[1], shading='nearest',cmap='turbo',norm='log')
axEEPAA_HF.set_title('ACESII 36359\n'+r'$\alpha=10^{\circ}$')
axEEPAA_HF.set_ylabel('Energy [eV]')
axEEPAA_HF.set_yscale('log')
axEEPAA_HF.set_ylim(30, 1500)


# HF eepaa PITCH data
axEEPAA_pitch_HF = fig.add_subplot(gs0[1,0])
axEEPAA_pitch_HF.get_xaxis().set_visible(False)
xData, yData = invariantLat_ptcle_high, Pitch
zData = np.transpose(data_dict_eepaa_high['Differential_Energy_Flux'][0][:, :, wPitch_engy_val])
cmap_EEPAA_pitch_HF = axEEPAA_pitch_HF.pcolormesh(xData, yData, zData, vmin=cbarLimits[0], vmax=cbarLimits[1], shading='nearest',cmap='turbo',norm='log')
axEEPAA_pitch_HF.set_ylabel('Pitch Angle [$^{\circ}$]\n'+f' Energy ={Energy[32]}eV')
axEEPAA_pitch_HF.set_ylim(0, 180)


# HF E data
axEB_HF_E = fig.add_subplot(gs0[2,0])
axEB_HF_E.get_xaxis().set_visible(False)
axEB_HF_E.annotate('Not Available',xy=(0.5,0.5),xytext=(0.0,0.0),ha='center',va='center',xycoords='axes fraction',textcoords='offset points',fontsize=20,color='tab:red')
axEB_HF_E.set_ylabel('$\delta E$r [mV/m]',color='tab:red')
axEB_HF_E.set_ylim(-12,12)
axEB_HF_E.spines['left'].set_color('tab:red')
axEB_HF_E.xaxis.label.set_color('tab:red')
axEB_HF_E.tick_params(axis='y',colors='tab:red')
axEB_HF_E.margins(x=0)

# HF B data
axEB_HF_B = fig.add_subplot(gs0[3,0])
axEB_HF_B.get_xaxis().set_visible(False)
axEB_HF_B.plot(invariantLat_wave_high, data_dict_mag_high['B_e'][0],color='tab:blue')
axEB_HF_B.set_ylabel('$\delta$Be [nT]',color='tab:blue')
axEB_HF_B.spines['left'].set_color('tab:blue')
axEB_HF_B.xaxis.label.set_color('tab:blue')
axEB_HF_B.tick_params(axis='y',colors='tab:blue')
axEB_HF_B.margins(x=0)

# # HF Poynting Flux data
# axS_HF = fig.add_subplot(gs0[4,0])
# axS_HF.plot(invariantLat_wave_high,(1E3)*data_dict_poynting_high['S_p'][0],label='Sp',color='black')
# axS_HF.set_ylabel('$\delta B_{\perp}^{2} V_{A}/2\mu_{0}$ \n [ergs/cm$^{2}$s]')
# axS_HF.set_ylim(-0.25,0.25)
# axS_HF.annotate('ILat [$^{\circ}$]',xy=(-0.15,-0.04),xytext=(0,0),ha='left',va='top',xycoords='axes fraction',textcoords='offset points')
# axS_HF.xaxis.set_major_locator(ticker.LinearLocator(N))
# axS_HF.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f')) # rounds the ILat axis to 2 decimal places
# axS_HF.margins(x=0)

#add additional x axes(ALT)
axS_HF.annotate('Alt [km]',xy=(-0.15,-0.16),xytext=(0,0),ha='left',va='top',xycoords='axes fraction',textcoords='offset points')
axS_HF_Alt = axS_HF.twiny() # create second axis
axS_HF_Alt.spines['bottom'].set_position(('axes',-0.14)) # move the axes to the bottom
axS_HF_Alt.xaxis.set_ticks_position("bottom") # move the ticks to the bottom
axS_HF_Alt.set_frame_on(True)
axS_HF_Alt.patch.set_visible(False)
axS_HF_Alt.spines['bottom'].set_visible(False)
axS_HF_Alt.tick_params(axis='x', length=0)
axS_HF_Alt.plot(data_dict_traj_high['geoAlt'][0], data_dict_traj_high['geomagAlt'][0], alpha=0)
axS_HF_Alt.xaxis.set_major_locator(ticker.LinearLocator(N))
axS_HF_Alt.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f')) # rounds the ILat axis to 2 decimal places


#add additional x axes(UT Time)
axS_HF.annotate('UT ', xy=(-0.15,-0.28), xytext=(0,0), ha='left', va='top', xycoords='axes fraction',textcoords='offset points')
axS_HF_UT = axS_HF.twiny() # create second axis
axS_HF_UT.spines['bottom'].set_position(('axes',-0.26)) # move the axes to the bottom
axS_HF_UT.xaxis.set_ticks_position("bottom") # move the ticks to the bottom
axS_HF_UT.set_frame_on(True)
axS_HF_UT.patch.set_visible(False)
axS_HF_UT.spines['bottom'].set_visible(False)
axS_HF_UT.tick_params(axis='x',length=0)
axS_HF_UT.plot(data_dict_traj_high['Epoch_esa'][0], data_dict_traj_high['geomagAlt'][0], alpha=0)
axS_HF_UT.xaxis.set_major_locator(ticker.LinearLocator(N))
axS_HF_UT.set_xticklabels(UTTicks_high)


# --- 1st Colorbar Column - HF ---
axCbar_HF = fig.add_subplot(gs0[0,1])
cbar_HF = fig.colorbar(mappable=cmap_EEPAA_HF,
                       cax=axCbar_HF,
                       norm='log')
cbar_HF.set_label('eV/cm$^{2}$-s-sr-eV',labelpad=-52)
l, b, w, h = axCbar_HF.get_position().bounds
axCbar_HF.set_position([(0.955)*l, b, w, h])

# --- 2nd Colorbar Column - HF ---
axCbar_pitch_HF = fig.add_subplot(gs0[1,1])
cbar_pitch_HF = fig.colorbar(mappable=cmap_EEPAA_pitch_HF,
                       cax=axCbar_pitch_HF,
                       norm='log')
cbar_pitch_HF.set_label('eV/cm$^{2}$-s-sr-eV',labelpad=-52)
l, b, w, h = axCbar_pitch_HF.get_position().bounds
axCbar_pitch_HF.set_position([(0.955)*l, b, w, h])


#%%%%%%%%%% LF %%%%%%%%%%%
# --- 2nd Data Column ---

# LF EEPAA Data
axEEPAA_LF = fig.add_subplot(gs0[0, 2])
axEEPAA_LF.get_xaxis().set_visible(False)
xData, yData = invariantLat_ptcle_low, Energy
zData = np.transpose(data_dict_eepaa_low['Differential_Energy_Flux'][0][:, 2, :])
cmap_EEPAA_LF = axEEPAA_LF.pcolormesh(xData, yData, zData, vmin=cbarLimits[0], vmax=cbarLimits[1], shading='nearest',cmap='turbo',norm='log')
axEEPAA_LF.set_title('ACESII 36364\n'+ r'$\alpha=10^{\circ}$')
axEEPAA_LF.set_ylabel('Energy [eV]')
axEEPAA_LF.set_yscale('log')
axEEPAA_LF.set_ylim(30, 1500)



# LF EEPAA pitch data
axEEPAA_pitch_LF = fig.add_subplot(gs0[1, 2])
axEEPAA_pitch_LF.get_xaxis().set_visible(False)
xData, yData = invariantLat_ptcle_low, Pitch
zData = np.transpose(data_dict_eepaa_low['Differential_Energy_Flux'][0][:, :, wPitch_engy_val])
cmap_EEPAA_pitch_LF = axEEPAA_pitch_LF.pcolormesh(xData, yData, zData, vmin=cbarLimits[0], vmax=cbarLimits[1], shading='nearest',cmap='turbo',norm='log')
axEEPAA_pitch_LF.set_ylabel('Pitch Angle [$^{\circ}$]\n'+f' Energy ={Energy[32]}eV')
axEEPAA_pitch_LF.set_ylim(0,180)



# LF E data
axEB_LF_E = fig.add_subplot(gs0[2, 2])
axEB_LF_E.get_xaxis().set_visible(False)
axEB_LF_E.plot(invariantLat_wave_low, data_dict_elec_low['E_r'][0], color='tab:red')
axEB_LF_E.set_ylabel('$\delta$Er [mV/m]',color='tab:red')
axEB_LF_E.set_ylim(-12, 12)
axEB_LF_E.spines['left'].set_color('tab:red')
axEB_LF_E.xaxis.label.set_color('tab:red')
axEB_LF_E.tick_params(axis='y',colors='tab:red')
axEB_LF_E.margins(x=0)


# LF B data
axEB_LF_B = fig.add_subplot(gs0[3, 2])
axEB_LF_B.plot(invariantLat_wave_low, data_dict_mag_low['B_e'][0], color='tab:blue')
axEB_LF_B.get_xaxis().set_visible(False)
axEB_LF_B.set_ylabel('$\delta$Be [nT]',color='tab:blue')
axEB_LF_B.set_ylim(-10, 10)
axEB_LF_B.spines['left'].set_color('tab:blue')
axEB_LF_B.xaxis.label.set_color('tab:blue')
axEB_LF_B.tick_params(axis='y',colors='tab:blue')
axEB_LF_B.margins(x=0)

axEB_LF_B.legend()




#%%%LF Poynting Flux data%%%
# axS_LF = fig.add_subplot(gs0[4, 2])
# axS_LF.plot(invariantLat_wave_low, (1E3)*data_dict_poynting_low['S_p'][0], label='Sp', color='black')
# axS_LF.plot(invariantLat_wave_low, (1E3)*data_dict_poynting_low['S_e'][0],label='Se', color='tab:red',alpha=0.7)
# axS_LF.plot(invariantLat_wave_low, (1E3)*data_dict_poynting_low['S_r'][0], label='Sr', color='tab:blue',alpha=0.7)
# axS_LF.set_ylabel('Energy Flux \n [ergs/cm$^2$s]')
# axS_LF.set_ylim(-0.025,0.025)
# axS_LF.annotate('ILat [$^{\circ}$]',xy=(-0.15,-0.04),xytext=(0,0),ha='left',va='top',xycoords='axes fraction',textcoords='offset points')
# axS_LF.xaxis.set_major_locator(ticker.LinearLocator(N))
# axS_LF.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f')) # rounds the ILat axis to 2 decimal places
# axS_LF.margins(x=0)
# axS_LF.legend(loc='upper right')

#add additional x axes(ALT)
axS_LF.annotate('Alt [km]',xy=(-0.15,-0.16),xytext=(0,0),ha='left',va='top',xycoords='axes fraction',textcoords='offset points')
axS_LF_Alt = axS_LF.twiny() # create second axis
axS_LF_Alt.spines['bottom'].set_position(('axes',-0.14)) # move the axes to the bottom
axS_LF_Alt.xaxis.set_ticks_position("bottom") # move the ticks to the bottom
axS_LF_Alt.set_frame_on(True)
axS_LF_Alt.patch.set_visible(False)
axS_LF_Alt.spines['bottom'].set_visible(False)
axS_LF_Alt.tick_params(axis='x',length=0)
axS_LF_Alt.plot(data_dict_traj_low['geoAlt'][0],data_dict_traj_low['geomagAlt'][0], alpha=0)
axS_LF_Alt.xaxis.set_major_locator(ticker.LinearLocator(N))
axS_LF_Alt.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f')) # rounds the ILat axis to 2 decimal places


#add additional x axes(UT Time)
axS_LF.annotate('UT ',xy=(-0.15, -0.28),xytext=(0,0),ha='left',va='top',xycoords='axes fraction',textcoords='offset points')
axS_LF_UT = axS_LF.twiny() # create second axis
axS_LF_UT.spines['bottom'].set_position(('axes',-0.26)) # move the axes to the bottom
axS_LF_UT.xaxis.set_ticks_position("bottom") # move the ticks to the bottom
axS_LF_UT.set_frame_on(True)
axS_LF_UT.patch.set_visible(False)
axS_LF_UT.spines['bottom'].set_visible(False)
axS_LF_UT.tick_params(axis='x',length=0)
axS_LF_UT.plot(data_dict_traj_low['Epoch_esa'][0], data_dict_traj_low['geomagAlt'][0], alpha=0)
axS_LF_UT.xaxis.set_major_locator(ticker.LinearLocator(N))
axS_LF_UT.set_xticklabels(UTTicks_low)

#add additional x axes(ALT)

# --- 1st Colorbar Column- LF ---
axCbar_LF = fig.add_subplot(gs0[0, 3])
cbar_LF = fig.colorbar(mappable=cmap_EEPAA_LF,
                       cax=axCbar_LF,
                       norm='log')
cbar_LF.set_label('eV/cm$^{2}$-s-sr-eV', labelpad=-52)
l, b, w, h = axCbar_LF.get_position().bounds
axCbar_LF.set_position([(0.975)*l, b, w, h])

# --- 2nd Colorbar Column - LF ---
axCbar_pitch_LF = fig.add_subplot(gs0[1, 3])
cbar_pitch_LF = fig.colorbar(mappable=cmap_EEPAA_pitch_LF,
                       cax=axCbar_pitch_LF,
                       norm='log')
cbar_pitch_LF.set_label('eV/cm$^{2}$-s-sr-eV', labelpad=-52)
l, b, w, h = axCbar_pitch_LF.get_position().bounds
axCbar_pitch_LF.set_position([(0.975)*l, b, w, h])




# --- SHOW PLOT ---
# plt.tight_layout()
plt.savefig(r'C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\Papers\ACESII_Alfvenic_Observations\Plots\\Plot2_Conjugacy.png')
# plt.show()