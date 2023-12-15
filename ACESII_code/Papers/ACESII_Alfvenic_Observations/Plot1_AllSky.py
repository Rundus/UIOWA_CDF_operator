# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: All Sky Imager data plot



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
from scipy.io import readsav
from math import isnan


print(color.UNDERLINE + f'Plot1_AllSky' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
architecture = 65535 # type of bit-architecture used for the allSky images. Nomially 16-bit
elevlimits = [20, 20] # value (in deg) of the cutoff elevation angle for the 5577A and 6300A line

# --- --- --- --- -
# --- LOAD DATA ---
# --- --- --- --- -

# trajectory
trajFolderPath = f'{ACES_data_folder}trajectories\\'
inputFilesTraj = [glob(trajFolderPath + rf'{fliers[0]}\\*_ILat_ILong*')[0], glob(trajFolderPath + rf'{fliers[1]}\\\\*_ILat_ILong*')[0]]

# ALl sky calibration
pathToSite = r'C:\Data\ACESII\all_sky\skibotn'
calFiles = [readsav(glob(pathToSite+'\\5577\\*.dat*')[0]), readsav(glob(pathToSite + '\\6300\\*.dat*')[0])]

# allsky image files
photoFiles = [glob(pathToSite + '\\5577\\*.png*'), glob(pathToSite + '\\6300\\*.png*')]


# --- GET TRAJECTORY DATA ---
prgMsg(f'Loading ACESII traj data')
data_dicts_traj = []
for i in range(2):
    data_dict_traj = loadDictFromFile(inputFilesTraj[i], input_data_dict={}, reduceData=False,targetTimes=[], wKeys=[])
    data_dict_traj['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch_esa'][0][i]) for i in (range(len(data_dict_traj['Epoch_esa'][0])))])
    data_dicts_traj.append(data_dict_traj)

# define some variables
EpochRocket = [data_dicts_traj[0]['Epoch_esa'][0], data_dicts_traj[1]['Epoch_esa'][0]]
geoAlt = [data_dicts_traj[0]['geoAlt'][0], data_dicts_traj[1]['geoAlt'][0]]
geoLat = [data_dicts_traj[0]['geoLat'][0], data_dicts_traj[1]['geoLat'][0]]
geoLong = [data_dicts_traj[0]['geoLong'][0], data_dicts_traj[1]['geoLong'][0]]
geoMagLat = [data_dicts_traj[0]['geomagLat'][0], data_dicts_traj[1]['geomagLat'][0]]
Done(start_time)


# --- GET ALL SKY CALIBRATION DATA ---
prgMsg('Collecting and processing Image data')
allGlats = [np.array(deepcopy(calFiles[0]['glats']))[::-1], np.array(deepcopy(calFiles[1]['glats']))[::-1]]
allGLongs = [np.array(deepcopy(calFiles[0]['glons']))[::-1], np.array(deepcopy(calFiles[1]['glons']))[::-1]]
allElevs = [np.array(deepcopy(calFiles[0]['elevs'])), np.array(deepcopy(calFiles[1]['elevs']))]


# --- GET ALL SKY IMAGE DATA ---
# get the image time series and the data itself into single variables
Epoch_AllSky = [[], []]
imageData = [[], []]
WLengths = ['5577', '6300']
prgMsg('Collecting Image Data')
for i in range(len(photoFiles)):
    for imageStr in photoFiles[i]:

        # get the timestamp
        strTimeStamp = imageStr.replace(f'{pathToSite}\\', '').replace(f'{WLengths[i]}\\', '').replace('skn4_', '').replace(f'_{WLengths[i]}_cal.png', '')
        year = int(strTimeStamp[0:4])
        month = int(strTimeStamp[4:6])
        day = int(strTimeStamp[6:8])
        hour = int(strTimeStamp[9:11])
        minute = int(strTimeStamp[11:13])
        second = int(strTimeStamp[13:15])
        dtTimeStamp = dt.datetime(year, month, day, hour, minute, second)
        Epoch_AllSky[i].append(dtTimeStamp)

        # get the grayscale data
        imageData[i].append(plt.imread(imageStr))

Done(start_time)


# -- collect and process all aurora image data ---
prgMsg('Converting All Sky Data to Rayleighs')
# remove nan values from data and replace with garbage AND convert all Images into Rayleighs
# description: the image data is a 16-bit counts value normalized by 65535.
# Invert this and use the calibration factor of 1 R/count given in the cal file
allImages = imageData

for i in range(2):  # wavelength

    # --- correct the calibration data ---
    for j in range(len(allGlats[i])):  # array rows
        for k in range(len(allGlats[i][j])):  # row values
            if isnan(allGlats[i][j][k]):
                allGlats[i][j][k] = 70
                for a in range(len(allImages[i])):  # correct this j,k point in all auroral images
                    allImages[i][a][j][k] = np.nan

            if isnan(allGLongs[i][j][k]):
                allGLongs[i][j][k] = 20
                for a in range(len(allImages[i])):  # correct this j,k point in all auroral images
                    allImages[i][a][j][k] = np.nan

            if allElevs[i][j][k] <= elevlimits[i] or isnan(allElevs[i][j][k]):
                for a in range(len(allImages[i])):  # correct this j,k point in all auroral images
                    allImages[i][a][j][k] = np.nan

    # --- convert images to rayleighs ---
    for a in range(len(allImages[i])):  # number of images
        for j in range(len(allImages[i][a])):  # arrays in particular image
            for k in range(len(allImages[i][a][j])):  # for each value in each array of image
                if not isnan(allImages[i][a][j][k]):
                    allImages[i][a][j][k] = int(allImages[i][a][j][k] * architecture)/1000

Done(start_time)



############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

# some pre-amble
trajColors = ['red','black']
cmapColor = 'viridis'
faceColorChoice = (156 / 255, 156 / 255, 156 / 255, 0.5) # in normalize RGBA

# --- figure info  ---
# plt.style.use('dark_background')
fig = plt.figure()
figure_height = 8
figure_width = 8

fig.set_figwidth(figure_width)
fig.set_figheight(figure_height)
fig.suptitle('Skibotn AllSky (5570 A)',fontsize=10,fontweight='bold')

# --- title ---
# fig.suptitle('ACES II')


##### Define the primary Gridspec #####
gs0 = gridspec.GridSpec(nrows=2, ncols=1, figure=fig, height_ratios=[0.99, 0.01], hspace=0.25) # splits figure between the plots and the colorbar at the very bottom
cbarVmin = 0
cbarVmax = 14



# --- Split Data-half into two columns ---
gsData = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, width_ratios=[1/5, 1/5], subplot_spec=gs0[0])

# -*- Altitude/AllSky Plots -*-
gs_altLat_BigAllSky = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, height_ratios=[3/10, 7/10], subplot_spec=gsData[0], hspace=0.3)

# --- AltLat plot ---
axAltLat = fig.add_subplot(gs_altLat_BigAllSky[0])
axAltLat.set_ylabel('Altitude [km]')
axAltLat.set_xlabel('Geomagnetic Lat [deg]')
axAltLat.set_xlim(67.25, 72.5)
axAltLat.set_ylim(0, 460)


# plot the psudo geomagnetic field line
slope = -1*(111/np.sin(np.radians(90 - 78.13))) # corresponds to line with -78.13deg inclination
for i in range(21):
    axAltLat.axline(xy1=(66+i*0.5, 0), slope=slope,color='tab:blue', linewidth=2, linestyle='-.', alpha=0.3)

# set the facecolor of the axAltLat plot
axAltLat.set_facecolor(faceColorChoice)


# plot the UTC labels
axGeographicLat = axAltLat.twiny()
axGeographicLat.plot(geoLat[0], geoAlt[0], color=trajColors[0], alpha=0) # High
axGeographicLat.plot(geoLat[1], geoAlt[1], color=trajColors[1], alpha=0) # Low
axGeographicLat.set_xlabel('Geographic Lat [deg]')
timeTargetsUTC_labels = [200, 250, 300, 350, 400, 450]
AltLat_vertical_Alignments = ['bottom' for tme in timeTargetsUTC_labels]
AltLat_horizontal_Alignments = ['right', 'right', 'right', 'left', 'left', 'left']
vertical_text_label_adjustments = [-0.05, -0.035, 0.01, 0.015, -0.04, -0.04]
horizontal_text_label_adjustments = [-0.002, -0.001, 0.000, 0.0005, 0.002, 0.001]

# plot the scatterpoint of each of the timeTargetUTC_labels. Plot the text itself only for the High Flyer and
# create a connecting line between the scatterpoints between the flyers
for i in range(2): # for each rocket
    for j, ttme in enumerate(timeTargetsUTC_labels):
        Index = np.abs(EpochRocket[i] - (ttme*1E9 + EpochRocket[i][0])).argmin()
        xPos = geoMagLat[i][Index]
        yPos = geoAlt[i][Index]

        if i == 0:
            # plot the text itself
            label = pycdf.lib.tt2000_to_datetime(EpochRocket[i][Index])
            deltaY = vertical_text_label_adjustments[j]*yPos
            deltaX = horizontal_text_label_adjustments[j]*xPos
            axAltLat.text(x=xPos+deltaX, y=yPos+deltaY, s=label.strftime("%H:%M:%S"), color=trajColors[i], va=AltLat_vertical_Alignments[j], ha=AltLat_horizontal_Alignments[j], size=10)

            # plot the connecting line
            Index_LF = Index = np.abs(EpochRocket[1] - (ttme*1E9 + EpochRocket[1][0])).argmin()
            xPos_LF = geoMagLat[1][Index_LF]
            yPos_LF = geoAlt[1][Index_LF]
            axAltLat.plot([xPos,xPos_LF],[yPos,yPos_LF], color='green', linestyle='--', alpha=0.5)
            # axAltLat.axline(xy1=(xPos, yPos), xy2=(xPos_LF, yPos_LF), color='green', linestyle='--', alpha=0.5)

        # plot a dot at the text label
        axAltLat.scatter(x=xPos, y=yPos, s=15, marker="o", color=trajColors[i])

# plot the trajectory over everything
axAltLat.plot(geoMagLat[0], geoAlt[0], color=trajColors[0], label='High Flyer') # High
axAltLat.plot(geoMagLat[1], geoAlt[1], color=trajColors[1], label='Low Flyer') # Low



# --- --- --- --- --- --
# --- BigAllSky plot ---
# --- --- --- --- --- --
projPC = ccrs.PlateCarree()  # MUST be kept on the set_extent, crs =, command AND pcolormesh transform command
lonW = 10
lonE = 23.5
latS = 68
latN = 74.5
cLat = (latN + latS) / 2
cLon = (lonW + lonE) / 2
res = '10m'

projType = ccrs.Stereographic(central_longitude=cLon, central_latitude=cLat)
axBigAllSky = fig.add_subplot(gs_altLat_BigAllSky[1], projection=projType)
alignment = ['left', 'right']
gl = axBigAllSky.gridlines(draw_labels=True, linewidth=1.5, alpha=0.6, linestyle='--')
gl.xlabel_style = {'size': 10, 'color': 'black', 'weight': 'bold'}
gl.ylabel_style = {'size': 10, 'color': 'black', 'weight': 'bold'}
gl.top_labels = False
axBigAllSky.set_extent([lonW, lonE, latS, latN], crs=projPC)  # controls lat/long axes display
axBigAllSky.coastlines(resolution=res, color='black', alpha=0.8)  # adds coastlines with resolution
# axBigAllSky.set_aspect(1)
cmapBigAllSky = axBigAllSky.pcolormesh(allGLongs[0], allGlats[0], allImages[0][0], cmap=cmapColor, transform=projPC, vmin=cbarVmin,vmax=cbarVmax)
axBigAllSky.set_facecolor(faceColorChoice)


# --- plot the rocket trajectory data on the large AllSky plot ---
axBigAllSky.plot(geoLong[0], geoLat[0], color=trajColors[0], transform=projPC) # High
axBigAllSky.plot(geoLong[1], geoLat[1], color=trajColors[1], transform=projPC) # Low

# plot specific UTC times on the trajectory lines

timeTargetsUTC_labels = [dt.datetime(2022, 11, 20, 17, 20, 10, 000),
                        dt.datetime(2022, 11, 20, 17, 21, 40, 000),
                        dt.datetime(2022, 11, 20, 17, 23, 10, 000),
                        dt.datetime(2022, 11, 20, 17, 24, 40, 000),
                        dt.datetime(2022, 11, 20, 17, 26, 10, 000),
                        dt.datetime(2022, 11, 20, 17, 27, 40, 000)
                        ]


for i in range(2): # for each rocket
    for ttme in timeTargetsUTC_labels:
        Index = np.abs(EpochRocket[i] - pycdf.lib.datetime_to_tt2000(ttme)).argmin()
        label = pycdf.lib.tt2000_to_datetime(EpochRocket[i][Index])

        xPos = geoLong[i][Index]
        yPos = geoLat[i][Index]

        # plot a dot at the text label
        axBigAllSky.scatter(x=xPos,y=yPos,s=15,marker="o",color=trajColors[i],transform=projPC)

        # plot the text itself
        deltaX = 0.03*xPos if i == 0 else -1*0.03*xPos
        axBigAllSky.text(x=xPos + deltaX, y=yPos, s =label.strftime("%H:%M:%S"),color=trajColors[i], ha=alignment[i], transform=projPC,size=12)





# --- --- --- --- --- --
# --- 5 AllSky Plots ---
# --- --- --- --- --- --
smallPlotTargetTimes = [dt.datetime(2022, 11, 20, 17, 20, 00, 000),
                        dt.datetime(2022, 11, 20, 17, 21, 30, 000),
                        dt.datetime(2022, 11, 20, 17, 23, 00, 000),
                        dt.datetime(2022, 11, 20, 17, 24, 30, 000),
                        dt.datetime(2022, 11, 20, 17, 26, 00, 000),
                        dt.datetime(2022, 11, 20, 17, 27, 30, 000)
                        ] # the dt corresponding to the AllSky image I want to use in the plot

nRows= int(len(smallPlotTargetTimes)/2)
nCols = 2
gs_smallerAllSkys = gridspec.GridSpecFromSubplotSpec(nrows=nRows, ncols=nCols, subplot_spec=gsData[1], hspace=0.15)

counter = 0
for i in range(nRows):
    for j in range(nCols):

        axSmallAllSky = fig.add_subplot(gs_smallerAllSkys[i, j], projection=projType)
        gl = axSmallAllSky.gridlines(draw_labels=True, linewidth=1, alpha=0.35, linestyle='--',color='k')
        axSmallAllSky.set_facecolor(faceColorChoice)
        gl.xlabel_style = {'size': 0, 'color': 'black', 'weight': 'bold'}
        gl.ylabel_style = {'size': 0, 'color': 'black', 'weight': 'bold'}
        gl.top_labels = False
        axSmallAllSky.set_extent([lonW, lonE, latS, latN], crs=projPC)  # controls lat/long axes display
        axSmallAllSky.coastlines(resolution=res, color='black', alpha=0.8)  # adds coastlines with resolution
        # axSmallAllSky.set_aspect(0.3)
        pltIndex = np.abs(np.array(Epoch_AllSky[0]) - smallPlotTargetTimes[counter]).argmin()
        axSmallAllSky.set_title(Epoch_AllSky[0][pltIndex].strftime("%H:%M:%S") + ' UTC', fontsize=10, weight='bold')
        axSmallAllSky.pcolormesh(allGLongs[0], allGlats[0], allImages[0][pltIndex], cmap=cmapColor, transform=projPC, vmin=cbarVmin, vmax=cbarVmax)


        # plot the rocket trajectories on all the small AllSky Plots
        # axSmallAllSky.plot(geoLong[0], geoLat[0], color=trajColors[0], transform=projPC)  # High
        # axSmallAllSky.plot(geoLong[1], geoLat[1], color=trajColors[1], transform=projPC)  # Low
        # for k in range(2):
        #     Index = np.abs(EpochRocket[k] - pycdf.lib.datetime_to_tt2000(smallPlotTargetTimes[counter])).argmin()
        #     xPos = geoLong[k][Index]
        #     yPos = geoLat[k][Index]
        #     axSmallAllSky.scatter(x=xPos, y=yPos, marker="o", color=trajColors[k], transform=projPC)

        counter += 1


# --- --- --- ----
# --- COLORBAR ---
# --- --- --- ----
gsColorBar = fig.add_subplot(gs0[1])
cbar = plt.colorbar(mappable=cmapBigAllSky, cax=gsColorBar, orientation='horizontal',fraction=0.046, pad=0.04)
cbar.set_label('Intensity [kR]', fontsize=15)
gsColorBar.tick_params(labelsize=20)

plt.tight_layout()
plt.savefig(r'C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\Papers\ACESII_Alfvenic_Observations\Plots\\Plot1_AllSky.png')
# plt.show()
