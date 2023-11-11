# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: All Sky Imager data plot



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
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


# --- INPUT FILES ---
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
    data_dict_traj = loadDictFromFile(inputFilesTraj[i], input_data_dict={},reduceData=False,targetTimes=[],wKeys=[])
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
        strTimeStamp = imageStr.replace(f'{pathToSite}\\','').replace(f'{WLengths[i]}\\','').replace('skn4_','').replace(f'_{WLengths[i]}_cal.png','')
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


# figure info
# plt.style.use('dark_background')
fig = plt.figure(dpi=100)
figure_height = 20
figure_width = 36

fig.set_figwidth(figure_width)
fig.set_figheight(figure_height)

# title
fig.suptitle('ACESII')

##### Define the primary Gridspec #####
gs0 = gridspec.GridSpec(nrows=2, ncols=1, figure=fig, height_ratios=[0.97, 0.03]) # splits figure between the plots and the colorbar at the very bottom
cbarVmin = 0
cbarVmax = 18


# --- Split Data-half into two columns ---
gsData = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, width_ratios=[2/3, 1/3], subplot_spec=gs0[0])


# -*- Altitude/AllSky Plots -*-
gs_altLat_BigAllSky = gridspec.GridSpecFromSubplotSpec(nrows=2,ncols=1,height_ratios=[1/3,2/3],subplot_spec=gsData[0])

# AltLat plot
axAltLat = fig.add_subplot(gs_altLat_BigAllSky[0])

axAltLat.plot(geoMagLat[0], geoAlt[0], color='tab:red') # High
axAltLat.plot(geoMagLat[1], geoAlt[1],color='tab:blue') # Low
axAltLat.set_ylabel('Altitude [km]')
axAltLat.set_xlabel('Geomagnetic Lat [deg]')
axAltLat.set_xlim(68, 72)

axGeographicLat = axAltLat.twiny()
axGeographicLat.plot(geoLat[0], geoAlt[0], color='tab:red') # High
axGeographicLat.plot(geoLat[1], geoAlt[1], color='tab:blue') # Low
axGeographicLat.set_xlabel('Geographic Lat [deg]')




# BigAllSky plot
projPC = ccrs.PlateCarree()  # MUST be kept on the set_extent, crs =, command AND pcolormesh transform command
lonW = 12
lonE = 18
latS = 68
latN = 74.5
cLat = (latN + latS) / 2
cLon = (lonW + lonE) / 2
res = '10m'

# projType = ccrs.Mollweide(central_longitude=cLon, central_latitude=cLat)
projType = ccrs.Stereographic(central_longitude=cLon, central_latitude=cLat)

axBigAllSky = fig.add_subplot(gs_altLat_BigAllSky[1],projection=projType)



# allSkyCmap = 'inferno'
allSkyCmap = 'turbo'
LowFlyerColor = 'royalblue'
HighFlyerColor = 'red'
LowFlyerProjectionColor = 'cyan'
HighFlyerProjectionColor = 'darkorange'

AltvsLatLineWidth = 4
LatvsLongLineWidth = 4
AllSkyLineWidth = 4

axBigAllSky.gridlines(draw_labels=True, linewidth=3, alpha=0.4, linestyle='--')
axBigAllSky.xlabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
axBigAllSky.ylabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
axBigAllSky.set_extent([lonW, lonE, latS, latN], crs=projPC)  # controls lat/long axes display
axBigAllSky.coastlines(resolution=res, color='white', alpha=0.8)  # adds coastlines with resolution
axBigAllSky.set_aspect(0.3)

cmapBigAllSky = axBigAllSky.pcolormesh(allGLongs[0], allGlats[0], allImages[0][0], cmap=allSkyCmap, transform=projPC,vmin=cbarVmin,vmax=cbarVmax)



# -*- 5 AllSky Plots -*-
gs_fiveAllSky = gridspec.GridSpecFromSubplotSpec(nrows=5,ncols=1,height_ratios=[1/5,1/5,1/5,1/5,1/5],subplot_spec=gsData[1])

for i in range(5):
    axSmallAllSky = fig.add_subplot(gs_fiveAllSky[i],projection=projType)
    axSmallAllSky.gridlines(draw_labels=True, linewidth=3, alpha=0.4, linestyle='--')
    axSmallAllSky.xlabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
    axSmallAllSky.ylabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
    axSmallAllSky.set_extent([lonW, lonE, latS, latN], crs=projPC)  # controls lat/long axes display
    axSmallAllSky.coastlines(resolution=res, color='white', alpha=0.8)  # adds coastlines with resolution
    axSmallAllSky.set_aspect(0.3)

    axSmallAllSky.pcolormesh(allGLongs[0], allGlats[0], allImages[0][0], cmap=allSkyCmap, transform=projPC, vmin=cbarVmin,vmax=cbarVmax)




# --- Define colorbar axis ---
cmap = 'turbo'
gsColorBar = fig.add_subplot(gs0[1])
plt.colorbar(mappable=cmapBigAllSky, cax=gsColorBar, orientation='horizontal',fraction=0.046, pad=0.04)
gsColorBar.set_title('Intensity [kR]', fontsize=20)
gsColorBar.tick_params(labelsize=20)

plt.savefig(r'C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\Papers\ACESII_Alfvenic_Observations\Plots\\Plot1_AllSky.png')
# plt.show()
