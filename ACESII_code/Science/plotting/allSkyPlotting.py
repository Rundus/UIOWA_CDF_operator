# --- allSkyPlotting.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Loads in the AllSky data, uses the calibration file to determine position
# finally loads in traj data to determine rocket trajectory


# assumes all light comes from these altitudes:
# 557.7nm: 150km
# 630nm: 250km


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools
import math
# --- --- --- --- ---

import time
from ACESII_code.class_var_func import Done, setupPYCDF

start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintSiteNames = False

# --- Select the Site ---
# 0 -> Skibotn
wSite = 0

modifier = ''
inputPath_modifier_AllSky = 'all_sky' # e.g. 'L1' or 'L1'. It's the name of the broader input folder inside data\ACESII
inputPath_modifier_traj = 'trajectories'
inputPath_modifier_magGeo = 'mag'
outputPath_modifier = '\\trajectories\\trajectory_plots' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder inside data\ACESII\ACESII_matlab

architecture = 65535 # type of bit-architecture used for the allSky images. Nomially 16-bit
elevlimits = [20, 20] # value (in deg) of the cutoff elevation angle for the 5577A and 6300A line

# --- AllSky MOVIE ---
createAllSkyMovie = True
fps = 1 # fps of the video
projectionAltitude = 100 # in km
useRealMagData = False # uses inflight mag data to determine projections for Alt vs Lat plot. Else uses IGRF

# plot specific locations
plotSpecificLocations = False
specificLocations = [i*10 for i in range(int(12371/10))]

# plot frame skips
doFrameSkips = True
frame_skips = 750 # how many Epoch_esa frames to skip

# Kenton wanted just the low flyer for a plot, so here's a toggle to do that
useOnlyLowFlyer = False

# scott didn't like the latlong plot, so here's something to remove it
removeLatLongPlot = True



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import datetime as dt
from matplotlib import pyplot as plt, animation
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
from copy import deepcopy
from ACESII_code.missionAttributes import ACES_mission_dicts
from ACESII_code.data_paths import ACES_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg, loadDictFromFile
from glob import glob
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)
from scipy.io import readsav
import pyIGRF


def allSkyIDL_to_py(wSite, justPrintSiteNames, rocketFolderPath):

    # --- load attributes for ACESII traj data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID

    allSkySiteFolder = f'{rocketFolderPath}\\{inputPath_modifier_AllSky}'
    sitePaths = glob(f'{allSkySiteFolder}\\*')
    sites = [path.replace(f'{allSkySiteFolder}\\','') for path in sitePaths]
    wSiteName = sites[wSite]
    pathToSite = sitePaths[wSite]
    WLengths = ['5577', '6300']

    if justPrintSiteNames:
        for i, file in enumerate(sites):
            print('[{:.0f}] {:80s}'.format(i, sites[i]))
    else:

        #############################
        # --- get the input files ---
        #############################

        # --- get cal files and convert to python---
        calFiles = [readsav(glob(pathToSite+'\\5577\\*.dat*')[0]),readsav(glob(pathToSite + '\\6300\\*.dat*')[0])]

        # --- get All Sky photos ---
        photoFiles = [glob(pathToSite + '\\5577\\*.png*'),glob(pathToSite + '\\6300\\*.png*')]

        # --- traj Data ---
        trajFolderPath = f'{ACES_data_folder}trajectories\\'
        inputFilesTraj = [glob(trajFolderPath + rf'{fliers[0]}\\*_ILat_ILong*')[0],
                         glob(trajFolderPath + rf'{fliers[1]}\\\\*_ILat_ILong*')[0]]

        # --- magGeo Data ---
        magGeoFolderPath = f'{ACES_data_folder}mag\\'
        inputFilesMagGeo = [glob(f'{magGeoFolderPath}\\{fliers[0]}\\*RingCore_Geo_despun*')[0],
                            glob(f'{magGeoFolderPath}\\{fliers[1]}\\*RingCore_Geo_despun*')[0]]

        ###############
        # --- START ---
        ###############
        print('\n')
        print(color.UNDERLINE + f'Plotting allSky Data for {wSiteName} data' + color.END)

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading ACESII traj data')
        data_dicts_traj = []
        for i in range(2):
            data_dict_traj = loadDictFromFile(inputFilesTraj[i], {})
            data_dict_traj['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch_esa'][0][i]) for i in (range(len(data_dict_traj['Epoch_esa'][0])))])
            data_dicts_traj.append(data_dict_traj)
        Done(start_time)

        # --- get data from the GeoMag files ---
        prgMsg(f'Loading ACESII RingCore data')
        data_dicts_geoMag = []
        for i in range(2):
            data_dict_geoMag = loadDictFromFile(inputFilesMagGeo[i],{})
            data_dict_geoMag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_geoMag['Epoch'][0][i]) for i in range(len(data_dict_geoMag['Epoch'][0]))])
            data_dicts_geoMag.append(data_dict_geoMag)


        # --- COLLECT IMAGE FILES AND TIMESTAMPS ---

        # get the image time series and the data itself into single variables
        Epoch_AllSky = [[], []]
        imageData = [[], []]
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

        ###################################
        # --- prepare data for plotting ---
        ###################################
        allGlats = [np.array(deepcopy(calFiles[0]['glats']))[::-1], np.array(deepcopy(calFiles[1]['glats']))[::-1]]
        allGLongs = [np.array(deepcopy(calFiles[0]['glons']))[::-1], np.array(deepcopy(calFiles[1]['glons']))[::-1]]
        allElevs = [np.array(deepcopy(calFiles[0]['elevs'])), np.array(deepcopy(calFiles[1]['elevs']))]
        allImages = imageData
        EpochRocket = [data_dicts_traj[0]['Epoch_esa'][0], data_dicts_traj[1]['Epoch_esa'][0]]
        geoAlt = [data_dicts_traj[0]['geoAlt'][0], data_dicts_traj[1]['geoAlt'][0]]
        geoLat = [data_dicts_traj[0]['geoLat'][0], data_dicts_traj[1]['geoLat'][0]]
        geoLong = [data_dicts_traj[0]['geoLong'][0], data_dicts_traj[1]['geoLong'][0]]
        geoILat = [data_dicts_traj[0]['geoILat'][0], data_dicts_traj[1]['geoILat'][0]]
        geoILong = [data_dicts_traj[0]['geoILong'][0], data_dicts_traj[1]['geoILong'][0]]



        # downsample the RingCore mag data to match the rocket data
        prgMsg('Downsampling B-Field')
        # create an array of indicies, format [[],[]]
        alignedMagIndicies = np.array([[np.abs(data_dicts_geoMag[i]['Epoch'][0] - EpochRocket[i][j]).argmin() for j in range(len(EpochRocket[i]))] for i in range(2)], dtype='object')

        # use those indices to downsample B_enu into format [[  [E,N,U],[E,N,U]....  ],[...]]
        B_ENU = np.array([[[data_dicts_geoMag[i]['B_east'][0][j], data_dicts_geoMag[i]['B_north'][0][j], data_dicts_geoMag[i]['B_up'][0][j]] for j in range(len(alignedMagIndicies[i]))] for i in range(2)], dtype='object')
        Done(start_time)

        prgMsg('Collecting and processing/cleaning Image data')

        # -- collect and process all aurora image data ---
        # remove nan values from data and replace with garbage AND convert all Images into Rayleighs
        # description: the image data is a 16-bit counts value normalized by 65535. Invert this and use the calibration factor of 1 R/count given in the cal file

        for i in range(2):  # wavelength

            # --- correct the calibration data ---
            for j in range(len(allGlats[i])):  # array rows
                for k in range(len(allGlats[i][j])):  # row values
                    if math.isnan(allGlats[i][j][k]):
                        allGlats[i][j][k] = 70
                        for a in range(len(allImages[i])):  # correct this j,k point in all auroral images
                            allImages[i][a][j][k] = np.nan

                    if math.isnan(allGLongs[i][j][k]):
                        allGLongs[i][j][k] = 20
                        for a in range(len(allImages[i])):  # correct this j,k point in all auroral images
                            allImages[i][a][j][k] = np.nan

                    if allElevs[i][j][k] <= elevlimits[i] or math.isnan(allElevs[i][j][k]):
                        for a in range(len(allImages[i])):  # correct this j,k point in all auroral images
                            allImages[i][a][j][k] = np.nan

            # --- convert images to rayleighs ---
            for a in range(len(allImages[i])): #number of images
                for j in range(len(allImages[i][a])): # arrays in particular image
                    for k in range(len(allImages[i][a][j])): # for each value in each array of image
                        if not math.isnan(allImages[i][a][j][k]):
                            allImages[i][a][j][k] = int(allImages[i][a][j][k]*architecture)

        Done(start_time)

        # --- determine the timestamps for each photo and rocket data ---
        Epoch_AllSky_tt2000 = [[pycdf.lib.datetime_to_tt2000(time) for time in Epoch_AllSky[0]],
                               [pycdf.lib.datetime_to_tt2000(time) for time in Epoch_AllSky[1]]]

        # for each image timestamp, find the index of the closest Epoch_esa in the High flyer's epoch
        imageIndiciesToRocket = [[], []]
        for i in range(len(Epoch_AllSky_tt2000)):
            for stamp in Epoch_AllSky_tt2000[i]:
                imageIndiciesToRocket[i].append(np.abs(data_dicts_traj[0]['Epoch_esa'][0] - stamp).argmin())

        imageIndiciesToRocket = np.array(imageIndiciesToRocket)

        # fill in imageIndiciesToRocket with indicies to match the size of High Flyer Epoch data
        imgIndicies = [[], []]

        for i in range(2):
            for j in range(len(imageIndiciesToRocket[i])):  # loop through the various images
                if j == 0:

                    for k in range(imageIndiciesToRocket[i][j + 1]):
                        imgIndicies[i].append(j)

                elif j == len(imageIndiciesToRocket[i]) - 1:  # if you're at the last image index
                    for k in range(len(EpochRocket[0]) - imageIndiciesToRocket[i][j]):
                        imgIndicies[i].append(j)

                else:
                    for k in range((imageIndiciesToRocket[i][j + 1] - imageIndiciesToRocket[i][j])):
                        imgIndicies[i].append(j)

        imgIndicies = np.array(imgIndicies)
        EpochMovie = np.array([pycdf.lib.tt2000_to_datetime(EpochRocket[0][i]) for i in range(len(EpochRocket[0]))])

        ###############################
        # --- EXTEND LOW FLYER DATA ---
        ###############################

        # --- extend Low Flyer Rocket data to be the same length as High flyer in the beginning and end ---
        highFlyerSize, lowFlyerSize = len(geoAlt[0]), len(geoAlt[1])

        # --- Append start Values to  ---
        no_of_points_start = int((EpochRocket[1][0] - EpochRocket[0][0]) / (rocketAttrs.MinorFrameTime))
        newAlt = [geoAlt[1][0] for i in range(no_of_points_start)]
        newLat = [geoLat[1][0] for i in range(no_of_points_start)]
        newILat = [geoILat[1][0] for i in range(no_of_points_start)]
        newLong = [geoLong[1][0] for i in range(no_of_points_start)]
        newILong = [geoILong[1][0] for i in range(no_of_points_start)]
        newB_ENU = [B_ENU[1][0] for i in range(no_of_points_start)]

        for i in range(len(geoAlt[1])):
            newAlt.append(geoAlt[1][i])
            newLat.append(geoLat[1][i])
            newLong.append(geoLong[1][i])
            newB_ENU.append(B_ENU[1][i])
            newILat.append(geoILat[1][i])
            newILong.append(geoILong[1][i])

        # --- Append the ending values ---
        remainingIndicies = highFlyerSize - (lowFlyerSize + no_of_points_start)

        for i in range(remainingIndicies):
            newAlt.append(geoAlt[1][-1])
            newLat.append(geoLat[1][-1])
            newLong.append(geoLong[1][-1])
            newB_ENU.append(B_ENU[1][-1])
            newILat.append(geoILat[1][-1])
            newILong.append(geoILong[1][-1])


        geoAlt[1], geoLat[1], geoLong[1], B_ENU[1],geoILat[1],geoILong[1] = np.array(newAlt), np.array(newLat), np.array(newLong), np.array(newB_ENU),np.array(newILat),np.array(newILong)

        # calculate the B-projection for the Alt vs Lat plot.
        # output is a set of x,y points for each flyer. Format: [[x, geolat] [100, geoAlt]]
        projectB = [[], []]

        lat_to_meter = 111.319488  # 1 deg latitude to kilometers on Earth

        if useRealMagData:
            for i in range(2):
                for j in range(len(B_ENU[i])):

                    projectB[i].append(
                        [
                            [geoLat[i][j] + (geoAlt[i][j]-projectionAltitude)*np.abs(B_ENU[i][j][1]/B_ENU[i][j][2])/lat_to_meter, geoLat[i][j]], [projectionAltitude, geoAlt[i][j]] # determine the project latitude by triangulization: x = geoAlt * |B_North/B_up|
                        ]
                    )
        else:

            for i in range(2):
                for j in range(len(EpochRocket[0])):
                    B = pyIGRF.igrf_value(geoLat[i][j], geoLong[i][j], geoAlt[i][j], 2022)
                    B_east = B[3]
                    B_north = B[4]
                    B_up = B[5]
                    projectB[i].append(
                        [
                            [geoLat[i][j] + (geoAlt[i][j] - projectionAltitude) * np.abs(B_north / B_up) / lat_to_meter, geoLat[i][j]], [projectionAltitude, geoAlt[i][j]]
                            # determine the project latitude by triangulization: x+x0 = geoAlt * |B_North/B_up|/lat_to_meters + geoLat
                        ]
                    )


        ###########################
        # --- PLOT ALLSKY MOVIE ---
        ###########################
        if createAllSkyMovie:

            ##########################################
            # --- ASSIGN IMAGE TIMESTAMP TO ROCKET ---
            ##########################################

            # --- determine the timestamps for each photo and rocket data ---
            Epoch_AllSky_tt2000 = [[pycdf.lib.datetime_to_tt2000(time) for time in Epoch_AllSky[0]],
                                   [pycdf.lib.datetime_to_tt2000(time) for time in Epoch_AllSky[1]]]

            # for each image timestamp, find the index of the closest Epoch_esa in the High flyer's epoch
            imageIndiciesToRocket = [[], []]
            for i in range(len(Epoch_AllSky_tt2000)):
                for stamp in Epoch_AllSky_tt2000[i]:
                    imageIndiciesToRocket[i].append(np.abs(data_dicts_traj[0]['Epoch_esa'][0] - stamp).argmin())

            imageIndiciesToRocket = np.array(imageIndiciesToRocket)

            # fill in imageIndiciesToRocket with indicies to match the size of High Flyer Epoch data
            imgIndicies = [[], []]

            for i in range(2):
                for j in range(len(imageIndiciesToRocket[i])): # loop through the various images
                    if j == 0:

                        for k in range(imageIndiciesToRocket[i][j+1]):
                            imgIndicies[i].append(j)

                    elif j == len(imageIndiciesToRocket[i])-1: # if you're at the last image index
                        for k in range( len(EpochRocket[0]) - imageIndiciesToRocket[i][j] ):
                            imgIndicies[i].append(j)

                    else:
                        for k in range((imageIndiciesToRocket[i][j+1] - imageIndiciesToRocket[i][j]) ):
                            imgIndicies[i].append(j)

            imgIndicies = np.array(imgIndicies)
            EpochMovie = np.array([ pycdf.lib.tt2000_to_datetime(EpochRocket[0][i]) for i in range(len(EpochRocket[0])) ])

            ###############################
            # --- EXTEND LOW FLYER DATA ---
            ###############################

            # --- extend Low Flyer Rocket data to be the same length as High flyer in the beginning and end ---
            highFlyerSize, lowFlyerSize = len(geoAlt[0]), len(geoAlt[1])

            # --- Append start Values to  ---
            no_of_points_start = int((EpochRocket[1][0] - EpochRocket[0][0]) / (rocketAttrs.MinorFrameTime))
            newAlt = [geoAlt[1][0] for i in range(no_of_points_start)]
            newLat = [geoLat[1][0] for i in range(no_of_points_start)]
            newLong = [geoLong[1][0] for i in range(no_of_points_start)]


            for i in range(len(geoAlt[1])):
                newAlt.append(geoAlt[1][i])
                newLat.append(geoLat[1][i])
                newLong.append(geoLong[1][i])

            # --- Append the ending values ---
            remainingIndicies = highFlyerSize - (lowFlyerSize + no_of_points_start)

            for i in range(remainingIndicies):
                newAlt.append(geoAlt[1][-1])
                newLat.append(geoLat[1][-1])
                newLong.append(geoLong[1][-1])

            geoAlt[1], geoLat[1], geoLong[1] = np.array(newAlt), np.array(newLat), np.array(newLong)

            ##########################
            # --- INITIALIZE PLOTS ---
            ##########################
            prgMsg('Initializing Allsky Plots')

            # --- prepare map information ---
            projPC = ccrs.PlateCarree()  # MUST be kept on the set_extent, crs =, command AND pcolormesh transform command
            lonW =8
            lonE = 22
            latS = 68
            latN = 74.5
            cLat = (latN + latS) / 2
            cLon = (lonW + lonE) / 2
            res = '10m'
            # projType = ccrs.Mollweide(central_longitude=cLon, central_latitude=cLat)
            projType = ccrs.Stereographic(central_longitude=cLon, central_latitude=cLat)

            allSkyCmap = 'inferno'
            LowFlyerColor = 'royalblue'
            HighFlyerColor = 'red'
            LowFlyerProjectionColor = 'cyan'
            HighFlyerProjectionColor = 'darkorange'

            AltvsLatLineWidth = 4
            LatvsLongLineWidth = 4
            AllSkyLineWidth = 4

            # ---------------------
            # --- START PLOTING ---
            # ---------------------

            # figure size
            figure_height = 20
            figure_width = 36
            plt.style.use('dark_background')
            fig = plt.figure(dpi=60)
            fig.set_figwidth(figure_width)
            fig.set_figheight(figure_height)

            # The overall plot shape
            nRows = 20
            nCols = 51 # MUST BE EVEN VALUE
            smallPlotSize = 4
            verticalSpacing = 3 #spacing between the top row and bottom row plots
            gs = gridspec.GridSpec(nrows=nRows, ncols=nCols, figure=fig)

            # define the axes
            if removeLatLongPlot:
                axAlt = fig.add_subplot(gs[0:smallPlotSize, 0:51])
            else:
                axAlt = fig.add_subplot(gs[0:smallPlotSize, 0:23+1])
                axLatLong = fig.add_subplot(gs[0:smallPlotSize, 27:51])

            ax5577 = fig.add_subplot(gs[smallPlotSize+verticalSpacing:nRows, 0:22+1], projection=projType)
            ax5577cbar = fig.add_subplot(gs[smallPlotSize+verticalSpacing:nRows, 23])
            ax6300 = fig.add_subplot(gs[smallPlotSize+verticalSpacing:nRows, 27:50], projection=projType)
            ax6300cbar = fig.add_subplot(gs[smallPlotSize+verticalSpacing:nRows, 50])

            # --- initialize the title ---
            supTitle = fig.suptitle(f'ACESII\n'
                                    f'{EpochMovie[0]}\n'
                                    f'Projection Altitude: {projectionAltitude} km', fontsize=40,color='white')

            # --- plot the norwegian map data ---
            # 5577A
            ax5577.set_title('557.7 nm - 150km', fontsize=30)
            gl5577 = ax5577.gridlines(draw_labels=True, linewidth=3, alpha=0.4, linestyle='--')
            gl5577.xlabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
            gl5577.ylabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
            ax5577.set_extent([lonW, lonE, latS, latN], crs=projPC)  # controls lat/long axes display
            ax5577.coastlines(resolution=res, color='white', alpha=0.8)  # adds coastlines with resolution
            ax5577.set_aspect(0.7)

            # 6300A
            ax6300.set_title('630.0 nm - 250km', fontsize=30)
            gl6300 = ax6300.gridlines(draw_labels=True, linewidth=3, color='white', alpha=0.4, linestyle='--')
            gl6300.xlabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
            gl6300.ylabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
            ax6300.set_extent([lonW, lonE, latS, latN], crs=projPC)  # controls lat/long axes display
            ax6300.coastlines(resolution=res, color='white', alpha=0.8)  # adds coastlines with resolution
            ax6300.set_aspect(0.7)

            # --- initialize the All Sky Image ---
            cbarVmin = 0
            cbarVmax = 16000
            # 5577A
            cmap5577 = ax5577.pcolormesh(allGLongs[0], allGlats[0], allImages[0][0], cmap=allSkyCmap, transform=projPC,vmin=cbarVmin,vmax=cbarVmax)
            plt.colorbar(mappable=cmap5577, cax=ax5577cbar, orientation='vertical',fraction=0.046, pad=0.04)
            ax5577cbar.set_ylabel('Brightness (Rayleigh)',fontsize=20)
            ax5577cbar.tick_params(labelsize=20)

            # 6300A
            cmap6300 = ax6300.pcolormesh(allGLongs[1], allGlats[1], allImages[1][0], cmap=allSkyCmap, transform=projPC,vmin=cbarVmin,vmax=cbarVmax)
            plt.colorbar(mappable=cmap6300, cax=ax6300cbar, orientation='vertical',fraction=0.046, pad=0.04)
            ax6300cbar.set_ylabel('Brightness (Rayleigh)',fontsize=20)
            ax6300cbar.tick_params(labelsize=20)

            # --- initialize total trajectory on the All Sky Image ---

            # plot projected B on 5577 Allsky
            ax5577.plot(geoILong[1], geoILat[1], color=LowFlyerProjectionColor, transform=projPC, linewidth=AllSkyLineWidth,label='Low Flyer Projected IGRF')
            ax5577.plot(geoILong[0], geoILat[0], color=HighFlyerProjectionColor, transform=projPC, linewidth=AllSkyLineWidth,label='High Flyer Projected IGRF') if not useOnlyLowFlyer else [0,]

            # plot Ilat/Ilong marker on 5577 on allsky
            AllSky5577_projected_marker_low = ax5577.scatter(geoILong[1][0], geoILat[1][0], color=LowFlyerProjectionColor, marker='x', linewidth=10, transform=projPC)
            AllSky5577_projected_marker_high = ax5577.scatter(geoILong[0][0], geoILat[0][0], color=HighFlyerProjectionColor, marker='x', linewidth=10, transform=projPC) if not useOnlyLowFlyer else [0,]

            # plot on 5577 Allsky
            ax5577.plot(geoLong[1], geoLat[1], color=LowFlyerColor,linewidth=4, transform=projPC, label='Low Flyer Trajectory')
            ax5577.plot(geoLong[0], geoLat[0], color=HighFlyerColor,linewidth=4, transform=projPC,label='High Flyer Trajectory') if not useOnlyLowFlyer else [0,]

            # plot marker on 5577 Allsky
            AllSky5577_marker_low = ax5577.scatter(geoLong[1][0], geoLat[1][0], color=LowFlyerColor, marker='x', linewidth=10, transform=projPC)
            AllSky5577_marker_high = ax5577.scatter(geoLong[0][0], geoLat[0][0], color=HighFlyerColor, marker='x', linewidth=10, transform=projPC) if not useOnlyLowFlyer else [0,]

            # plot projected B on 6300 Allsky
            ax6300.plot(geoILong[1], geoILat[1], color=LowFlyerProjectionColor, transform=projPC, linewidth=AllSkyLineWidth)
            ax6300.plot(geoILong[0], geoILat[0], color=HighFlyerProjectionColor, transform=projPC, linewidth=AllSkyLineWidth) if not useOnlyLowFlyer else [0,]

            # plot Ilat/Ilong marker on 6300 on allsky
            AllSky6300_projected_marker_low = ax6300.scatter(geoILong[1][0], geoILat[1][0], color=LowFlyerProjectionColor, marker='x', linewidth=10, transform=projPC)
            AllSky6300_projected_marker_high = ax6300.scatter(geoILong[0][0], geoILat[0][0], color=HighFlyerProjectionColor, marker='x', linewidth=10, transform=projPC) if not useOnlyLowFlyer else [0,]

            # plot on 6300 Allsky
            ax6300.plot(geoLong[1], geoLat[1], color=LowFlyerColor,linewidth=4, transform=projPC)
            ax6300.plot(geoLong[0], geoLat[0], color=HighFlyerColor,linewidth=4, transform=projPC) if not useOnlyLowFlyer else [0,]

            # plot marker on 6300 Allsky
            AllSky6300_marker_low = ax6300.scatter(geoLong[1][0], geoLat[1][0], color=LowFlyerColor, linewidth=10, marker='x',transform=projPC)
            AllSky6300_marker_high = ax6300.scatter(geoLong[0][0], geoLat[0][0], color=HighFlyerColor, linewidth=10, marker='x',transform=projPC) if not useOnlyLowFlyer else [0,]

            # -----------------------------------------------------------------------------------
            # stylization parameters of text labels in alt vs lat and lat vs long plots
            textOffset_alt = 0.15
            textOffset_latlong = [[-0.15, -0.15], [0.15, 0.25]]
            text_alignment = ['left', 'right']
            textUTC_style_lat = [dict(size=20, color=HighFlyerColor), dict(size=20, color=LowFlyerColor)]
            textUTC_style_project_lat = [dict(size=20, color=HighFlyerProjectionColor), dict(size=20, color=LowFlyerProjectionColor)]
            textUTC_style_alt = [dict(size=20, color=HighFlyerColor), dict(size=20, color=LowFlyerColor)]
            altrounding = 1  # how many decimals to round to
            # -----------------------------------------------------------------------------------

            # --- initialize Altitude vs Lat Plot ---
            axAlt.set_ylabel('Altitude [km]', fontsize=20)
            axAlt.set_xlabel('Geo Lat', fontsize=20)
            axAlt.tick_params(axis='x', which='both', labelsize=25)
            axAlt.tick_params(axis='y', which='both', labelsize=25)
            axAlt.axhline(projectionAltitude, color='limegreen', linestyle='--', alpha=0.5,linewidth=AltvsLatLineWidth)
            axAlt.text(69.05, 120, '100 km', color='limegreen', fontsize=20)

            # trajectory
            axAlt.plot(geoLat[0], geoAlt[0], color=HighFlyerColor, linewidth=AltvsLatLineWidth) if not useOnlyLowFlyer else [0,]
            axAlt.plot(geoLat[1], geoAlt[1], color=LowFlyerColor, linewidth=AltvsLatLineWidth)

            # B-projection
            if useRealMagData:
                legendAltLatLabel = f'Flight B-Field Data {projectionAltitude} km'
            else:
                legendAltLatLabel = f'Projected IGRF {projectionAltitude} km'
            BprojectionHigh, = axAlt.plot(projectB[0][0][0], projectB[0][0][1], color='white', alpha=0.7, linestyle='--') if not useOnlyLowFlyer else [0,]  # high flyer
            BprojectionLow, = axAlt.plot(projectB[1][0][0], projectB[1][0][1], color='white', alpha=0.7, linestyle='--',label=legendAltLatLabel)  # low flyer

            # legend for AltvsLat projection
            axAlt.legend(loc="upper left",fontsize=15)

            # marker
            Altlat_marker_high = axAlt.scatter(geoLat[0][0], geoAlt[0][0], color=HighFlyerColor, marker='x', linewidth=15) if not useOnlyLowFlyer else [0,]
            Altlat_marker_low = axAlt.scatter(geoLat[1][0], geoAlt[1][0], color=LowFlyerColor, marker='x', linewidth=15)

            # text labels
            Altlat_text_high = axAlt.text(geoLat[0][0] - textOffset_alt, geoAlt[0][0],
                                       f'{round(geoAlt[0][0], altrounding)} km', ha=text_alignment[1],
                                       **textUTC_style_alt[0]) if not useOnlyLowFlyer else [0,]
            Altlat_text_low = axAlt.text(geoLat[1][0] - textOffset_alt, geoAlt[1][0],
                                      f'{round(geoAlt[1][0], altrounding)} km', ha=text_alignment[1],
                                      **textUTC_style_alt[1])



            # --- intialize Trajectory ILat vs Ilong projection ---
            if not removeLatLongPlot:

                # full trajectory
                axLatLong.plot(geoILong[0], geoILat[0], color=HighFlyerProjectionColor, label=f'Projected IGRF ({projectionAltitude} km)',linewidth=LatvsLongLineWidth) if not useOnlyLowFlyer else [0,]
                axLatLong.plot(geoILong[1], geoILat[1], color=LowFlyerProjectionColor, label=f'Projected IGRF ({projectionAltitude} km)',linewidth=LatvsLongLineWidth)

                # marker
                IlatILong_marker_high = axLatLong.scatter(geoILong[0][0], geoILat[0][0], color=HighFlyerProjectionColor, marker='x',linewidth=10) if not useOnlyLowFlyer else [0,]
                IlatILong_marker_low = axLatLong.scatter(geoILong[1][0], geoILat[1][0], color=LowFlyerProjectionColor, marker='x',linewidth=10)

                # text labels
                # IlatILong_text_high = axLatLong.text(geoILong[0][0] - textOffset_latlong[0][0],
                #                                            geoILat[0][0] - textOffset_latlong[0][1],
                #                                            f'({round(geoILong[0][0], altrounding)}' + '$^{\circ}$' + f', {round(geoILat[0][0], altrounding)}' + '$^{\circ}$)',
                #                                            ha=text_alignment[0], **textUTC_style_project_lat[0])
                # IlatILong_text_low = axLatLong.text(geoILong[1][0] - textOffset_latlong[1][0],
                #                                           geoILat[1][0] - textOffset_latlong[1][1],
                #                                           f'({round(geoILong[1][0], altrounding)}' + '$^{\circ}$' + f', {round(geoILat[1][0], altrounding)}' + '$^{\circ}$)',
                #                                           ha=text_alignment[1], **textUTC_style_project_lat[1])

                IlatILong_text_high = axLatLong.text(14.65,
                                                     70.2,
                                                     f'({round(geoILong[0][0], altrounding)}' + '$^{\circ}$' + f', {round(geoILat[0][0], altrounding)}' + '$^{\circ}$)',
                                                     ha='center', **textUTC_style_project_lat[0]) if not useOnlyLowFlyer else [0,]
                IlatILong_text_low = axLatLong.text(15.1,
                                                    70.2,
                                                    f'({round(geoILong[1][0], altrounding)}' + '$^{\circ}$' + f', {round(geoILat[1][0], altrounding)}' + '$^{\circ}$)',
                                                    ha='center', **textUTC_style_project_lat[1])

                # --- initialize Longitude vs Latitude Plot ---
                axLatLong.set_xlim(12.5, 17.5)
                axLatLong.set_ylim(69, 74)
                axLatLong.set_ylabel('Geo Lat', fontsize=30)
                axLatLong.set_xlabel('Geo Long', fontsize=30)
                axLatLong.tick_params(axis='x', which='both', labelsize=25)
                axLatLong.tick_params(axis='y', which='both', labelsize=25)

                # trajectory
                axLatLong.plot(geoLong[0], geoLat[0], color=HighFlyerColor,linewidth=LatvsLongLineWidth, label='High Trajectory') if not useOnlyLowFlyer else [0,]
                axLatLong.plot(geoLong[1], geoLat[1], color=LowFlyerColor,linewidth=LatvsLongLineWidth, label='Low Trajectory')

                # marker
                latLong_marker_high = axLatLong.scatter(geoLong[0][0], geoLat[0][0], color=HighFlyerColor, marker='x', linewidth=LatvsLongLineWidth) if not useOnlyLowFlyer else [0,]
                latLong_marker_low = axLatLong.scatter(geoLong[1][0], geoLat[1][0], color=LowFlyerColor, marker='x', linewidth=LatvsLongLineWidth)

                # text labels
                # latLong_text_high = axLatLong.text(geoLong[0][0] - textOffset_latlong[0][0],
                #                                    geoLat[0][0] - textOffset_latlong[0][1],
                #                                    f'({round(geoLong[0][0], altrounding)}' + '$^{\circ}$' + f', {round(geoLat[0][0], altrounding)}' + '$^{\circ}$)',
                #                                    ha=text_alignment[0], **textUTC_style_lat[0])
                # latLong_text_low = axLatLong.text(geoLong[1][0] - textOffset_latlong[1][0],
                #                                   geoLat[1][0] - textOffset_latlong[1][1],
                #                                   f'({round(geoLong[1][0], altrounding)}' + '$^{\circ}$' + f', {round(geoLat[1][0], altrounding)}' + '$^{\circ}$)',
                #                                   ha=text_alignment[1], **textUTC_style_lat[1])

                latLong_text_high = axLatLong.text(14.65,
                                                   69.5,
                                                   f'({round(geoLong[0][0], altrounding)}' + '$^{\circ}$' + f', {round(geoLat[0][0], altrounding)}' + '$^{\circ}$)',
                                                   ha='center', **textUTC_style_lat[0]) if not useOnlyLowFlyer else [0,]
                latLong_text_low = axLatLong.text(15.1,
                                                  69.5,
                                                  f'({round(geoLong[1][0], altrounding)}' + '$^{\circ}$' + f', {round(geoLat[1][0], altrounding)}' + '$^{\circ}$)',
                                                  ha='center', **textUTC_style_lat[1])

                # legend for LatvsLong projection
                axLatLong.legend(loc="lower left", fontsize=17)
            else:
                # plot the lat/long coordinates that update on the 5577 graph

                IlatILong_text_low = ax5577.text(8.4, 73.25, f'({round(geoILong[1][0], altrounding)}' + '$^{\circ}$' + f', {round(geoILat[1][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_project_lat[1], transform=projPC)
                IlatILong_text_high = ax5577.text(8.5, 73, f'({round(geoILong[0][0], altrounding)}' + '$^{\circ}$' + f', {round(geoILat[0][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_project_lat[0],transform=projPC) if not useOnlyLowFlyer else [0, ]
                latLong_text_low = ax5577.text(8.6, 72.75, f'({round(geoLong[1][0], altrounding)}' + '$^{\circ}$' + f', {round(geoLat[1][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_lat[1], transform=projPC)
                latLong_text_high = ax5577.text(8.7, 72.5, f'({round(geoLong[0][0], altrounding)}' + '$^{\circ}$' + f', {round(geoLat[0][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_lat[0],transform=projPC) if not useOnlyLowFlyer else [0, ]

                ax5577.legend(loc='upper left',prop={'size': 14})

            Done(start_time)


            ###########################
            # --- ANIMATE THE MOVIE ---
            ###########################
            def animatePlot(i):

                # update the title
                supTitle.set_text(f'ACESII\n'
                                  f'{EpochMovie[i]}\n'
                                  f'Projection Alt: {projectionAltitude} km')

                # --- ALL SKY ---
                # update all sky images
                wImage5577 = imgIndicies[0][i]
                wImage6300 = imgIndicies[1][i]

                cmap5577.set_array(allImages[0][imgIndicies[0][i]].ravel())
                cmap6300.set_array(allImages[1][imgIndicies[1][i]].ravel())

                # update the allsky5577 rocket marker positions
                AllSky5577_marker_high.set_offsets([geoLong[0][i], geoLat[0][i]]) if not useOnlyLowFlyer else [0,]
                AllSky5577_marker_low.set_offsets([geoLong[1][i], geoLat[1][i]])
                AllSky5577_projected_marker_high.set_offsets([geoILong[0][i], geoILat[0][i]]) if not useOnlyLowFlyer else [0,]
                AllSky5577_projected_marker_low.set_offsets([geoILong[1][i], geoILat[1][i]])

                # update the allsky6300 rocket marker positions
                AllSky6300_marker_high.set_offsets([geoLong[0][i], geoLat[0][i]]) if not useOnlyLowFlyer else [0,]
                AllSky6300_marker_low.set_offsets([geoLong[1][i], geoLat[1][i]])
                AllSky6300_projected_marker_high.set_offsets([geoILong[0][i], geoILat[0][i]]) if not useOnlyLowFlyer else [0,]
                AllSky6300_projected_marker_low.set_offsets([geoILong[1][i], geoILat[1][i]])

                # --- ALT VS LAT ---
                # update Alt vs lat marker
                Altlat_marker_high.set_offsets([geoLat[0][i], geoAlt[0][i]]) if not useOnlyLowFlyer else [0,]
                Altlat_marker_low.set_offsets([geoLat[1][i], geoAlt[1][i]])

                # update B-project lines on Alt vs lat plot
                BprojectionHigh.set_xdata(projectB[0][i][0]) if not useOnlyLowFlyer else [0,]
                BprojectionHigh.set_ydata(projectB[0][i][1]) if not useOnlyLowFlyer else [0,]
                BprojectionLow.set_xdata(projectB[1][i][0])
                BprojectionLow.set_ydata(projectB[1][i][1])

                # update Alt vs lat text
                Altlat_text_high.set_x(geoLat[0][i]) if not useOnlyLowFlyer else [0,]
                Altlat_text_high.set_y(geoAlt[0][i]) if not useOnlyLowFlyer else [0,]
                Altlat_text_high.set_text(f'{round(geoAlt[0][i], altrounding)} km') if not useOnlyLowFlyer else [0,]

                Altlat_text_low.set_x(geoLat[1][i])
                Altlat_text_low.set_y(geoAlt[1][i])
                Altlat_text_low.set_text(f'{round(geoAlt[1][i], altrounding)} km')

                # --- LAT VS LONG ---
                if not removeLatLongPlot:

                    # --- update marker ---
                    # only this code is in the "removeLatLongPlot" since if we remove the latlong plot these
                    # same markers go to the 5577 plot to be updated
                    latLong_marker_high.set_offsets([geoLong[0][i], geoLat[0][i]]) if not useOnlyLowFlyer else [0,]
                    latLong_marker_low.set_offsets([geoLong[1][i], geoLat[1][i]])

                    # update Ilat vs Ilong marker
                    IlatILong_marker_high.set_offsets([geoILong[0][i], geoILat[0][i]]) if not useOnlyLowFlyer else [0, ]
                    IlatILong_marker_low.set_offsets([geoILong[1][i], geoILat[1][i]])

                # update lat vs long text
                latLong_text_high.set_text(f'({round(geoLong[0][i], altrounding)}' + '$^{\circ}$' +
                                           f', {round(geoLat[0][i], altrounding)}' + '$^{\circ}$)') if not useOnlyLowFlyer else [0,]
                latLong_text_low.set_text(f'({round(geoLong[1][i], altrounding)}' + '$^{\circ}$' +
                                           f', {round(geoLat[1][i], altrounding)}' + '$^{\circ}$)')

                # update Ilat vs Ilong text
                IlatILong_text_high.set_text(f'({round(geoILong[0][i], altrounding)}' + '$^{\circ}$' +
                                           f', {round(geoILat[0][i], altrounding)}' + '$^{\circ}$)') if not useOnlyLowFlyer else [0,]

                IlatILong_text_low.set_text(f'({round(geoILong[1][i], altrounding)}' + '$^{\circ}$' +
                                          f', {round(geoILat[1][i], altrounding)}' + '$^{\circ}$)')


            prgMsg('Creating AllSky Movie')

            if plotSpecificLocations:
                locations = [i for i in specificLocations]
            elif doFrameSkips:
                locations = [i for i in range(0, len(EpochMovie), frame_skips)]  # NEEDS TO BE THE HIGH FLYER LENGTH
            else:
                locations = [i for i in range(len(EpochMovie))]


            anim = animation.FuncAnimation(fig=fig, func=animatePlot, interval=1000 / fps, frames=locations)

            # need a .mp4 writer, following code points to it
            writervideo = animation.FFMpegWriter(fps=fps)
            import matplotlib as mpl
            mpl.rcParams['animation.ffmpeg_path'] = r'C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\supportCode\ffmpeg\bin\ffmpeg.exe'
            anim.save(r'C:\Data\ACESII\trajectories\trajectory_plots\movies\ACESII_AllSky.mp4',writer=writervideo)


            Done(start_time)








# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder

if justPrintSiteNames:
    allSkyIDL_to_py(wSite, justPrintSiteNames,rocketFolderPath)
else:
    allSkyIDL_to_py(wSite, justPrintSiteNames,rocketFolderPath)