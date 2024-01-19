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
import cartopy.crs as ccrs
from netCDF4 import Dataset as netcdf_dataset


print(color.UNDERLINE + f'Plot1_AllSky' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# -------------GENERAL PLOT TOGGLES-------------------
cmapColor = 'viridis'
faceColorChoice = (156 / 255, 156 / 255, 156 / 255, 0.5)  # in normalize RGBA
# --------------ALTvsLAT------------------
altLatPlot = True
altLat_wScatterPoints = np.linspace(200, 400, 6)
trajColors = ['tab:red', 'tab:orange']
altLat_labelFontSize = 60
altLat_textSize = 45
altLat_tickSize = 55
altLat_scatterSize = 400
altLat_lineThickness = 8
# ---------------BigAllSky-----------------
BigAllSkyPlot = True
BigAllSky_wScatterPoints = np.linspace(150, 500, 7)
# lonW = 10
# lonE = 23.5
lonW = 11.5
lonE = 18.5
latS = 68
latN = 75
res = '10m'
wImage = 10
cbarVmin,cbarVmax = 0, 16 # in kRayleigh
BigAllSky_textSize = 35
BigAllSky_tickSize = 40
BigAllSky_scatterSize = 425
BigAllSky_lineThickness = 8
BigAllSky_GridSize = 5
BigAllSky_TitleSize = 50
# --------------------------------
SmallAllSkyPlot = False
lonW = 11.5
lonE = 18.5
latS = 68
latN = 75
res = '10m'
# --------------------------------
makeColorbarPlot = True


# --- --- --- --- --- --- -
# --- LOAD ALL THE DATA ---
# --- --- --- --- --- --- -

# trajectory
trajFolderPath = f'{ACES_data_folder}trajectories\\'
inputFilesTraj = [glob(trajFolderPath + rf'{fliers[0]}\\*_ILat_ILong*')[0],
                  glob(trajFolderPath + rf'{fliers[1]}\\\\*_ILat_ILong*')[0]]

# --- GET TRAJECTORY DATA ---
prgMsg(f'Loading ACESII traj data')
data_dicts_traj = []
for i in range(2):
    data_dict_traj = loadDictFromFile(inputFilesTraj[i], input_data_dict={}, reduceData=False, targetTimes=[], wKeys=[])
    data_dict_traj['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch_esa'][0][i]) for i in (range(len(data_dict_traj['Epoch_esa'][0])))])
    data_dicts_traj.append(data_dict_traj)

# define some variables
EpochRocket = [data_dicts_traj[0]['Epoch_esa'][0], data_dicts_traj[1]['Epoch_esa'][0]]
geoAlt = [data_dicts_traj[0]['geoAlt'][0], data_dicts_traj[1]['geoAlt'][0]]
geoLat = [data_dicts_traj[0]['geoLat'][0], data_dicts_traj[1]['geoLat'][0]]
geoLong = [data_dicts_traj[0]['geoLong'][0], data_dicts_traj[1]['geoLong'][0]]
geoMagLat = [data_dicts_traj[0]['geomagLat'][0], data_dicts_traj[1]['geomagLat'][0]]
Done(start_time)


# Load AllSky data
prgMsg('Loading Allsky Data')
data_dict_allSky5577 = loadDictFromFile(glob(r'C:\Data\ACESII\all_sky\skibotn\5577\\*.cdf')[0])
data_dict_allSky6300 = loadDictFromFile(glob(r'C:\Data\ACESII\all_sky\skibotn\6300\\*.cdf')[0])
Done(start_time)




############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################


# --- ALTITUDE VS LATITUDE PLOT ---
if altLatPlot:
    prgMsg('Plotting AltLat')

    # --- --- --- --- ---
    # --- AltLat plot ---
    # --- --- --- --- ---
    # axAltLat = fig.add_subplot(gs_altLat_BigAllSky[0])
    fig, axAltLat = plt.subplots()
    figure_height = 10
    figure_width = 43
    fig.set_figwidth(figure_width)
    fig.set_figheight(figure_height)

    axAltLat.set_ylabel('Altitude [km]', fontsize=altLat_labelFontSize,weight='bold')
    axAltLat.set_xlabel('Geomagnetic Lat [deg]',fontsize=altLat_labelFontSize,weight='bold')
    axAltLat.set_xlim(67.25, 72.5)
    axAltLat.set_ylim(0, 460)

    # plot the pseudo geomagnetic field line
    slope = -1 * (111 / np.sin(np.radians(90 - 78.13)))  # corresponds to line with -78.13deg inclination
    for i in range(21):
        axAltLat.axline(xy1=(66 + i * 0.5, 0), slope=slope, color='tab:blue', linewidth=altLat_lineThickness, linestyle='-.', alpha=0.3)
    axAltLat.legend(['B$_{Geo}$'], loc='upper right',fontsize=65)

    # set the facecolor of the axAltLat plot
    axAltLat.set_facecolor(faceColorChoice)

    # plot the UTC labels
    axGeographicLat = axAltLat.twiny()
    axGeographicLat.plot(geoLat[0], geoAlt[0], color=trajColors[0], alpha=0)  # High
    axGeographicLat.plot(geoLat[1], geoAlt[1], color=trajColors[1], alpha=0)  # Low
    axGeographicLat.set_xlabel('Geographic Lat [deg]',fontsize=altLat_labelFontSize,weight='bold')
    timeTargetsUTC_labels = [int(num) for num in altLat_wScatterPoints]
    AltLat_vertical_Alignments = ['bottom' for tme in timeTargetsUTC_labels]
    AltLat_horizontal_Alignments = ['right', 'right', 'right', 'center', 'left', 'left']
    vertical_text_label_adjustments = [-0.05, 0.01, 0.025, 0.025, 0.01, -0.03]
    horizontal_text_label_adjustments = [-0.002, -0.001, 0.000, 0.0, -0.001, 0.001]

    # plot the scatterpoint of each of the timeTargetUTC_labels. Plot the text itself only for the High Flyer and
    # create a connecting line between the scatterpoints between the flyers
    for i in range(2):  # for each rocket
        for j, ttme in enumerate(timeTargetsUTC_labels):
            Index = np.abs(EpochRocket[i] - (ttme * 1E9 + EpochRocket[i][0])).argmin()
            xPos = geoMagLat[i][Index]
            yPos = geoAlt[i][Index]

            if i == 0:
                # plot the text itself
                label = pycdf.lib.tt2000_to_datetime(EpochRocket[i][Index])
                deltaY = vertical_text_label_adjustments[j] * yPos
                deltaX = horizontal_text_label_adjustments[j] * xPos
                axAltLat.text(x=xPos + deltaX, y=yPos + deltaY, s=label.strftime("%H:%M:%S"), color=trajColors[i],
                              va=AltLat_vertical_Alignments[j], ha=AltLat_horizontal_Alignments[j], size=altLat_textSize)

                # plot the connecting line
                Index_LF = Index = np.abs(EpochRocket[1] - (ttme * 1E9 + EpochRocket[1][0])).argmin()
                xPos_LF = geoMagLat[1][Index_LF]
                yPos_LF = geoAlt[1][Index_LF]
                axAltLat.plot([xPos, xPos_LF], [yPos, yPos_LF], color='green', linestyle='--', alpha=0.5,linewidth=altLat_lineThickness)
                # axAltLat.axline(xy1=(xPos, yPos), xy2=(xPos_LF, yPos_LF), color='green', linestyle='--', alpha=0.5)

            # plot a dot at the text label
            axAltLat.scatter(x=xPos, y=yPos, s=altLat_scatterSize, marker="o", color=trajColors[i])

    # adjust the tick label size
    axAltLat.tick_params(axis='both',labelsize=altLat_tickSize)
    axGeographicLat.tick_params(axis='both',labelsize=altLat_tickSize)


    # plot the trajectory over everything
    axAltLat.plot(geoMagLat[0], geoAlt[0], color=trajColors[0], label='High Flyer',linewidth=altLat_lineThickness)  # High
    axAltLat.plot(geoMagLat[1], geoAlt[1], color=trajColors[1], label='Low Flyer',linewidth=altLat_lineThickness)  # Low
    plt.tight_layout()
    plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Paper_Photos\Plot1\\AltLat.png')
    Done(start_time)
    # plt.show()


# --- BIG ALLSKYIMAGER PLOT ---
if BigAllSkyPlot:
    prgMsg('Plotting BigAllSky')

    # --- --- --- --- --- --
    # --- BigAllSky plot ---
    # --- --- --- --- --- --
    for i in range(2):
        # --- PLOT MAP OF NORWAY ---
        projTransform = ccrs.PlateCarree()
        fig, axBigAllSky = plt.subplots(1,subplot_kw=dict(projection=projTransform))
        figure_height = 20
        figure_width = 20
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)

        # get the elevation map data and plot it
        fname = glob(r"C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\Papers\ACESII_Alfvenic_Observations\ElevationData\*.nc*")[0]
        dataset = netcdf_dataset(fname) # Load data into separate arrays,if variables are not known, print(dataset.variables) to check them
        elev = dataset.variables['elevation'][:]
        lats = dataset.variables['lat'][:]
        lons = dataset.variables['lon'][:]
        vmin, vmax = -8000, 3000
        v = np.linspace(vmin, vmax, 100, endpoint=True) # Set how many contour lines to display
        axBigAllSky.contourf(lons, lats, elev, v, cmap="gray", vmin=vmin, vmax=vmax, transform=projTransform)

        # gridlines
        gl = axBigAllSky.gridlines(draw_labels=True, linewidth=BigAllSky_GridSize,
                                   alpha=0.4,
                                   linestyle='--',
                                   color='black')
        gl.xlabel_style = {'size': BigAllSky_tickSize, 'color': 'black', 'weight': 'bold'}
        gl.ylabel_style = {'size': BigAllSky_tickSize, 'color': 'black', 'weight': 'bold'}
        gl.top_labels = False

        # extent of map
        axBigAllSky.set_extent([lonW, lonE, latS, latN])  # controls lat/long axes display

        # coastlines
        axBigAllSky.coastlines(resolution=res, color='black',  alpha=1,linewidth=3)  # adds coastlines with resolution

        if i == 0:
            #--- Plot the Big AllSky image ---
            cmapBigAllSky = axBigAllSky.pcolormesh(data_dict_allSky5577['GLongs'][0], data_dict_allSky5577['GLats'][0], data_dict_allSky5577['AllSkyImages'][0][wImage],
                                                   cmap=cmapColor,
                                                   transform=projTransform,
                                                   vmin=cbarVmin,
                                                   vmax=cbarVmax,
                                                   alpha=1)
            BigAllSky_outputPath = r'C:\Users\cfelt\OneDrive\Desktop\Paper_Photos\Plot1\\BigAllSky_5570.png'
            fig.suptitle('Skibton 5577A - 150 km\n' + data_dict_allSky5577['Epoch'][0][wImage].strftime("%Y-%B-%d %H:%M:%S") + ' UTC',fontsize=BigAllSky_TitleSize,weight='bold')
        elif i == 1:
            cmapBigAllSky = axBigAllSky.pcolormesh(data_dict_allSky6300['GLongs'][0], data_dict_allSky6300['GLats'][0],
                                                   data_dict_allSky6300['AllSkyImages'][0][wImage],
                                                   cmap=cmapColor,
                                                   transform=projTransform,
                                                   vmin=cbarVmin,
                                                   vmax=cbarVmax,
                                                   alpha=1)
            BigAllSky_outputPath = r'C:\Users\cfelt\OneDrive\Desktop\Paper_Photos\Plot1\\BigAllSky_6300.png'
            fig.suptitle('Skibton 6300A - 250 km\n' + data_dict_allSky6300['Epoch'][0][wImage].strftime("%Y-%B-%d %H:%M:%S")+ ' UTC',fontsize=BigAllSky_TitleSize,weight='bold')
        axBigAllSky.set_facecolor(faceColorChoice)

        # --- plot the rocket trajectory data on the large AllSky plot ---
        axBigAllSky.plot(geoLong[0], geoLat[0], color=trajColors[0], transform=projTransform,linewidth=BigAllSky_lineThickness) # High
        axBigAllSky.plot(geoLong[1], geoLat[1], color=trajColors[1], transform=projTransform,linewidth=BigAllSky_lineThickness) # Low

        # plot specific UTC times on the trajectory lines
        timeTargetsUTC_labels = [int(tme*1E9 + pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 000))) for tme in BigAllSky_wScatterPoints]

        alignment = ['left', 'right']
        for i in range(2): # for each rocket
            for ttme in timeTargetsUTC_labels:

                Index = np.abs(EpochRocket[i] - ttme).argmin()
                label = pycdf.lib.tt2000_to_datetime(EpochRocket[i][Index])

                xPos = geoLong[i][Index]
                yPos = geoLat[i][Index]

                # plot a dot at the text label
                axBigAllSky.scatter(x=xPos,y=yPos,s=BigAllSky_scatterSize,marker="o",color=trajColors[i],transform=projTransform)

                # plot the text itself
                deltaX = 0.0075*xPos if i == 0 else -1*0.0075*xPos
                axBigAllSky.text(x=xPos + deltaX, y=yPos, s =label.strftime("%H:%M:%S"),color=trajColors[i], ha=alignment[i], transform=projTransform,size=BigAllSky_textSize)

        plt.tight_layout()
        plt.savefig(BigAllSky_outputPath)
        # plt.show()
        Done(start_time)

# --- --- --- --- --- --
# --- 5 AllSky Plots ---
# --- --- --- --- --- --
if SmallAllSkyPlot:
    prgMsg('Creating Smaller Plots')
    imageDicts = [data_dict_allSky5577,data_dict_allSky6300]
    wColor = ['5577','6300']

    for j, data_dict in enumerate(imageDicts):
        for i, tme in enumerate(data_dict['Epoch'][0]):

            # --- PLOT MAP OF NORWAY ---
            projTransform = ccrs.PlateCarree()
            fig, axSmallAllSky = plt.subplots(1,subplot_kw=dict(projection=projTransform))
            figure_height = 20
            figure_width = 20
            fig.set_figwidth(figure_width)
            fig.set_figheight(figure_height)

            # get the elevation map data and plot it
            # fname = glob(r"C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\Papers\ACESII_Alfvenic_Observations\ElevationData\*.nc*")[0]
            # dataset = netcdf_dataset(fname) # Load data into separate arrays,if variables are not known, print(dataset.variables) to check them
            # elev = dataset.variables['elevation'][:]
            # lats = dataset.variables['lat'][:]
            # lons = dataset.variables['lon'][:]
            # vmin, vmax = -6000, 3000
            # v = np.linspace(vmin, vmax, 100, endpoint=True) # Set how many contour lines to display
            # axSmallAllSky.contourf(lons, lats, elev, v, cmap="gray", vmin=vmin, vmax=vmax, transform=projTransform)

            # gridlines
            gl = axSmallAllSky.gridlines(draw_labels=True, linewidth=BigAllSky_GridSize,
                                       alpha=0.1,
                                       linestyle='--',
                                       color='black')
            gl.xlabel_style = {'size': BigAllSky_tickSize, 'color': 'black', 'weight': 'bold'}
            gl.ylabel_style = {'size': BigAllSky_tickSize, 'color': 'black', 'weight': 'bold'}
            gl.top_labels = False

            # extent of map
            axSmallAllSky.set_extent([lonW, lonE, latS, latN])  # controls lat/long axes display

            # coastlines
            axSmallAllSky.coastlines(resolution=res, color='black',  alpha=1,linewidth=3)  # adds coastlines with resolution

            #--- Plot the Big AllSky image ---
            cmapBigAllSky = axSmallAllSky.pcolormesh(data_dict['GLongs'][0], data_dict['GLats'][0], data_dict['AllSkyImages'][0][i],
                                                   cmap=cmapColor,
                                                   transform=projTransform,
                                                   vmin=cbarVmin,
                                                   vmax=cbarVmax,
                                                   alpha=1)
            BigAllSky_outputPath = rf'C:\Users\cfelt\OneDrive\Desktop\Paper_Photos\Plot1\\SmallAllSky\\smallAllSky_{wColor[j]}A_{tme.strftime("%H%M%S")}.png'
            fig.suptitle(data_dict['Epoch'][0][i].strftime("%Y-%B-%d %H:%M")+ ' UTC',fontsize=BigAllSky_TitleSize,weight='bold')
            axSmallAllSky.set_facecolor(faceColorChoice)
            plt.tight_layout()
            plt.savefig(BigAllSky_outputPath)
            # plt.show()
            Done(start_time)

    Done(start_time)

# --- --- --- ----
# --- COLORBAR ---
# --- --- --- ----
if makeColorbarPlot:

    # --- PLOT MAP OF NORWAY ---
    projTransform = ccrs.PlateCarree()
    fig, axColorbarPlot = plt.subplots(1, subplot_kw=dict(projection=projTransform))
    figure_height = 8
    figure_width = 40
    fig.set_figwidth(figure_width)
    fig.set_figheight(figure_height)

    # gridlines
    gl = axColorbarPlot.gridlines(draw_labels=True, linewidth=BigAllSky_GridSize,
                               alpha=0.0,
                               linestyle='--',
                               color='black')
    gl.xlabel_style = {'size': BigAllSky_tickSize, 'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'size': BigAllSky_tickSize, 'color': 'black', 'weight': 'bold'}
    gl.top_labels = False
    gl.bottom_labels = False

    # extent of map
    axColorbarPlot.set_extent([lonW, lonE, latS, latN])  # controls lat/long axes display

    # --- Plot the Big AllSky image ---
    cmapBigAllSky = axColorbarPlot.pcolormesh(data_dict_allSky5577['GLongs'][0], data_dict_allSky5577['GLats'][0],
                                           data_dict_allSky5577['AllSkyImages'][0][wImage],
                                           cmap=cmapColor,
                                           transform=projTransform,
                                           vmin=cbarVmin,
                                           vmax=cbarVmax,
                                           alpha=1)

    cbar = plt.colorbar(mappable=cmapBigAllSky, orientation='horizontal',fraction=0.8, pad=2)
    cbar.set_label('Intensity [kR]', fontsize=65, weight='bold')
    cbar.ax.tick_params(labelsize=60)

    # plt.tight_layout()
    plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Paper_Photos\Plot1\colorbar.png')
    # plt.show()
    Done(start_time)
