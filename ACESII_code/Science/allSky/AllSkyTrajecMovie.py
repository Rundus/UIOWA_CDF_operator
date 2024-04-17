# --- AllSkyTrajecMovie.py ---
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
from ACESII_code.myImports import *

start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintSiteNames = False

# --- Select the Site ---
inputPath_modifier_AllSky = 'all_sky\skibotn' # e.g. 'L1' or 'L1'. It's the name of the broader input folder inside data\ACESII
inputPath_modifier_attitude = 'attitude'

# --- AllSky MOVIE ---
createAllSkyMovie = True
fps = 20 # fps of the video
projectionAltitude = [150, 250] # in km. Format: [green, red]

# plot specific locations
plotSpecificLocations = False
specificLocations = [i*10 for i in range(int(12371/10))]

# plot frame skips
doFrameSkips = True
frame_skips = 5 # how many Epoch_esa frames to skip

# Kenton wanted just the low flyer for a plot, so here's a toggle to do that
useOnlyLowFlyer = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import CHAOS
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
from matplotlib import pyplot as plt, animation
from ACESII_code.data_paths import ACES_data_folder, fliers
from my_matplotlib_Assets.colorbars.matlab_parula import matlab_parula_cmap



# --- --- --- ---  --
# --- PLOT PARAMS ---
# --- --- --- ---  --
dpi = 60
figure_height = 10
figure_width = 10
lonW = 8
lonE = 27
latS = 68
latN = 74
cLat = (latN + latS) / 2
cLon = (lonW + lonE) / 2
res = '110m'
projProjection = ccrs.Orthographic(central_longitude=15, central_latitude=70)
projTransform = ccrs.PlateCarree()

# ---------------------
allSkyCmap = matlab_parula_cmap()
# 'tab:red', 'tab:orange'
# LowFlyerColor = 'royalblue'
# HighFlyerColor = 'red'
LowFlyerColor = 'tab:orange'
HighFlyerColor = 'tab:red'
LowFlyerProjectionColor = 'cyan'
HighFlyerProjectionColor = 'darkorange'

# ---------------------
AltvsLatLineWidth = 4
LatvsLongLineWidth = 4
AllSkyLineWidth = 4

# stylization parameters of text labels in alt vs lat and lat vs long plots
textOffset_alt = 0.15
textOffset_latlong = [[-0.15, -0.15], [0.15, 0.25]]
text_alignment = ['left', 'right']
textUTC_style_lat = [dict(size=20, color=HighFlyerColor), dict(size=20, color=LowFlyerColor)]
textUTC_style_project_lat = [dict(size=20, color=HighFlyerProjectionColor), dict(size=20, color=LowFlyerProjectionColor)]
textUTC_style_alt = [dict(size=20, color=HighFlyerColor), dict(size=20, color=LowFlyerColor)]
altrounding = 1  # how many decimals to round to

# --- cbar ---
cbarVmin = 0
cbarVmax = 16000




# --- FUNCTION ---
def AllSkyTrajecMovie(justPrintSiteNames, rocketFolderPath):

    # --- load attributes for ACESII traj data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID
    allSkySiteFolder = f'{rocketFolderPath}\\{inputPath_modifier_AllSky}'
    sites = [path.replace(f'{allSkySiteFolder}\\','') for path in glob(f'{allSkySiteFolder}\\*')]


    if justPrintSiteNames:
        for i, file in enumerate(sites):
            print('[{:.0f}] {:80s}'.format(i, sites[i]))
        return


    #############################
    # --- get the input files ---
    #############################

    # --- get All Sky data ---
    inputFiles_allsky = [glob(f'{ACES_data_folder}\\{inputPath_modifier_AllSky}\\5577\\*.cdf*')[0], glob(f'{ACES_data_folder}\\{inputPath_modifier_AllSky}\\\\6300\\*.cdf*')[0]]
    data_dicts_allSky = [loadDictFromFile(inputFiles_allsky[0]), loadDictFromFile(inputFiles_allsky[1])]

    Epoch_AllSky = [data_dicts_allSky[0]['Epoch'][0], data_dicts_allSky[1]['Epoch'][0]] # get the image time series and the data itself into single variables
    allImages = [data_dicts_allSky[0]['AllSkyImages'][0], data_dicts_allSky[1]['AllSkyImages'][0]]

    # --- traj Data ---
    inputFilesTraj = [glob(f'{rocketFolderPath}\\{inputPath_modifier_attitude}\\{fliers[0]}\\*.cdf')[0], glob(f'{rocketFolderPath}\\{inputPath_modifier_attitude}\\{fliers[1]}\\*.cdf')[0]]

    print('\n')
    print(color.UNDERLINE + f'Creating AllSky Movie' + color.END)

    # --- get the data from the tmCDF file ---
    prgMsg(f'Loading attitude data')
    data_dicts_attitude = []
    for i in range(2):
        data_dict_attitude = loadDictFromFile(inputFilesTraj[i])
        data_dict_attitude['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_attitude['Epoch'][0][i]) for i in (range(len(data_dict_attitude['Epoch'][0])))])
        data_dicts_attitude.append(data_dict_attitude)
    Done(start_time)

    ###################################
    # --- prepare data for plotting ---
    ###################################
    allGlats = [data_dicts_allSky[0]['GLats'][0], data_dicts_allSky[1]['GLats'][0]]
    allGLongs = [data_dicts_allSky[0]['GLongs'][0], data_dicts_allSky[1]['GLongs'][0]]
    EpochRocket = [data_dicts_attitude[0]['Epoch'][0], data_dicts_attitude[1]['Epoch'][0]]
    Alt = [data_dicts_attitude[0]['Alt'][0]/1000, data_dicts_attitude[1]['Alt'][0]/1000]
    Lat = [data_dicts_attitude[0]['Lat'][0], data_dicts_attitude[1]['Lat'][0]]
    Long = [data_dicts_attitude[0]['Long'][0], data_dicts_attitude[1]['Long'][0]]

    # --- determine the timestamps for each photo and rocket data ---
    Epoch_AllSky_tt2000 = [[pycdf.lib.datetime_to_tt2000(time) for time in Epoch_AllSky[0]],
                           [pycdf.lib.datetime_to_tt2000(time) for time in Epoch_AllSky[1]]]

    # for each image timestamp, find the index of the closest Epoch in the High flyer's epoch
    imageIndiciesToRocket = [[], []]
    for i in range(len(Epoch_AllSky_tt2000)):
        for stamp in Epoch_AllSky_tt2000[i]:
            imageIndiciesToRocket[i].append(np.abs(data_dicts_attitude[0]['Epoch'][0] - stamp).argmin())

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
    highFlyerSize, lowFlyerSize = len(Alt[0]), len(Alt[1])

    # --- Append start Values to  ---

    no_of_points_start = np.abs(EpochRocket[0] - EpochRocket[1][0]).argmin()

    newAlt = [Alt[1][0] for i in range(no_of_points_start)]
    newLat = [Lat[1][0] for i in range(no_of_points_start)]
    newLong = [Long[1][0] for i in range(no_of_points_start)]

    for i in range(len(Alt[1])):
        newAlt.append(Alt[1][i])
        newLat.append(Lat[1][i])
        newLong.append(Long[1][i])

    # --- Append the ending values ---
    remainingIndicies = highFlyerSize - (lowFlyerSize + no_of_points_start)

    for i in range(remainingIndicies):
        newAlt.append(Alt[1][-1])
        newLat.append(Lat[1][-1])
        newLong.append(Long[1][-1])

    Alt[1], Lat[1], Long[1] = np.array(newAlt), np.array(newLat), np.array(newLong)


    ############################################
    # --- calculate the ILat/ILong variables ---
    ############################################

    prgMsg('Loading CHAOS Data')

    # re-create the geoILat,geoIlong data
    data_dict_IonoProj = {
        'ILat_RED': [deepcopy(data_dicts_attitude[i]['ILat'][0]) for i in range(2)],
        'ILong_RED': [deepcopy(data_dicts_attitude[i]['ILong'][0]) for i in range(2)],
        'ILat_GREEN': [deepcopy(data_dicts_attitude[i]['ILat'][0]) for i in range(2)],
        'ILong_GREEN': [deepcopy(data_dicts_attitude[i]['ILong'][0]) for i in range(2)],
        'B_slope':[[],[]],
    }

    for rktIdx in range(2): # ROCKET
        B = np.array(CHAOS(Lat[rktIdx], Long[rktIdx], Alt[rktIdx], [pycdf.lib.tt2000_to_datetime(data_dicts_attitude[0]['Epoch'][0][i]) for i in range(len(data_dicts_attitude[0]['Epoch'][0]))]))
        data_dict_IonoProj['B_slope'][rktIdx].append(B[:, 2]/B[:, 1] )

    Done(start_time)

    ##########################################
    # --- ASSIGN IMAGE TIMESTAMP TO ROCKET ---
    ##########################################

    # --- determine the timestamps for each photo and rocket data ---
    # Note: The image timestamps are taken at the beginning of the image's collection period (30seconds),
    # so we will adjust the time tag to be the middle of the integration period by adding 15seconds
    Epoch_AllSky_tt2000 = [[pycdf.lib.datetime_to_tt2000(time) + 15E9 for time in Epoch_AllSky[0]],
                           [pycdf.lib.datetime_to_tt2000(time) + 15E9 for time in Epoch_AllSky[1]]]

    # for each image timestamp, find the index of the closest Epoch_esa in the High flyer's epoch
    imageIndiciesToRocket = [[], []]
    for i in range(len(Epoch_AllSky_tt2000)):
        for stamp in Epoch_AllSky_tt2000[i]:
            imageIndiciesToRocket[i].append(np.abs(data_dicts_attitude[0]['Epoch'][0] - stamp).argmin())

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


    ###########################
    # --- PLOT ALLSKY MOVIE ---
    ###########################
    if createAllSkyMovie:

        ##########################
        # --- INITIALIZE PLOTS ---
        ##########################
        prgMsg('Initializing Allsky Plots')

        # ----------------------
        # --- START PLOTTING ---
        # ----------------------

        # figure size
        plt.style.use('dark_background')
        fig = plt.figure(dpi=dpi)
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)

        # The overall plot shape
        gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1, 4])

        # define the axes
        axAlt = fig.add_subplot(gs[0, :])
        ax5577 = fig.add_subplot(gs[1, 0], projection=projProjection)
        ax6300 = fig.add_subplot(gs[1, 1], projection=projProjection)

        # --- initialize the title ---
        supTitle = fig.suptitle(f'ACESII\n'
                                f'{EpochMovie[0]}\n', fontsize=40,color='white')

        # --- plot the norwegian map data ---
        # 5577A
        # ax5577.set_title('GREEN 557.7 nm (150 km)', fontsize=30)
        gl5577 = ax5577.gridlines(draw_labels=True, linewidth=3, alpha=0.4, linestyle='--')
        gl5577.xlabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
        gl5577.ylabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
        gl5577.right_labels=False
        ax5577.set_extent([lonW, lonE, latS, latN])  # controls lat/long axes display
        ax5577.coastlines(resolution=res, color='white', alpha=0.8)  # adds coastlines with resolution
        ax5577.set_aspect(0.3)


        # 6300A
        # ax6300.set_title('RED 630.0 nm (250 km)', fontsize=30)
        gl6300 = ax6300.gridlines(draw_labels=True, linewidth=3, color='white', alpha=0.4, linestyle='--')
        gl6300.xlabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
        gl6300.ylabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
        ax6300.set_extent([lonW, lonE, latS, latN])  # controls lat/long axes display
        ax6300.coastlines(resolution=res, color='white', alpha=0.8)  # adds coastlines with resolution
        ax6300.set_aspect(0.3)


        # --- --- --- --- --- --- --- --- ----
        # --- initialize the All Sky Image ---
        # --- --- --- --- --- --- --- --- ----

        # 5577A
        allSky5577 = ax5577.pcolormesh(allGLongs[0], allGlats[0], allImages[0][2], cmap=allSkyCmap, transform=projTransform, vmin=cbarVmin, vmax=cbarVmax)

        # 6300A
        allSky6300 = ax6300.pcolormesh(allGLongs[1], allGlats[1], allImages[1][2], cmap=allSkyCmap, transform=projTransform, vmin=cbarVmin, vmax=cbarVmax)

        # --- initialize total trajectory on the All Sky Image ---
        # --- GREEN ---
        # plot on 5577 Allsky
        ax5577.plot(Long[1], Lat[1], color=LowFlyerColor, linewidth=4, transform=projTransform)
        ax5577.plot(Long[0], Lat[0], color=HighFlyerColor, linewidth=4, transform=projTransform) if not useOnlyLowFlyer else [0,]

        # plot marker on 5577 Allsky
        AllSky5577_marker_low = ax5577.scatter(Long[1][0], Lat[1][0], color=LowFlyerColor, marker='x', linewidth=10, transform=projTransform)
        AllSky5577_marker_high = ax5577.scatter(Long[0][0], Lat[0][0], color=HighFlyerColor, marker='x', linewidth=10, transform=projTransform) if not useOnlyLowFlyer else [0, ]

        # --- RED ---
        # plot on 6300 Allsky
        ax6300.plot(Long[1], Lat[1], color=LowFlyerColor,linewidth=4, transform=projTransform)
        ax6300.plot(Long[0], Lat[0], color=HighFlyerColor,linewidth=4, transform=projTransform) if not useOnlyLowFlyer else [0,]

        # plot marker on 6300 Allsky
        AllSky6300_marker_low = ax6300.scatter(Long[1][0], Lat[1][0], color=LowFlyerColor, linewidth=10, marker='x', transform=projTransform)
        AllSky6300_marker_high = ax6300.scatter(Long[0][0], Lat[0][0], color=HighFlyerColor, linewidth=10, marker='x', transform=projTransform) if not useOnlyLowFlyer else [0,]

        # --- --- --- --- --- --- --- --- --- ---
        # --- initialize Altitude vs Lat Plot ---
        # --- --- --- --- --- --- --- --- --- ---
        axAlt.set_ylabel('Altitude [km]', fontsize=20)
        axAlt.set_xlabel('Geo Lat', fontsize=20)
        axAlt.tick_params(axis='x', which='both', labelsize=25)
        axAlt.tick_params(axis='y', which='both', labelsize=25)
        axAlt.axhline(projectionAltitude[0], color='limegreen', linestyle='--', alpha=0.5, linewidth=AltvsLatLineWidth)
        axAlt.axhline(projectionAltitude[1], color='limegreen', linestyle='--', alpha=0.5, linewidth=AltvsLatLineWidth)
        axAlt.text(69.05, 170, '150 km', color='limegreen', fontsize=20)
        axAlt.text(69.05, 270, '250 km', color='limegreen', fontsize=20)

        # trajectory
        axAlt.plot(Lat[0], Alt[0], color=HighFlyerColor, linewidth=AltvsLatLineWidth) if not useOnlyLowFlyer else [0,]
        axAlt.plot(Lat[1], Alt[1], color=LowFlyerColor, linewidth=AltvsLatLineWidth)

        # B-projection
        legendAltLatLabel = f'Projected IGRF'
        BprojectionHigh, = axAlt.axline(xy1=Alt[0][0], slope=data_dict_IonoProj['B_slope'][0][0], color='white', alpha=0.7, linestyle='--')
        BprojectionLow, = axAlt.axline(xy1=Alt[1][0], slope=data_dict_IonoProj['B_slope'][1][0], color='white', alpha=0.7, linestyle='--')
        # BprojectionHigh, = axAlt.plot(projectB[0][0][0], projectB[0][0][1], color='white', alpha=0.7, linestyle='--') if not useOnlyLowFlyer else [0,]  # high flyer
        # BprojectionLow, = axAlt.plot(projectB[1][0][0], projectB[1][0][1], color='white', alpha=0.7, linestyle='--', label=legendAltLatLabel)  # low flyer

        # legend for AltvsLat projection
        axAlt.legend(loc="upper left", fontsize=15)

        # marker
        Altlat_marker_high = axAlt.scatter(Lat[0][0], Alt[0][0], color=HighFlyerColor, marker='x', linewidth=15) if not useOnlyLowFlyer else [0,]
        Altlat_marker_low = axAlt.scatter(Lat[1][0], Alt[1][0], color=LowFlyerColor, marker='x', linewidth=15)

        # text labels
        Altlat_text_high = axAlt.text(Lat[0][0] - textOffset_alt, Alt[0][0],
                                   f'{round(Alt[0][0], altrounding)} km', ha=text_alignment[1],
                                   **textUTC_style_alt[0]) if not useOnlyLowFlyer else [0,]
        Altlat_text_low = axAlt.text(Lat[1][0] - textOffset_alt, Alt[1][0],
                                  f'{round(Alt[1][0], altrounding)} km', ha=text_alignment[1],
                                  **textUTC_style_alt[1])

        # --- intialize Trajectory ILat vs Ilong projection ---

        # plot the lat/long coordinates that update on the 5577 graph
        # IlatILong_text_low_GREEN = ax5577.text(12, 73.25, f'({round(data_dict_IonoProj["geoILong_GREEN"][1][0], altrounding)}' + '$^{\circ}$' + f', {round(data_dict_IonoProj["geoILat_GREEN"][1][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_project_lat[1], transform=projPC)
        # IlatILong_text_high_GREEN = ax5577.text(12.05, 73, f'({round(data_dict_IonoProj["geoILong_GREEN"][0][0], altrounding)}' + '$^{\circ}$' + f', {round(data_dict_IonoProj["geoILat_GREEN"][0][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_project_lat[0],transform=projPC) if not useOnlyLowFlyer else [0, ]
        # latLong_text_low_GREEN = ax5577.text(12.1, 72.75, f'({round(Long[1][0], altrounding)}' + '$^{\circ}$' + f', {round(Lat[1][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_lat[1], transform=projPC)
        # latLong_text_high_GREEN = ax5577.text(12.15, 72.5, f'({round(Long[0][0], altrounding)}' + '$^{\circ}$' + f', {round(Lat[0][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_lat[0],transform=projPC) if not useOnlyLowFlyer else [0, ]

        ax5577.legend(loc='upper left', prop={'size': 14})

        # plot the lat/long coordinates that update on the 6300 graph
        # IlatILong_text_low_RED = ax6300.text(12.5, 70.25, f'({round(data_dict_IonoProj["geoILong_RED"][1][0], altrounding)}' + '$^{\circ}$' + f', {round(data_dict_IonoProj["geoILat_RED"][1][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_project_lat[1], transform=projPC)
        # IlatILong_text_high_RED = ax6300.text(12.55, 70, f'({round(data_dict_IonoProj["geoILong_RED"][0][0], altrounding)}' + '$^{\circ}$' + f', {round(data_dict_IonoProj["geoILat_RED"][0][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_project_lat[0], transform=projPC) if not useOnlyLowFlyer else [0, ]
        # latLong_text_low_RED = ax6300.text(12.6, 69.75, f'({round(Long[1][0], altrounding)}' + '$^{\circ}$' + f', {round(Lat[1][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_lat[1], transform=projPC)
        # latLong_text_high_RED = ax6300.text(12.65, 69.5, f'({round(Long[0][0], altrounding)}' + '$^{\circ}$' + f', {round(Lat[0][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_lat[0], transform=projPC) if not useOnlyLowFlyer else [0, ]

        # ax6300.legend(loc='upper left', prop={'size': 14})

        # --- Add the cbar ---
        cax = fig.add_axes([0.932, 0.0325, 0.0155, 0.5])
        cbar = plt.colorbar(mappable=allSky5577, cax=cax, orientation='vertical', fraction=0.046, pad=0.04)
        cbar.set_label('Brightness (Rayleigh)',fontsize=20)
        # cbar.tick_params(labelsize=20)

        Done(start_time)

        fig.savefig(r'C:\Data\ACESII\Movies_and_Media\movies\garbage.png')



        ###########################
        # --- ANIMATE THE MOVIE ---
        # ###########################
        def animatePlot(i):
            print('', end='\r' + color.RED + f'{round(100 * i / len(EpochMovie), 1)} %' + color.END)

            # update the title
            supTitle.set_text(f'ACESII\n'
                              f'{EpochMovie[i]}\n')

            # --- ALL SKY ---
            # update all sky images
            wImage5577 = imgIndicies[0][i]
            wImage6300 = imgIndicies[1][i]

            cmap5577.set_array(allImages[0][imgIndicies[0][i]].ravel())
            cmap6300.set_array(allImages[1][imgIndicies[1][i]].ravel())

            # update the allsky5577 rocket marker positions
            AllSky5577_marker_high.set_offsets([Long[0][i], Lat[0][i]]) if not useOnlyLowFlyer else [0, ]
            AllSky5577_marker_low.set_offsets([Long[1][i], Lat[1][i]])
            # AllSky5577_projected_marker_high.set_offsets([data_dict_IonoProj['geoILong_GREEN'][0][i], data_dict_IonoProj['geoILat_GREEN'][0][i]]) if not useOnlyLowFlyer else [0, ]
            # AllSky5577_projected_marker_low.set_offsets([data_dict_IonoProj['geoILong_GREEN'][1][i], data_dict_IonoProj['geoILat_GREEN'][1][i]])

            # update the allsky6300 rocket marker positions
            AllSky6300_marker_high.set_offsets([Long[0][i], Lat[0][i]]) if not useOnlyLowFlyer else [0, ]
            AllSky6300_marker_low.set_offsets([Long[1][i], Lat[1][i]])
            # AllSky6300_projected_marker_high.set_offsets([data_dict_IonoProj['geoILong_RED'][0][i], data_dict_IonoProj['geoILat_RED'][0][i]]) if not useOnlyLowFlyer else [0, ]
            # AllSky6300_projected_marker_low.set_offsets([data_dict_IonoProj['geoILong_RED'][1][i], data_dict_IonoProj['geoILat_RED'][1][i]])

            # --- ALT VS LAT ---
            # update Alt vs Lat marker
            Altlat_marker_high.set_offsets([Lat[0][i], Alt[0][i]]) if not useOnlyLowFlyer else [0, ]
            Altlat_marker_low.set_offsets([Lat[1][i], Alt[1][i]])

            # update B-project lines on Alt vs lat plot
            BprojectionHigh.set_xdata(projectB[0][i][0]) if not useOnlyLowFlyer else [0, ]
            BprojectionHigh.set_ydata(projectB[0][i][1]) if not useOnlyLowFlyer else [0, ]
            BprojectionLow.set_xdata(projectB[1][i][0])
            BprojectionLow.set_ydata(projectB[1][i][1])

            # update Alt vs lat text
            Altlat_text_high.set_x(Lat[0][i]) if not useOnlyLowFlyer else [0, ]
            Altlat_text_high.set_y(Alt[0][i]) if not useOnlyLowFlyer else [0, ]
            Altlat_text_high.set_text(f'{round(Alt[0][i], altrounding)} km') if not useOnlyLowFlyer else [0, ]

            Altlat_text_low.set_x(Lat[1][i])
            Altlat_text_low.set_y(Alt[1][i])
            Altlat_text_low.set_text(f'{round(Alt[1][i], altrounding)} km')

            # --- GREEN ---
            # update lat vs long text
            # latLong_text_high_GREEN.set_text(f'({round(Long[0][i], altrounding)}' + '$^{\circ}$' +
            #                            f', {round(Lat[0][i], altrounding)}' + '$^{\circ}$)') if not useOnlyLowFlyer else [0, ]
            # latLong_text_low_GREEN.set_text(f'({round(Long[1][i], altrounding)}' + '$^{\circ}$' +
            #                            f', {round(Lat[1][i], altrounding)}' + '$^{\circ}$)')

            # update Ilat vs Ilong text
            # IlatILong_text_high_GREEN.set_text(f'({round(data_dict_IonoProj["geoILong_GREEN"][0][i], altrounding)}' + '$^{\circ}$' +
            #                            f', {round(data_dict_IonoProj["geoILat_GREEN"][0][i], altrounding)}' + '$^{\circ}$)') if not useOnlyLowFlyer else [0,]
            #
            # IlatILong_text_low_GREEN.set_text(f'({round(data_dict_IonoProj["geoILong_GREEN"][1][i], altrounding)}' + '$^{\circ}$' +
            #                           f', {round(data_dict_IonoProj["geoILat_GREEN"][1][i], altrounding)}' + '$^{\circ}$)')

            # --- RED ---
            # update lat vs long text
            # latLong_text_high_RED.set_text(f'({round(Long[0][i], altrounding)}' + '$^{\circ}$' +
            #                                  f', {round(Lat[0][i], altrounding)}' + '$^{\circ}$)') if not useOnlyLowFlyer else [0, ]
            # latLong_text_low_RED.set_text(f'({round(Long[1][i], altrounding)}' + '$^{\circ}$' +
            #                                 f', {round(Lat[1][i], altrounding)}' + '$^{\circ}$)')

            # update Ilat vs Ilong text
            # IlatILong_text_high_RED.set_text(f'({round(data_dict_IonoProj["geoILong_RED"][0][i], altrounding)}' + '$^{\circ}$' +
            #                                    f', {round(data_dict_IonoProj["geoILat_RED"][0][i], altrounding)}' + '$^{\circ}$)') if not useOnlyLowFlyer else [0, ]
            #
            # IlatILong_text_low_RED.set_text(f'({round(data_dict_IonoProj["geoILong_RED"][1][i], altrounding)}' + '$^{\circ}$' +
            #                                   f', {round(data_dict_IonoProj["geoILat_RED"][0][i], altrounding)}' + '$^{\circ}$)')

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
        anim.save(rf'C:\Data\ACESII\Movies_and_Media\movies\ACESII_AllSky_Movie.mp4',writer=writervideo)


        Done(start_time)








# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder

if justPrintSiteNames:
    AllSkyTrajecMovie(justPrintSiteNames,rocketFolderPath)
else:
    AllSkyTrajecMovie(justPrintSiteNames,rocketFolderPath)
