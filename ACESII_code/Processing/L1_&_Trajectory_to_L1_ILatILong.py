# --- L1_&_Trajectory_to_L1_ILatILong.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Takes in the ACESII trajectory data (in geodetic and geomagnetic) as well
# as the ESA data (individual Instruments). Interpolates the Geomag and Alt data then produces the Ionospheric Projected Lattitude and
# Longitude


# TODO: This file wont work if there's multple ESA datafiles with "eepaa" or "iepaa" etc in the L1 folder! Must only be one.Should fix this

# --- --- --- --- ---
import time
from ACESII_code.class_var_func import Done, setupPYCDF,prgMsg
start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

modifier = '_ILatILong'

targetProjectionAltitude = 100 # altitude you want to project B (in km). Should be ~100km

# select which files to convert
wInstr = 'eepaa' # valid inputs are strings of eepaa, iepaa and leesa


##########################
# --- PLOTTING TOGGLES ---
##########################

# --- One Flyer on plot ---
plot_BFieldVectors = False
plot_B_Projection_single = False

# --- Both Flyers ---
plot_B_Projection_both = False
useKM = True # default is to use LAT/LONG

# sets the reference point of the origin to Andoya Space
useAndoya = True
plotConnectingLines = True

# --- ESA Overlay Plots ---
plot_ESA_Overlay = True
wPitch = 1 # values from 0,1,2 ... 20 ---> -10deg, 0deg, 10deg ...



# --- output data? ---
outputDataFile = False




# --- --- --- ---
# --- import ---
# --- --- --- ---
import numpy as np
import pyIGRF
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from os import remove,path
from glob import glob
from ACESII_code.data_paths import fliers, ACES_data_folder, ACES_csv_trajectories, TRICE_data_folder, ACES_L0_files, TRICE_L0_files
from ACESII_code.class_var_func import color, L1_TRICE_Quick
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
setupPYCDF()
from spacepy import pycdf
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field

print(color.BOLD + color.CYAN + 'L1_&_Trajectory_to_L1_ILatILong.py' + color.END + color.END)
def Trajectory_to_ESA_ILatILong(wInstr, rocketFolderPath):

    if wInstr == 'leesa':
        rangelen = 1
    else:
        rangelen = 2

    # --- ACES II Flight/Integration Data ---
    rocketAttrs,b,c = ACES_mission_dicts()
    globalAttrsMod = rocketAttrs.globalAttributes[0]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L1'
    L2ModelData = L1_TRICE_Quick(0)

    # Set the paths for the file names
    TrajectoryFiles = []
    L1Files = []
    L1ILatILongFiles = []
    L1_names = []
    L1ILatILong_names = []
    dataFile_name = []
    fileoutName = []

    for i in range(rangelen):
        TrajectoryFiles.append(glob(f"{rocketFolderPath}trajectories\{fliers[i]}\*.cdf"))

        L1Files.append(glob(f'{rocketFolderPath}L1\{fliers[i]}\*{wInstr}*'))
        L1ILatILongFiles.append(glob(f'{rocketFolderPath}L1\{fliers[i]}{modifier}\*{wInstr}*'))

        L1_names.append([ifile.replace(f'{rocketFolderPath}L1\{fliers[i]}\\', '') for ifile in L1Files[i]])
        L1ILatILong_names.append([ofile.replace(f'{rocketFolderPath}L1\{fliers[i]}{modifier}\\', '') for ofile in L1ILatILongFiles[i]])

        dataFile_name.append(L1Files[i][0].replace(f'{rocketFolderPath}L1\{fliers[i]}\\', ''))
        fileoutName.append(dataFile_name[i].replace('l1', 'l1_ILat_ILong'))

    if wInstr == 'lp':
        print(color.RED + 'Cannot processes Langmuir Probe File' + color.END)
    else:
        print(color.UNDERLINE + f'Processing to L1 data for {wInstr.upper()} instrument' + color.END)

        ######################
        # --- LOAD IN DATA ---
        ######################
        prgMsg('Loading data from L1 and Trajectory Files')
        data_dicts = []
        data_dicts_traj = []

        for i in range(rangelen):

            data_dict = {}
            with pycdf.CDF(L1Files[i][0]) as L1DataFile:
                for key, val in L1DataFile.items():
                    data_dict = {**data_dict, **{key : [L1DataFile[key][...] , {key:val for key,val in L1DataFile[key].attrs.items()  }  ]  }  }

            data_dict['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_esa'][0][i]) for i in range(len(data_dict['Epoch_esa'][0]))])

            data_dicts.append(data_dict)

            data_dict_traj = {}

            with pycdf.CDF(TrajectoryFiles[i][0]) as TrajDataFile:
                for key,val in TrajDataFile.items():
                    data_dict_traj = {**data_dict_traj, **{key: [TrajDataFile[key][...], {key: val for key, val in TrajDataFile[key].attrs.items()}]}}

            data_dict_traj['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch'][0][i]) for i in range(len(data_dict_traj['Epoch'][0]))])

            data_dicts_traj.append(data_dict_traj)

        Done(start_time)


        ##############################################
        # --- DownSample GPS Alt, Lat, Long, Epoch ---
        ##############################################

        prgMsg('Downsampling GPS data')
        # Using the Trajectory data's Epoch, downsample the traj data to align with the ESA data
        AltDS = [[],[]]
        LatDS = [[],[]]
        LongDS = [[],[]]
        EpochDS = [[],[]]

        for j in range(rangelen):

            for i in range(len(data_dicts[j]['Epoch_esa'][0])):
                targetIndex = np.abs(data_dicts_traj[j]['Epoch'][0] -  data_dicts[j]['Epoch_esa'][0][i]).argmin()
                AltDS[j].append(data_dicts_traj[j]['Alt'][0][targetIndex])
                LatDS[j].append(data_dicts_traj[j]['Lat'][0][targetIndex])
                LongDS[j].append(data_dicts_traj[j]['Long'][0][targetIndex])
                EpochDS[j].append(data_dicts_traj[j]['Epoch'][0][targetIndex])

            AltDS[j] = np.array(AltDS[j])
            LatDS[j] = np.array(LatDS[j])
            LongDS[j] = np.array(LongDS[j])
            EpochDS[j] = np.array(EpochDS[j])

        Done(start_time)

        ################################
        # --- Calculate IGRF B-Field ---
        ################################
        prgMsg('Getting IGRF Field')
        # -- Output order forpyIGRF.igrf_value ---
        # Declination (+ E | - W)
        # Inclination (+ D | - U), should be ~78deg for what we're doing
        # Horizontal Intensity
        # North Comp (+ N | - S)
        # East Comp (+ E | - W)
        # Vertical Comp (+ D | - U)
        # Total Field

        IGRF = [[],[]]
        date = 2022 + 323 / 365  # Corresponds to 11/20/2022
        for i in range(rangelen):
            for j in range(len(data_dicts[i]['Epoch_esa'][0])):
                IGRF[i].append(pyIGRF.igrf_value(LatDS[i][j], LongDS[i][j], AltDS[i][j], date))

        Done(start_time)

        #################################
        # --- I-LAT I-LONG PROJECTION ---
        #################################

        prgMsg('Projecting B-Fields')

        lat_to_meter = 111.319488 # 1 deg latitude to kilometers on Earth
        def long_to_meter(lat):
            return 111.319488 * math.cos(lat*math.pi/180)

        intersectionPoints = [[],[]]
        BFieldDirNorm = [[],[]]
        BFieldLoc = [[],[]]


        for j in range(rangelen):

            for i in range(len(data_dicts[j]['Epoch_esa'][0])):
                #Coordiantes reported in (Long (x) , Lat (y), Alt (z))
                vLoc = np.array([LongDS[j][i],LatDS[j][i],AltDS[j][i]]) # IGRF vector Location, should be rocket coordinates
                vDir = np.array([IGRF[j][i][4],IGRF[j][i][3],-1*IGRF[j][i][5]]) # IGRF vector Direction. -1 added in third elemental due to down being positive in IGRF given
                vDirNorm = vDir / np.linalg.norm(vDir) # Normalize IGRF to get its direction only. This will make t larger, but that's fine
                BFieldDirNorm[j].append(vDirNorm)
                BFieldLoc[j].append(vLoc)

                # Determine the Delta-Latitude and Longitutde
                Theta_dec = IGRF[j][i][0]
                Theta_in = IGRF[j][i][1]
                h = vLoc[2] - targetProjectionAltitude
                deltaLat = (1/lat_to_meter) * h*math.tan((math.pi/180) *(90 - Theta_in))
                deltaLong = (1/long_to_meter(LatDS[j][i])) * h*math.tan((math.pi/180) * (Theta_dec))
                intersectionPoints[j].append([vLoc[0]+deltaLong, vLoc[1] + deltaLat, 0])

        Done(start_time)

        ############################################
        # --- CONVERT TO GEOMAGNETIC COORDINATES ---
        ############################################
        prgMsg('Converting Coordinates')
        geodeticLongIntersects = [[],[]]
        geodeticLatIntersects = [[],[]]
        geodeticAltIntersects = [[],[]]
        geodetic = [[],[]]
        geoMagFootPrint_lat = [[],[]]
        geoMagFootPrint_long = [[],[]]

        for j in range(rangelen):
            geodeticLongIntersects[j] = np.array([intersectionPoints[j][i][0] for i in range(len(intersectionPoints[j]))]) # Long
            geodeticLatIntersects[j] = np.array([intersectionPoints[j][i][1] for i in range(len(intersectionPoints[j]))])  # Lat
            geodeticAltIntersects[j] = np.array([intersectionPoints[j][i][2] for i in range(len(intersectionPoints[j]))])  # Alt
            geodetic[j] = np.array([geodeticAltIntersects[j], geodeticLatIntersects[j], geodeticLongIntersects[j]]).transpose()

            ISOtime = [ pycdf.lib.tt2000_to_datetime(EpochDS[j][i]).isoformat() for i in range(len(EpochDS[j]))]
            cvals_GDZ = coord.Coords(geodetic[j], 'GDZ', 'sph')
            cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
            cvals_GDZ_MAG = cvals_GDZ.convert('MAG', 'sph')

            geoMagFootPrint_lat[j] = cvals_GDZ_MAG.lati
            geoMagFootPrint_long[j] = cvals_GDZ_MAG.long

        Done(start_time)


        #####################
        # --- APPEND DATA ---
        #####################
        for i in range(rangelen):

            dataNEW = {
                'mapped_Iono_lattitude_geomagnetic': np.array(geoMagFootPrint_lat[i]),
                'mapped_Iono_longitude_geomagnetic': np.array(geoMagFootPrint_long[i]),
                'mapped_Iono_longitude_geodetic': np.array(geodeticLongIntersects[i]),
                'mapped_Iono_lattitude_geodetic': np.array(geodeticLatIntersects[i])}

            # Add the intersection latitude and longitude to the data_dict
            for key, val in dataNEW.items():
                data_dicts[i] = {**data_dicts[i], **{key: [val, {'LABLAXIS': key,
                                                                 'DEPEND_0': 'Epoch_esa', 'DEPEND_1': None,
                                                                 'DEPEND_2': None,
                                                                 'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                                 'UNITS': 'deg',
                                                                 'VALIDMIN': val.min(), 'VALIDMAX': val.max(),
                                                                 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}


        ##################
        # --- PLOTTING ---
        ##################

        if plot_BFieldVectors:

            for j in range(rangelen):
                start = 0
                end = len(BFieldDirNorm[j])
                skip = 350

                V = 10*np.array(BFieldDirNorm[j][start:end:skip])
                origin = np.array(BFieldLoc[j][start:end:skip])  # origin point

                xDirData = np.array([V[i][0] for i in range(len(V)) ])
                yDirData = np.array([V[i][1] for i in range(len(V)) ])
                zDirData = np.array([V[i][2] for i in range(len(V)) ])

                xLocData = np.array([origin[i][0] for i in range(len(origin))])
                yLocData = np.array([origin[i][1] for i in range(len(origin))])
                zLocData = np.array([origin[i][2] for i in range(len(origin))])

                fig = plt.figure(figsize=(10, 10))
                ax1 = fig.add_subplot(121)
                ax1.quiver(xLocData, zLocData, xDirData, zDirData, scale=41)
                ax1.set_xlabel('Long')
                ax1.set_ylabel('Alt')

                ax2 = fig.add_subplot(122)
                ax2.quiver(yLocData, zLocData, yDirData, zDirData, scale=41)
                ax2.set_xlabel('Lat')
                ax2.set_ylabel('Alt')

                plt.title(rf'{fliers[j]} flyer B-Field Direction')
                fig.savefig(rf'D:\Data\ACESII\trajectories\trajectory_plots\{fliers[j]}\{fliers[j]}_BVectors.png')

        if plot_B_Projection_single:

            for i in range(rangelen):

                # ticks
                spacing_of_majorticks = 1
                tick_params_major1 = dict(labelsize=14, which='major', size=10, pad=3)
                tick_params_major2 = dict(labelsize=14, which='major', size=10, pad=3)
                tick_params_minor1 = dict(labelsize=10, which='minor', size=5)
                tick_params_minor2 = dict(labelsize=10, which='minor', size=5)
                labelsize = 12

                # geographic long vs  geodetic/geomag lat
                fig = plt.figure(figsize=(9,9))
                ax1 = fig.add_subplot(111)
                ax1.minorticks_on()
                ax1.xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
                ax1.xaxis.set_major_formatter('{x:.0f}')
                ax1.xaxis.set_minor_locator(AutoMinorLocator())
                ax1.plot(geodeticLatIntersects[i], geodeticLongIntersects[i], color='blue',label = 'B-Field Projection - Geographic')
                ax1.plot(LatDS[i],LongDS[i],color='red',label='Rocket Ground Track')
                ax1.set_xlabel('Geographic Lat',size=labelsize)
                ax1.set_ylabel('Geographic Long',size=labelsize)
                ax1.tick_params(**tick_params_major1)
                ax1.tick_params(**tick_params_minor1)
                ax1.legend()

                # Get the specific time stamps and plot some vlines
                step = 2000
                epochStepped = EpochDS[i][step:len(EpochDS[j]):step]
                latStepped = LatDS[i][step:len(LatDS[j]):step]
                longStepped =LongDS[i][step:len(LongDS[j]):step]
                epochStepped_dt = [ pycdf.lib.tt2000_to_datetime(epochStepped[j]).time() for j in range(len(epochStepped))]
                epochStepped_text = [time.strftime("%H:%M:%S") for time in epochStepped_dt]


                for k,text in enumerate(epochStepped_text):
                    ax1.text(latStepped[k]-0.1, longStepped[k]-0.1, text,ha='right')
                    ax1.axvline(x = latStepped[k], color = 'black',alpha=0.1,linestyle='--')

                ax1.scatter(latStepped, longStepped, marker = 'x',color = 'red',s = 100)
                ax2 = ax1.twiny()
                ax2.minorticks_on()
                ax2.plot(geoMagFootPrint_lat[i], geodeticLongIntersects[i], color='blue',alpha = 0.0)
                ax2.set_xlabel('Geomagnetic Lat',size=labelsize)
                ax2.minorticks_on()
                ax2.xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
                ax2.xaxis.set_major_formatter('{x:.0f}')
                ax2.xaxis.set_minor_locator(AutoMinorLocator())
                ax2.tick_params(**tick_params_major2)
                ax2.tick_params(**tick_params_minor2)

                rflyers = ['HighFlyer_','LowFlyer_']

                plt.title(f'{wInstr.upper()} B-Field Projection Altitude: {targetProjectionAltitude} km ')

                fig.savefig(rf'D:\Data\ACESII\trajectories\trajectory_plots\{fliers[i]}\{rflyers[i]}BProjection.png')

        if plot_B_Projection_both:

            ######################
            # --- PLOT TOGGLES ---
            ######################

            # figure
            figureSize = (9,9)

            # ticks
            tick_params_major1 = dict(labelsize=14, which='major', size=10, pad=3) # geodetic
            tick_params_major2 = dict(labelsize=14, which='major', size=10, pad=3) # geomag
            tick_params_minor1 = dict(labelsize=10, which='minor', size=5) # geodetic
            tick_params_minor2 = dict(labelsize=10, which='minor', size=5) # geomag
            labelsize = 12

            step = 1250 # index spacing between "X" timestamp marks
            step_adjust = 200 # intial adjustment to the stepsize above
            colors = ['red', 'blue'] # high/low flyer groundtrack colors
            ha = ['left', 'right'] # high/low flyer position of timestamp labels for
            projectcolors = ['salmon', 'deepskyblue']

            # Andoya Long/Lat
            if useAndoya:
                refLat = rocketAttrs.Andoya_Space_Lat_Long[0]
                refLong = rocketAttrs.Andoya_Space_Lat_Long[1]
            else:
                refLat = 0
                refLong = 0

            if useKM:
                adjusts = [[1, 0], [-1 * 2, -1 * 8]] # km adjustments to the timestamp labels
                mod_file_ID = 'km'
                spacing_of_majorticks = 50
            else:
                adjusts = [[0.05, 0.05], [-0.05, -0.05]] # long/lat adjustments to the timestamp labels
                mod_file_ID = 'latlong'
                spacing_of_majorticks = 1

            # --- --- --- --- --- --- --- --- ---
            # --- COLLECT DATA TO BE PLOTTED ---
            # --- --- --- --- --- --- --- --- ---

            # Get the specific time stamps and plot some vlines
            latStepped = [[], []]
            longStepped = [[], []]
            latStepped_km = [[], []]
            longStepped_km = [[], []]
            geodeticLatIntersects_km = [[], []]
            geodeticLongIntersects_km = [[], []]
            LatDS_km = [[], []]
            LongDS_km = [[], []]
            latIntersectsStepped = [[], []]
            longIntersectsStepped = [[], []]
            latIntersectsStepped_km = [[], []]
            longIntersectsStepped_km = [[], []]
            epochStepped_text = [[],[]]

            for i in range(rangelen):

                # Collect STEPPED data
                epochStepped = EpochDS[i][step + step_adjust + int((step / 2)):len(EpochDS[i]):step]
                epochStepped_dt = [pycdf.lib.tt2000_to_datetime(epochStepped[j]).time() for j in range(len(epochStepped))]
                epochStepped_text[i] = [time.strftime("%H:%M:%S") for time in epochStepped_dt]
                latStepped[i] = LatDS[i][step + step_adjust + int((step / 2)):len(LatDS[i]):step]
                longStepped[i] = LongDS[i][step + step_adjust + int((step / 2)):len(LongDS[i]):step]

                # Convert lat/long data to KM
                for j in range(len(geodeticLatIntersects[i])):
                    geodeticLatIntersects_km[i].append(lat_to_meter * (geodeticLatIntersects[i][j] - refLat))
                    geodeticLongIntersects_km[i].append(long_to_meter(refLat) * (geodeticLongIntersects[i][j] - refLong))
                    LatDS_km[i].append(lat_to_meter * (LatDS[i][j] - refLat))
                    LongDS_km[i].append(long_to_meter(refLat) * LongDS[i][j] - long_to_meter(refLat) * refLong)

                # collect STEPPED data (again)
                latIntersectsStepped[i] = geodeticLatIntersects[i][step + step_adjust + int((step / 2)):len(geodeticLatIntersects[i]):step]
                longIntersectsStepped[i] = geodeticLongIntersects[i][step + step_adjust + int((step / 2)):len(geodeticLongIntersects[i]):step]
                latIntersectsStepped_km[i] = geodeticLatIntersects_km[i][step + step_adjust + int((step / 2)):len(geodeticLatIntersects_km[i]):step]
                longIntersectsStepped_km[i] = geodeticLongIntersects_km[i][step + step_adjust + int((step / 2)):len(geodeticLongIntersects_km[i]):step]
                latStepped_km[i] = [lat_to_meter*(latStepped[i][h] - refLat) for h in range(len(latStepped[i]))]
                longStepped_km[i] = [long_to_meter(refLat)*(longStepped[i][h] - refLong) for h in range(len(longStepped[i]))]

            # --- --- --- ---
            # --- PLOTTING ---
            # --- --- --- ---
            fig = plt.figure(figsize=figureSize)
            ax1 = fig.add_subplot(111)
            ax1.minorticks_on()
            ax1.xaxis.set_major_formatter('{x:.0f}')

            ax1.tick_params(**tick_params_major1)
            ax1.tick_params(**tick_params_minor1)

            plt.grid(True, which='Major', alpha=0.5)
            plt.grid(True, which='Minor', alpha=0.2)

            # Plot the Trajectory Data
            if useKM:
                # AXIS LABELS
                ax1.set_xlabel(f'(+E|-W) [km]', size=labelsize)
                ax1.set_ylabel(f'(+N|-S) [km]', size=labelsize)

                # TICKS
                major_xticks = np.arange(-50, 50, spacing_of_majorticks)
                major_yticks = np.arange(0, 550, spacing_of_majorticks)
                ax1.set_xticks(major_xticks)
                ax1.set_yticks(major_yticks)
                ax1.xaxis.set_minor_locator(AutoMinorLocator())

                # GROUND TRACK TRAJECTS
                ax1.plot(geodeticLongIntersects_km[0], geodeticLatIntersects_km[0], color=projectcolors[0],label='High Flyer - Geograhpic Projection')
                ax1.plot(geodeticLongIntersects_km[1], geodeticLatIntersects_km[1], color=projectcolors[1],label='Low Flyer - Geographic Projection')
                ax1.plot(LongDS_km[0], LatDS_km[0], color='red', label='High Flyer - Ground Track')
                ax1.plot(LongDS_km[1], LatDS_km[1], color='blue', label='Low Flyer - Ground Track')
            else:
                # AXIS LABELS
                ax1.set_xlabel('Geographic Long', size=labelsize)
                ax1.set_ylabel('Geographic Lat', size=labelsize)

                # TICKS
                ax1.xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
                ax1.xaxis.set_minor_locator(AutoMinorLocator())

                # GROUND TRACK TRAJECTS
                ax1.plot(geodeticLongIntersects[0], geodeticLatIntersects[0], color='salmon',label='High Flyer - Geograhpic Projection')
                ax1.plot(geodeticLongIntersects[1], geodeticLatIntersects[1], color='deepskyblue', label='Low Flyer - Geographic Projection')
                ax1.plot(LongDS[0], LatDS[0], color='red', label='High Flyer - Ground Track')
                ax1.plot(LongDS[1], LatDS[1], color='blue', label='Low Flyer - Ground Track')

                # GEOMAG AXIS
                ax2 = ax1.twinx()
                ax2.minorticks_on()
                ax2.plot(geodeticLongIntersects[0], geoMagFootPrint_lat[0], color='blue', alpha=0.0)
                ax2.plot(geodeticLongIntersects[1], geoMagFootPrint_lat[1], color='red', alpha=0.0)
                ax2.set_ylabel('Geomagnetic Lat', size=labelsize)
                ax2.minorticks_on()
                ax2.xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
                ax2.xaxis.set_major_formatter('{x:.0f}')
                ax2.xaxis.set_minor_locator(AutoMinorLocator())
                ax2.tick_params(**tick_params_major2)
                ax2.tick_params(**tick_params_minor2)


            if useKM:
                # Plot X marks
                ax1.scatter(longStepped_km[0], latStepped_km[0], marker='x', color=colors[0], s=100)
                ax1.scatter(longStepped_km[1], latStepped_km[1], marker='x', color=colors[1], s=100)

                for i in range(rangelen):
                    for k in range(len(epochStepped_text[i])):
                        ax1.text(longStepped_km[i][k] + adjusts[i][0], latStepped_km[i][k] + adjusts[i][1], epochStepped_text[i][k], ha=ha[i], color=colors[i])

                if plotConnectingLines:

                    for i in range(len(latIntersectsStepped[1])):
                        # plot the time tags for the High Flyer Only
                        ax1.text(longIntersectsStepped_km[0][i] + adjusts[0][0], latIntersectsStepped_km[0][i] + adjusts[0][1], epochStepped_text[0][i], ha=ha[0], color=projectcolors[0])
                        ax1.scatter(longIntersectsStepped_km[0][i], latIntersectsStepped_km[0][i], marker='x', color=projectcolors[0], s=100)

                        # Plot the connecting lines
                        x_values = [longIntersectsStepped_km[0][i], longIntersectsStepped_km[1][i]]
                        y_values = [latIntersectsStepped_km[0][i], latIntersectsStepped_km[1][i]]
                        ax1.plot(x_values, y_values, color='green', linestyle="--", alpha=0.7)

            else:
                # Plot X marks
                ax1.scatter(longStepped[0], latStepped[0], marker='x', color=colors[0], s=100)
                ax1.scatter(longStepped[1], latStepped[1], marker='x', color=colors[1], s=100)

                for i in range(rangelen):
                    for k in range(len(epochStepped_text[i])):
                        ax1.text(longStepped[i][k] + adjusts[i][0], latStepped[i][k] + adjusts[i][1], epochStepped_text[i][k], ha=ha[i], color=colors[i])

                if plotConnectingLines:

                    for i in range(len(latIntersectsStepped[1])):
                        # plot the time tags for the High Flyer Only
                        ax1.text(longIntersectsStepped[0][i] + adjusts[0][0], latIntersectsStepped[0][i] + adjusts[0][1], epochStepped_text[0][i], ha=ha[0], color=projectcolors[0])
                        ax1.scatter(longIntersectsStepped[0][i], latIntersectsStepped[0][i], marker='x', color=projectcolors[0], s=100)

                        # Plot the connecting lines
                        x_values = [longIntersectsStepped[0][i], longIntersectsStepped[1][i]]
                        y_values = [latIntersectsStepped[0][i], latIntersectsStepped[1][i]]
                        ax1.plot(x_values, y_values, color='green', linestyle="--", alpha=0.7)


            ax1.legend()
            plt.title(f'{wInstr} \n B-Field Projection Altitude: {targetProjectionAltitude} km \n Origin: {refLat}$^\circ$N {refLong}$^\circ$E')
            fig.savefig(rf'D:\Data\ACESII\trajectories\trajectory_plots\BProjection_{mod_file_ID}.png')

        if plot_ESA_Overlay:

            prgMsg('Plotting ESA overlay')

            ######################
            # --- GET ESA DATA ---
            ######################

            ESAData =[[],[]]
            EpochData = [[],[]]

            for i in range(rangelen):
                esaData_temp = np.array(data_dicts[i][wInstr][0])
                ESAData[i] = esaData_temp[:,wPitch,:]
                EpochData[i] = np.array(data_dicts[i]['Epoch_esa'][0])

            energies = data_dicts[0]['Energy'][0]
            pitches = data_dicts[0]['Pitch_Angle'][0]

            print(len(EpochData[0]), len(EpochData[1]))

            #############################
            # --- "Zoom" into Dataset ---
            #############################

            start = 1750
            end = 1500

            lowFlyer_adjust = 1000

            LatDS = [ LatDS[0][start:-end]  ,LatDS[1][(start+lowFlyer_adjust):(-end)]   ]
            AltDS = [ AltDS[0][start:-end]  ,AltDS[1][(start+lowFlyer_adjust):(-end)]   ]
            ESAData =[ ESAData[0][start:-end]  , ESAData[1][(start+lowFlyer_adjust):(-end)]   ]
            EpochData = [EpochData[0][start:-end], EpochData[1][(start+lowFlyer_adjust):(-end)]   ]

            ESAData = np.array(ESAData, dtype='object')
            EpochData = np.array(EpochData, dtype='object')

            print(len(EpochData[0]),len(EpochData[1]))

            ##################
            # --- PLOTTING ---
            ##################

            # ticks
            tick_params_major1 = dict(labelsize=14, which='major', size=10, pad=3)  # geodetic
            tick_params_major2 = dict(labelsize=14, which='major', size=10, pad=3)  # geomag
            tick_params_minor1 = dict(labelsize=10, which='minor', size=5)  # geodetic
            tick_params_minor2 = dict(labelsize=10, which='minor', size=5)  # geomag
            labelsize = 12

            # Generate the plot
            cdict = {'red': ((0.0, 0.0, 0.0),
                             (0.1, 0.5, 0.5),
                             (0.2, 0.0, 0.0),
                             (0.4, 0.2, 0.2),
                             (0.6, 0.0, 0.0),
                             (0.8, 1.0, 1.0),
                             (1.0, 1.0, 1.0)),
                     'green': ((0.0, 0.0, 0.0),
                               (0.1, 0.0, 0.0),
                               (0.2, 0.0, 0.0),
                               (0.4, 1.0, 1.0),
                               (0.6, 1.0, 1.0),
                               (0.8, 1.0, 1.0),
                               (1.0, 0.0, 0.0)),
                     'blue': ((0.0, 0.0, 0.0),
                              (0.1, 0.5, 0.5),
                              (0.2, 1.0, 1.0),
                              (0.4, 1.0, 1.0),
                              (0.6, 0.0, 0.0),
                              (0.8, 0.0, 0.0),
                              (1.0, 0.0, 0.0))}
            from matplotlib import colors
            my_cmap = colors.LinearSegmentedColormap('my_colormap', cdict, 256)

            x = LatDS[0]
            y = energies
            Z = ESAData[0].transpose()

            fig, ax = plt.subplots(2)
            fig.set_figheight(10)
            fig.set_figwidth(15)

            ####################
            # --- HIGH FLYER ---
            ####################

            # PCOLORMESH
            cmap = ax[0].pcolormesh(x, y, Z,vmin=2,vmax=15,cmap=my_cmap,alpha=0.8)
            cbar = plt.colorbar(cmap,pad=0.1)
            cbar.ax.tick_params(labelsize=20)
            cbar.set_label('Counts',size=20)
            ax[0].set_ylabel('Energy [eV]', size=18)
            ax[0].set_xlabel('Geodetic Lat [deg]', size=18)
            ax[0].set_ylim(bottom = 20, top = 1300)

            # TICKS
            spacing_of_majorticks = 1
            ax[0].xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
            ax[0].xaxis.set_minor_locator(AutoMinorLocator())
            ax[0].minorticks_on()
            ax[0].xaxis.set_major_formatter('{x:.0f}')
            ax[0].tick_params(**tick_params_major1)
            ax[0].tick_params(**tick_params_minor1)

            # Alt AXIS
            spacing_of_majorticks = 1
            ax1 = ax[0].twinx()
            ax1.minorticks_on()
            ax1.plot(LatDS[0], AltDS[0], color='red')
            ax1.set_ylabel('Alt [km]', size=15)
            ax1.minorticks_on()
            ax1.xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
            ax1.xaxis.set_major_formatter('{x:.0f}')
            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(**tick_params_major2)
            ax1.tick_params(**tick_params_minor2)

            # 2ND PLOT
            x = LatDS[1]
            y = energies
            Z = ESAData[1].transpose()

            # PCOLORMESH
            cmap = ax[1].pcolormesh(x, y, Z, vmin=2, vmax=15, cmap=my_cmap, alpha=0.8)
            cbar = plt.colorbar(cmap, pad=0.1)
            cbar.ax.tick_params(labelsize=20)
            cbar.set_label('Counts', size=20)
            ax[1].set_ylabel('Energy [eV]', size=18)
            ax[1].set_xlabel('Geodetic Lat [deg]', size=18)
            ax[1].set_ylim(bottom=20, top=1300)

            # TICKS
            spacing_of_majorticks = 1
            ax[1].xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
            ax[1].xaxis.set_minor_locator(AutoMinorLocator())
            ax[1].minorticks_on()
            ax[1].xaxis.set_major_formatter('{x:.0f}')
            ax[1].tick_params(**tick_params_major1)
            ax[1].tick_params(**tick_params_minor1)

            # Alt AXIS
            spacing_of_majorticks = 1
            ax2 = ax[1].twinx()
            ax2.minorticks_on()
            ax2.plot(LatDS[1], AltDS[1], color='blue')
            ax2.set_ylabel('Alt [km]', size=15)
            ax2.minorticks_on()
            ax2.xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
            ax2.xaxis.set_major_formatter('{x:.0f}')
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(**tick_params_major2)
            ax2.tick_params(**tick_params_minor2)

            ax[0].set_title(f'{wInstr} \n Pitch {pitches[wPitch]}$^\circ$', size=30)
            plt.tight_layout()

            fig.savefig(rf'D:\Data\ACESII\trajectories\trajectory_plots\EEPAA_overlay.png')

            Done(start_time)






        if outputDataFile:
            prgMsg('Writing out Data')

            for i in range(rangelen):

                #####################
                # --- Output Data ---
                #####################

                outputPath = f'{rocketFolderPath}L1\{fliers[i]}{modifier}\\{fileoutName[i]}'

                # --- delete output file if it already exists ---
                if path.exists(outputPath):
                    remove(outputPath)
                pycdf.lib.set_backward(False)

                # --- open the output file ---
                with pycdf.CDF(outputPath, '') as cdfFile:
                    cdfFile.readonly(False)

                    # --- write out global attributes ---
                    globalAttrsMod = rocketAttrs.globalAttributes[i]
                    ModelData = L1_TRICE_Quick(i)
                    inputGlobDic = ModelData.cdfFile.globalattsget()

                    for key, val in inputGlobDic.items():
                        if key == 'Descriptor':
                            globalAttrsMod[key] = 'None'
                        if key in globalAttrsMod:
                            cdfFile.attrs[key] = globalAttrsMod[key]
                        else:
                            cdfFile.attrs[key] = val

                    # --- WRITE OUT DATA ---
                    for varKey, varVal in data_dicts[i].items():
                        if 'Epoch' in varKey:  # epoch data
                            cdfFile.new(varKey, data=varVal[0], type=33)
                        else:  # other data
                            cdfFile.new(varKey, data=varVal[0])

                        # --- Write out the attributes and variable info ---
                        for attrKey, attrVal in data_dicts[i][varKey][1].items():
                            if attrKey == 'VALIDMIN':
                                cdfFile[varKey].attrs[attrKey] = varVal[0].min()
                            elif attrKey == 'VALIDMAX':
                                cdfFile[varKey].attrs[attrKey] = varVal[0].max()
                            elif attrVal != None:
                                cdfFile[varKey].attrs[attrKey] = attrVal

            Done(start_time)






# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder

if len(glob(f'{rocketFolderPath}L1\{fliers[0]}\*.cdf')) == 0 :
    print(color.RED + 'There are no .cdf files in the specified directory:' + color.END)
    print(f'{rocketFolderPath}L1\{fliers[0]}\*.cdf')
elif len(glob(f'{rocketFolderPath}L1\{fliers[1]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory:' + color.END)
    print(f'{rocketFolderPath}L1\{fliers[1]}\*.cdf')
else:
    Trajectory_to_ESA_ILatILong(wInstr, rocketFolderPath)

