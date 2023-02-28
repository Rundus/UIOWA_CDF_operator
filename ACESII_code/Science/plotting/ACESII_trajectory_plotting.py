# --- csv_to_cdf_attitude.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Turn the .cdf files of the TRICE attitude data into cdf files
import datetime
# --- --- --- --- ---
import time as stopwatch
from ACESII_code.class_var_func import Done, setupPYCDF,prgMsg
start_time = stopwatch.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Alt vs Lat/geoLat
plotALTvsLAT = False

# Lat vs Long
plotLATvsLONG = False

# B-Field Projection
plotILATvsILONG = False
useKM = False
wInstr_projection = 'eepaa'

# ESA overlay
plotESAoverlay = False
wRocket_overlay = 5
wPitches_overlay = [1,10,19] # 0deg, 90deg, 180deg
wInstr_overlay = 'eepaa'

# ESA Movie
plotESAmovie = True
plotDistFunc = True # If false, it plots esa L1 counts
polar = False
create_movieESAplots = True # =True if you need to create new plots for the movie, else it uses the plots in plots_for_movies
frame_skips = 50 # 1 is no skip i.e. all frames used, 5 is use only every 5th frame, 10 etc...
wInstr_movie = 'eepaa'

# Use connecting lines in all the relevant plots?
connectingLines = True

# --- --- --- ---
# --- import ---
# --- --- --- ---
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from tqdm import tqdm
from glob import glob
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from ACESII_code.data_paths import fliers, ACES_data_folder
from ACESII_code.class_var_func import color
from ACESII_code.missionAttributes import ACES_mission_dicts
setupPYCDF()
from spacepy import pycdf

print(color.BOLD + color.CYAN + 'ACESII_trajectory_plotting.py' + color.END + color.END)


def storeData(inputDict, dataPaths):
    parentDict = deepcopy(inputDict)
    for pkey, pval in parentDict.items():  # for each instrument, 3 loops

        size = 1 if pkey == 'leesa' else 2  # adjust for leesa only being on High Flyer

        for i in range(size):  # one for each flyer
            data_dict = {}

            try: # Case for if the datafile doesn't exist yet
                with pycdf.CDF(dataPaths[pkey][i][0]) as dataFile:
                    for key, val in dataFile.items():
                        data_dict = {**data_dict, **{key: [dataFile[key][...], {key: val for key, val in dataFile[key].attrs.items()}]}}

                data_dict['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_esa'][0][i]) for i in range(len(data_dict['Epoch_esa'][0]))])
                parentDict[pkey].append(data_dict)
            except:
                parentDict[pkey].append(data_dict)

    return parentDict

def ACESIIplotting():

    rocketAttrs,b,c = ACES_mission_dicts()

    # Trajectory Data
    trajFolderPath = f'{ACES_data_folder}trajectories\\'
    dataPath_traj = [glob(trajFolderPath + rf'{fliers[0]}\\*_ILat_ILong*'),
                     glob(trajFolderPath + rf'{fliers[1]}\\\\*_ILat_ILong*')]

    # L1 ESA data
    FolderPath = f'{ACES_data_folder}L1\\'
    dataPath_L1ESA = {'eepaa':[glob(FolderPath + rf'{fliers[0]}\\*eepaa_*'), glob(FolderPath + rf'{fliers[1]}\\\\*eepaa_*')],
                      'iepaa':[glob(FolderPath + rf'{fliers[0]}\\*iepaa*'),glob(FolderPath + rf'{fliers[1]}\\\\*iepaa*')],
                      'leesa':[glob(FolderPath + rf'{fliers[0]}\\*leesa*'),glob(FolderPath + rf'{fliers[1]}\\\\*leesa*')]}

    # Distribution Function Data
    FolderPath = f'{ACES_data_folder}\science\DistFunc\\'
    dataPath_Dist = {'eepaa':[glob(FolderPath + rf'{fliers[0]}\\*eepaa*'), glob(FolderPath + rf'{fliers[1]}\\\\*eepaa*')],
                      'iepaa':[glob(FolderPath + rf'{fliers[0]}\\*iepaa*'),glob(FolderPath + rf'{fliers[1]}\\\\*iepaa*')],
                      'leesa':[glob(FolderPath + rf'{fliers[0]}\\*leesa*'),glob(FolderPath + rf'{fliers[1]}\\\\*leesa*')]}

    prgMsg('Collecting Trajectory data')

    # --- Create TRAJECTORY Data Dicts ---
    data_dicts_traj = []

    for i in range(len(dataPath_traj)):
        data_dict_traj = {}
        with pycdf.CDF(dataPath_traj[i][0]) as trajFile:
            for key, val in trajFile.items():
                data_dict_traj = {**data_dict_traj, **{key: [trajFile[key][...], {key: val for key, val in trajFile[key].attrs.items()}]}}

        data_dict_traj['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch'][0][i]) for i in range(len(data_dict_traj['Epoch'][0]))])
        data_dicts_traj.append(data_dict_traj)

    Done(start_time)

    # --- Create ESA Data Dicts ---
    prgMsg('Collecting ESA L1 data')
    data_dicts_template = {'eepaa':[],'iepaa':[],'leesa':[]}
    data_dicts_ESA = storeData(data_dicts_template,dataPath_L1ESA)
    Done(start_time)

    # --- Create Distribution_Func Data Dicts ---
    prgMsg('Collecting Distribution Function data')
    data_dicts_template = {'eepaa': [], 'iepaa': [], 'leesa': []}
    data_dicts_dist = storeData(data_dicts_template,dataPath_Dist)
    Done(start_time)


    ####################################################################
    # --- find the points 100s, 200s ... 600s in the Trajectory data ---
    ####################################################################
    prgMsg('Finding Target Times')
    timeTargets = [[100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600], [100, 150, 200, 250, 300, 350, 400]]

    timeTargetsIndices = [[],[]]
    timeTargetsEpoch = []
    for i in range(len(timeTargets)):
        for timeTarg in timeTargets[i]:
            timeTargetsIndices[i].append(np.abs(data_dicts_traj[i]['Epoch'][0] - (rocketAttrs.Launch_Times[i] + timeTarg*(10**(9)))  ).argmin())

        target_in_Epoch = [pycdf.lib.tt2000_to_datetime(data_dicts_traj[i]['Epoch'][0][index]) for index in timeTargetsIndices[i]]
        timeTargetsEpoch.append([targetTime.strftime("%H:%M:%S") for targetTime in target_in_Epoch])

    timeTargetsData = {
        'Epoch' : [[pycdf.lib.tt2000_to_datetime(data_dicts_traj[0]['Epoch'][0][index]).time().strftime("%H:%M:%S") for index in timeTargetsIndices[0]],[pycdf.lib.tt2000_to_datetime(data_dicts_traj[1]['Epoch'][0][index]).time().strftime("%H:%M:%S") for index in timeTargetsIndices[1]]],
        'geoAlt': [[data_dicts_traj[0]['geoAlt'][0][index] for index in timeTargetsIndices[0]], [data_dicts_traj[1]['geoAlt'][0][index] for index in timeTargetsIndices[1]]],
        'geoLat': [[data_dicts_traj[0]['geoLat'][0][index] for index in timeTargetsIndices[0]], [data_dicts_traj[1]['geoLat'][0][index] for index in timeTargetsIndices[1]]],
        'geoLong': [[data_dicts_traj[0]['geoLong'][0][index] for index in timeTargetsIndices[0]], [data_dicts_traj[1]['geoLong'][0][index] for index in timeTargetsIndices[1]]],
        'geomagLat': [[data_dicts_traj[0]['geomagLat'][0][index] for index in timeTargetsIndices[0]], [data_dicts_traj[1]['geomagLat'][0][index] for index in timeTargetsIndices[1]]],
        'geomagLong': [[data_dicts_traj[0]['geomagLong'][0][index] for index in timeTargetsIndices[0]], [data_dicts_traj[1]['geomagLong'][0][index] for index in timeTargetsIndices[1]]],
        'geoLat_km': [[data_dicts_traj[0]['geoLat_km'][0][index] for index in timeTargetsIndices[0]], [data_dicts_traj[1]['geoLat_km'][0][index] for index in timeTargetsIndices[1]]],
        'geoLong_km': [[data_dicts_traj[0]['geoLong_km'][0][index] for index in timeTargetsIndices[0]], [data_dicts_traj[1]['geoLong_km'][0][index] for index in timeTargetsIndices[1]]],
        'geoILat': [[data_dicts_traj[0]['geoILat'][0][index] for index in timeTargetsIndices[0]],[data_dicts_traj[1]['geoILat'][0][index] for index in timeTargetsIndices[1]]],
        'geoILong': [[data_dicts_traj[0]['geoILong'][0][index] for index in timeTargetsIndices[0]], [data_dicts_traj[1]['geoILong'][0][index] for index in timeTargetsIndices[1]]],
        'geoILat_km':[[data_dicts_traj[0]['geoILat_km'][0][index] for index in timeTargetsIndices[0]],[data_dicts_traj[1]['geoILat_km'][0][index] for index in timeTargetsIndices[1]]],
        'geoILong_km':[[data_dicts_traj[0]['geoILong_km'][0][index] for index in timeTargetsIndices[0]],[data_dicts_traj[1]['geoILong_km'][0][index] for index in timeTargetsIndices[1]]]
    }
    Done(start_time)

    if plotALTvsLAT:

        from plottingStyles import ALTvsLAT

        ##################
        # --- PLOTTING ---
        ##################

        # --- Plot the Data ---

        prgMsg('Plotting ALTvsLAT Data')
        fig = plt.figure(figsize=ALTvsLAT.figure_size)

        #########################
        # --- geographic data ---
        #########################

        ax1 = fig.add_subplot(111)

        # Plot trajectories and UTC labels
        for i in range(len(dataPath_traj)):
            ax1.plot(data_dicts_traj[i]['geoLat'][0], data_dicts_traj[i]['geoAlt'][0], **ALTvsLAT.trajectory_style[i]) # Trajectories
            ax1.scatter(timeTargetsData['geoLat'][i], timeTargetsData['geoAlt'][i], **ALTvsLAT.marker_style[i]) # UTC Text

        # ticks major/minor
        ax1.minorticks_on()
        ax1.set_ylabel('Altitude (km)', **ALTvsLAT.axes_label_style[0])
        ax1.set_xlabel('Geographic Latitude', **ALTvsLAT.axes_label_style[0])
        ax1.xaxis.set_major_locator(MultipleLocator(ALTvsLAT.spacing_of_majorticks))
        ax1.xaxis.set_major_formatter('{x:.0f}')
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(**ALTvsLAT.majorTicks_style[0])
        ax1.tick_params(**ALTvsLAT.minorTicks_style[0])

        # Plot the text
        for j in range(2):
            for i, txt in enumerate(timeTargetsEpoch[j]):
                ax1.text(timeTargetsData['geoLat'][j][i] + ALTvsLAT.scatter_text_alignment_offsets[j][i],timeTargetsData['geoAlt'][j][i], txt, ha=ALTvsLAT.scatter_text_alignment[j][i], **ALTvsLAT.textUTC_style[j])

        # plot connecting lines
        if connectingLines:
            highFlyerOffset = 2
            for i in range(len(timeTargetsData['geoLat'][1])):
                x_values = [timeTargetsData['geoLat'][1][i], timeTargetsData['geoLat'][0][i+highFlyerOffset]]
                y_values = [timeTargetsData['geoAlt'][1][i], timeTargetsData['geoAlt'][0][i+highFlyerOffset]]
                ax1.plot(x_values, y_values, **ALTvsLAT.connecting_lines_style)

        ##########################
        # --- geomagnetic data ---
        ##########################
        ax2 = ax1.twiny()
        ax2.plot(data_dicts_traj[0]['geomagLat'][0], data_dicts_traj[0]['geoAlt'][0], **ALTvsLAT.plot_style_MAG) # High Flyer
        ax2.plot(data_dicts_traj[1]['geomagLat'][0], data_dicts_traj[1]['geoAlt'][0], **ALTvsLAT.plot_style_MAG)  # Low Flyer

        # ticks major/minor
        ax2.minorticks_on()
        ax2.xaxis.set_major_locator(MultipleLocator(ALTvsLAT.spacing_of_majorticks))
        ax2.xaxis.set_major_formatter('{x:.0f}')
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(**ALTvsLAT.majorTicks_style[1])
        ax2.tick_params(**ALTvsLAT.minorTicks_style[1])
        ax2.set_xlabel('Geomagnetic Latitude',**ALTvsLAT.axes_label_style[1])

        ax3 = ax1.twinx()
        ax3.get_yaxis().set_visible(False)
        slope = -2.378842439667985558915 # corresponds to line with -78.13deg inclination
        for i in range(11):
            ax3.axline(xy1=(69+i*0.5,-1*slope), xy2=(69 + (0.5*(i+1)), 0),**ALTvsLAT.mag_background_style)

        fig.savefig(r'D:\Data\ACESII\trajectories\trajectory_plots\Altitude_vs_lat_vs_maglat.png')
        Done(start_time)

    if plotLATvsLONG:
        from plottingStyles import LATvsLONG

        prgMsg('Plotting Lattitude vs Longitude')

        #######################
        # --- Lat/Long data ---
        #######################
        fig = plt.figure(figsize=LATvsLONG.figure_size)

        ax1 = fig.add_subplot(111)
        ax1.grid(True)

        # lat vs long
        for i in range(len(dataPath_traj)):
            ax1.plot(data_dicts_traj[i]['geoLat'][0],data_dicts_traj[i]['geoLong'][0],**LATvsLONG.trajectory_style[i])

        # ticks
        ax1.minorticks_on()
        ax1.xaxis.set_major_locator(MultipleLocator(LATvsLONG.spacing_of_majorticks))
        ax1.xaxis.set_major_formatter('{x:.0f}')
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.set_xlabel('Geographic Latitude', **LATvsLONG.axes_label_style[0])
        ax1.set_ylabel('Geographic Longitude', **LATvsLONG.axes_label_style[0])
        ax1.legend(['High Flyer','Low Flyer'])
        ax1.tick_params(**LATvsLONG.majorTicks_style[0])
        ax1.tick_params(**LATvsLONG.minorTicks_style[0])

        # Plot the Text
        for j in range(2):

            for i, txt in enumerate(timeTargetsEpoch[j]):
                # Markers
                ax1.scatter(timeTargetsData['geoLat'][j][i], timeTargetsData['geoLong'][j][i],**LATvsLONG.marker_style[j])
                # UTC labels
                ax1.text(timeTargetsData['geoLat'][j][i] + LATvsLONG.scatter_text_alignment_offsets[j][i],
                         timeTargetsData['geoLong'][j][i], txt, ha=LATvsLONG.scatter_text_alignment[j][i],
                         **LATvsLONG.textUTC_style[j])

        # ax 2
        ax2 = ax1.twiny()

        # geomagnetic lat vs long
        for i in range(len(dataPath_traj)):
            ax2.plot(data_dicts_traj[i]['geomagLat'][0], data_dicts_traj[i]['geoLong'][0], **LATvsLONG.plot_style_MAG)

        ax2.set_xlabel('Geomagnetic Latitude',**LATvsLONG.axes_label_style[1])

        # ticks major/minor
        ax2.minorticks_on()
        ax2.xaxis.set_major_locator(MultipleLocator(LATvsLONG.spacing_of_majorticks))
        ax2.xaxis.set_major_formatter('{x:.0f}')
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(**LATvsLONG.majorTicks_style[1])
        ax2.tick_params(**LATvsLONG.minorTicks_style[1])



        fig.savefig(r'D:\Data\ACESII\trajectories\trajectory_plots\Long_vs_lat_vs_maglat.png')

    if plotILATvsILONG:

        from plottingStyles import ILATvsILONG

        fig = plt.figure(figsize=ILATvsILONG.figureSize)
        ax1 = fig.add_subplot(111)
        ax1.minorticks_on()
        ax1.xaxis.set_major_formatter('{x:.0f}')

        ax1.tick_params(**ILATvsILONG.majorTicks_style[0])
        ax1.tick_params(**ILATvsILONG.minorTicks_style[0])

        plt.grid(True, **ILATvsILONG.gridMajor_style)
        plt.grid(True, **ILATvsILONG.gridMinor_style )

        # Plot the Trajectory Data
        if useKM:
            mod_file_ID = 'km'

            # AXIS LABELS
            ax1.set_xlabel(f'(+E|-W) [km]', **ILATvsILONG.axes_label_style)
            ax1.set_ylabel(f'(+N|-S) [km]', **ILATvsILONG.axes_label_style)

            # TICKS
            major_xticks = np.arange(-50, 50, ILATvsILONG.spacing_of_majorticks_km)
            major_yticks = np.arange(0, 550, ILATvsILONG.spacing_of_majorticks_km)
            ax1.set_xticks(major_xticks)
            ax1.set_yticks(major_yticks)
            ax1.xaxis.set_minor_locator(AutoMinorLocator())

            # GROUND TRACK TRAJECTS
            for i in range(len(dataPath_traj)):
                ax1.plot(data_dicts_traj[i]['geoLong_km'][0], data_dicts_traj[i]['geoLat_km'][0], **ILATvsILONG.trajectory_style[i])
                ax1.plot(data_dicts_traj[i]['geoILong_km'][0], data_dicts_traj[i]['geoILat_km'][0],**ILATvsILONG.intersections_style[i])

                # Plot X marks
                ax1.scatter(timeTargetsData['geoLong_km'][i], timeTargetsData['geoLat_km'][i], **ILATvsILONG.marker_style[i])

                # UTC Text
                for k in range(len(timeTargetsEpoch[i])):
                    ax1.text(timeTargetsData['geoLong_km'][i][k] + ILATvsILONG.adjusts_km[i][0],
                             timeTargetsData['geoLat_km'][i][k] + ILATvsILONG.adjusts_km[i][1],
                             timeTargetsEpoch[i][k], **ILATvsILONG.textUTC_style[i])

            if connectingLines:

                # plot the time tags for the High Flyer Only but plot the lines on the low flyer only
                for k in range(len(timeTargetsEpoch[1])):
                    ax1.text(timeTargetsData['geoILong_km'][0][k] + ILATvsILONG.adjusts_km[0][0],timeTargetsData['geoILat_km'][0][k] + ILATvsILONG.adjusts_km[0][1],timeTargetsEpoch[1][k], **ILATvsILONG.textUTC_intersect_style)

                    # Plot X marks
                    ax1.scatter(timeTargetsData['geoILong_km'][0][k], timeTargetsData['geoILat_km'][0][k], **ILATvsILONG.marker_intersect_style)

                    # Plot the connecting lines
                    x_values = [timeTargetsData['geoILong_km'][0][k],timeTargetsData['geoILong_km'][1][k]]
                    y_values = [timeTargetsData['geoILat_km'][0][k], timeTargetsData['geoILat_km'][1][k]]

                    ax1.plot(x_values, y_values, **ILATvsILONG.connecting_lines_style)

        # TYPICALLY LAT/LONG
        else:
            mod_file_ID = 'latlong'
            ax1.set_xlabel('Geographic Long', **ILATvsILONG.axes_label_style)
            ax1.set_ylabel('Geographic Lat', **ILATvsILONG.axes_label_style)

            # TICKS
            major_xticks = np.arange(-50, 50, ILATvsILONG.spacing_of_majorticks)
            major_yticks = np.arange(0, 550, ILATvsILONG.spacing_of_majorticks)
            ax1.set_xticks(major_xticks)
            ax1.set_yticks(major_yticks)
            ax1.xaxis.set_minor_locator(AutoMinorLocator())

            # GROUND TRACK TRAJECTS
            for i in range(len(dataPath_traj)):
                ax1.plot(data_dicts_traj[i]['geoLong'][0], data_dicts_traj[i]['geoLat'][0],
                         **ILATvsILONG.trajectory_style[i])
                ax1.plot(data_dicts_traj[i]['geoILong'][0], data_dicts_traj[i]['geoILat'][0],
                         **ILATvsILONG.intersections_style[i])

                # Plot X marks
                ax1.scatter(timeTargetsData['geoLong'][i], timeTargetsData['geoLat'][i],
                            **ILATvsILONG.marker_style[i])

                # UTC Text
                for k in range(len(timeTargetsEpoch[i])):
                    ax1.text(timeTargetsData['geoLong'][i][k] + ILATvsILONG.adjusts[i][0],
                             timeTargetsData['geoLat'][i][k] + ILATvsILONG.adjusts[i][1],
                             timeTargetsEpoch[i][k], **ILATvsILONG.textUTC_style[i])

                # GEOMAGNETIC AXIS
                ax2 = ax1.twinx()
                ax2.plot(data_dicts_traj[0]['geoLong'][0], data_dicts_traj[0]['geomagLat'][0],**ILATvsILONG.plot_style_MAG)
                # ticks major/minor
                ax2.minorticks_on()
                ax2.yaxis.set_major_locator(MultipleLocator(ILATvsILONG.spacing_of_majorticks))
                ax2.yaxis.set_major_formatter('{x:.0f}')
                ax2.yaxis.set_minor_locator(AutoMinorLocator())
                ax2.tick_params(**ILATvsILONG.majorTicks_style[1])
                ax2.tick_params(**ILATvsILONG.minorTicks_style[1])
                ax2.set_ylabel('Geomagnetic Lat', **ILATvsILONG.axes_label_style)

            if connectingLines:

                # plot the time tags for the High Flyer Only but plot the lines on the low flyer only
                for k in range(len(timeTargetsEpoch[1])):
                    ax1.text(timeTargetsData['geoILong'][0][k] + ILATvsILONG.adjusts[0][0],
                             timeTargetsData['geoILat'][0][k] + ILATvsILONG.adjusts[0][1], timeTargetsEpoch[1][k],
                             **ILATvsILONG.textUTC_intersect_style)

                    # Plot X marks
                    ax1.scatter(timeTargetsData['geoILong'][0][k], timeTargetsData['geoILat'][0][k],
                                **ILATvsILONG.marker_intersect_style)

                    # Plot the connecting lines
                    x_values = [timeTargetsData['geoILong'][0][k], timeTargetsData['geoILong'][1][k]]
                    y_values = [timeTargetsData['geoILat'][0][k], timeTargetsData['geoILat'][1][k]]

                    ax1.plot(x_values, y_values, **ILATvsILONG.connecting_lines_style)

        ax1.legend()
        plt.title(f'{wInstr_projection} \n B-Field Projection Altitude: {ILATvsILONG.targetProjectionAltitude} km \n Origin: {ILATvsILONG.refLat}$^\circ$N {ILATvsILONG.refLong}$^\circ$E')
        fig.savefig(rf'D:\Data\ACESII\trajectories\trajectory_plots\BProjection_{mod_file_ID}.png')

    if plotESAoverlay:

        from plottingStyles import ESAoverlay

        data_dicts_ESA = data_dicts_ESA[wInstr_overlay]



        prgMsg('Plotting ESA overlay')

        #############################
        # --- "Zoom" into Dataset ---
        #############################

        # index locations of specific latitude
        start = np.abs(data_dicts_traj[wRocket_overlay-4]['geoLat'][0] -ESAoverlay.start[wRocket_overlay-4]).argmin()
        end = np.abs(data_dicts_traj[wRocket_overlay-4]['geoLat'][0] -ESAoverlay.end[wRocket_overlay-4]).argmin()

        # Get a subset of the data
        Energy = data_dicts_ESA[0]['Energy'][0]
        Pitch = data_dicts_ESA[0]['Pitch_Angle'][0]
        EpochData = data_dicts_ESA[wRocket_overlay-4]['Epoch_esa'][0][start:end]
        ESAData = data_dicts_ESA[wRocket_overlay-4][wInstr_overlay][0][start:end]

        ##########################
        # --- PLOT Spectrogram ---
        ##########################

        fig, ax = plt.subplots(3)

        fig.set_figheight(ESAoverlay.figure_height)
        fig.set_figwidth(ESAoverlay.figure_width)

        for i in range(len(wPitches_overlay)):

            # Trajectory and ESA data
            X = data_dicts_traj[wRocket_overlay-4]['geoLat'][0][start:end]
            Y = Energy
            Z = np.array([ESAData[j][wPitches_overlay[i]] for j in range(len(ESAData))]).transpose() # Get a slice in pitch data

            # PCOLORMESH
            cmap = ax[i].pcolormesh(X, Y, Z, **ESAoverlay.pcolormesh_style)
            cbar = plt.colorbar(cmap, **ESAoverlay.colorbar_style,ax = ax[i])
            cbar.ax.tick_params(**ESAoverlay.majorTicks_style)
            cbar.set_label('Counts',**ESAoverlay.cbar_label_style)
            ax[i].set_ylabel('Energy [eV]', **ESAoverlay.axes_label_style)

            if i == 2:
                ax[i].set_xlabel('Geodetic Lat [deg]', **ESAoverlay.axes_label_style)

            ax[i].set_ylim(**ESAoverlay.yaxis_lim)

            # title
            ax[i].set_title(f'{wInstr_overlay}, Pitch = {Pitch[wPitches_overlay[i]]}$^\circ$', **ESAoverlay.title_style)

            # margins
            ax[i].margins(**ESAoverlay.margins_style)


            # TICKS
            spacing_of_majorticks = 1
            ax[i].xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
            ax[i].xaxis.set_minor_locator(AutoMinorLocator())
            ax[i].minorticks_on()
            ax[i].xaxis.set_major_formatter('{x:.0f}')
            ax[i].tick_params(**ESAoverlay.majorTicks_style)
            ax[i].tick_params(**ESAoverlay.minorTicks_style)

            # Alt AXIS
            spacing_of_majorticks = 1
            ax1 = ax[i].twinx()
            ax1.minorticks_on()
            ax1.plot(data_dicts_traj[wRocket_overlay-4]['geoLat'][0][start:end], data_dicts_traj[wRocket_overlay-4]['geoAlt'][0][start:end], **ESAoverlay.trajectory_style[wRocket_overlay - 4])
            ax1.set_ylabel('Alt [km]', **ESAoverlay.axes_label_style)
            ax1.minorticks_on()
            ax1.xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
            ax1.xaxis.set_major_formatter('{x:.0f}')
            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(**ESAoverlay.minorTicks_style)
            ax1.tick_params(**ESAoverlay.minorTicks_style)
            plt.tight_layout()

        fig.savefig(rf'D:\Data\ACESII\trajectories\trajectory_plots\{fliers[wRocket_overlay - 4]}\{wInstr_overlay}_overlay_ALTvsLAT.png')
        Done(start_time)

    if plotESAmovie:
        wInstr_movie = 'eepaa'

        if create_movieESAplots:

            prgMsg('Creating ESA plots for movie \n')

            if plotDistFunc:
                data_dicts_movie = data_dicts_dist[wInstr_movie]
                wInstr_movie = 'Distribution_Function'
            else:
                wInstr_movie = wInstr_movie
                data_dicts_movie = data_dicts_ESA[wInstr_movie]



            from plottingStyles import ESAmovie

            # Get the data
            Energy = data_dicts_movie[0]['Energy'][0]
            Pitch = data_dicts_movie[0]['Pitch_Angle'][0]
            geoAlt = [data_dicts_traj[0]['geoAlt'][0],data_dicts_traj[1]['geoAlt'][0]]
            geoLat = [data_dicts_traj[0]['geoLat'][0],data_dicts_traj[1]['geoLat'][0]]
            geoLong = [data_dicts_traj[0]['geoLong'][0],data_dicts_traj[1]['geoLong'][0]]
            EpochData = [np.array(data_dicts_movie[0]['Epoch_esa'][0]),np.array(data_dicts_movie[1]['Epoch_esa'][0])]
            EpochData_dates = [np.array([pycdf.lib.tt2000_to_datetime(EpochData[0][i]).strftime("%H:%M:%S:%f") for i in range(len(EpochData[0]))]),
                               np.array([pycdf.lib.tt2000_to_datetime(EpochData[1][i]).strftime("%H:%M:%S:%f") for i in range(len(EpochData[1]))]),]
            ESAData =  [data_dicts_movie[0][wInstr_movie][0],data_dicts_movie[1][wInstr_movie][0]]

            ########################################
            # --- ADD FILLVALS TO LOW FLYER DATA ---
            ########################################
            highFlyerSize = len(geoAlt[0])
            lowFlyerSize = len(geoAlt[1])

            # --- Append start Values ---
            no_of_points_start = int((EpochData[1][0] - EpochData[0][0]) / (rocketAttrs.MinorFrameTime))
            newAlt = [geoAlt[1][0] for i in range(no_of_points_start)]
            newLat = [geoLat[1][0] for i in range(no_of_points_start)]
            ESAfillVal = [[0 for engy in range(len(Energy))] for ptch in range(len(Pitch)) ]
            newESA = [ESAfillVal for i in range(no_of_points_start)]

            for i in range(len(geoAlt[1])):
                newAlt.append(geoAlt[1][i])
                newLat.append(geoLat[1][i])
                newESA.append(ESAData[1][i])

            # --- Append the ending values ---
            remainingIndicies = highFlyerSize - (lowFlyerSize + no_of_points_start)

            for i in range(remainingIndicies):
                newAlt.append(geoAlt[1][-1])
                newLat.append(geoLat[1][-1])
                newESA.append(ESAfillVal)

            geoAlt[1] = np.array(newAlt)
            geoLat[1] = np.array(newLat)
            ESAData[1] = np.array(newESA)


            ################################
            # --- GENERATE ALL THE PLOTS ---
            ################################

            locations = [i for i in range(0,len(EpochData_dates[0]),frame_skips)]
            # locations = [6700] # NEEDS TO BE THE HIGH FLYER LENGTH

            for location in tqdm(locations):

                fig = plt.figure()
                fig.set_figwidth(ESAmovie.figure_width)
                fig.set_figheight(ESAmovie.figure_height)

                gs = fig.add_gridspec(nrows = 2, ncols =1,height_ratios=[1,6])
                gsAlt = gs[0].subgridspec(nrows =1, ncols = 1)
                gsESA = gs[1].subgridspec(nrows = 2, ncols= 2,hspace=0.1,height_ratios=[1,1],    width_ratios=[90,1])

                # Altitude Plot
                ax1 = fig.add_subplot(gsAlt[0,:])

                # # ESA Data plots
                ax2 = fig.add_subplot(gsESA[0,0])
                ax3 = fig.add_subplot(gsESA[1,0],sharex=ax2,sharey=ax2)
                ax_empty = fig.add_subplot(gsESA[0,1])

                # Give colorbar its own axis but size 2
                ax_colorbar = fig.add_subplot(gsESA[0:,1])

                ax_empty.remove()# Doing this resizes the colorbar to be worth 2 sizes in length

                # updating title
                fig.suptitle(f'{wInstr_movie.upper()}\n{EpochData_dates[0][location]}',**ESAmovie.title_style)

                # --- --- --- --- --- --- -
                # --- PLOT1: TRAJECTORY ---
                # --- --- --- --- --- --- -

                for i in range(2):
                    # Plot the trajectory
                    ax1.plot(geoLat[i], geoAlt[i], **ESAmovie.trajectory_style[i])  # Trajectories

                    # Plot the specific location we're at based on index
                    ax1.scatter(geoLat[i][location], geoAlt[i][location], **ESAmovie.marker_style[i])

                    # Plot guiding latitude

                    # Plot the altitude text
                    ax1.text(geoLat[i][location] - ESAmovie.textOffset, geoAlt[i][location], f'{round(geoAlt[i][location],ESAmovie.altrounding)} km' , ha=ESAmovie.text_alignment, **ESAmovie.textUTC_style[i])

                ax1.set_ylabel('Altitude (km)', **ESAmovie.axes_label_traj_style)
                ax1.set_xlabel('Geographic Latitude', **ESAmovie.axes_label_traj_style)
                ax1.xaxis.set_major_locator(MultipleLocator(ESAmovie.spacing_of_majorticks_Traj))
                ax1.xaxis.set_major_formatter('{x:.0f}')
                ax1.xaxis.set_minor_locator(AutoMinorLocator())
                ax1.tick_params(**ESAmovie.majorTicks_style)
                ax1.tick_params(**ESAmovie.minorTicks_style)

                ####################################
                # --- PLOT2: HIGH FLYER ESA DATA ---
                ####################################
                if polar:
                    Theta = [ptch*np.pi/180 for ptch in Pitch]
                    Radii = Energy[ESAmovie.EnergyStart:ESAmovie.EnergyEnd]
                    Z = np.array([ESAData[0][location][ptch][ESAmovie.EnergyStart:ESAmovie.EnergyEnd] for ptch in range(len(Pitch))]).transpose()  # Get a slice in pitch data

                    plt.grid(False)
                    rtick_locs = Radii
                    rtick_labels = ['%.f' % r for r in rtick_locs]
                    ax2.set_rgrids(rtick_locs, rtick_labels)
                    ax2.set_thetalim(thetamin=-10, thetamax=190)
                    cmap = ax2.pcolormesh(Pitch, Radii, Z,**ESAmovie.pcolormesh_style)
                    ax2.set_thetagrids([theta * 10 - 10 for theta in range(210 // 10)])
                    plt.grid(True)
                else:
                    X = Pitch
                    Y = Energy[ESAmovie.EnergyStart:ESAmovie.EnergyEnd]
                    Z = np.array([ESAData[0][location][ptch][ESAmovie.EnergyStart:ESAmovie.EnergyEnd] for ptch in range(len(Pitch))]).transpose()  # Get a slice in pitch data
                    cmap = ax2.pcolormesh(X, Y, Z,**ESAmovie.pcolormesh_style)
                    ax2.tick_params(**ESAmovie.majorTicks_style)
                    ax2.tick_params(**ESAmovie.minorTicks_style)
                    ax2.set_ylabel('Energy [eV]', **ESAmovie.axes_label_ESA_style)
                    ax2.xaxis.set_major_locator(MultipleLocator(ESAmovie.spacing_of_majorticks_ESA))
                    ax2.xaxis.set_major_formatter('{x:.0f}')

                ############################
                # --- LOW FLYER ESA DATA ---
                ############################

                if polar:
                    Theta = [ptch * np.pi / 180 for ptch in Pitch]
                    Radii = Energy[ESAmovie.EnergyStart:ESAmovie.EnergyEnd]
                    Z = np.array([ESAData[1][location][ptch][ESAmovie.EnergyStart:ESAmovie.EnergyEnd] for ptch in range(len(Pitch))]).transpose()  # Get a slice in pitch data

                    plt.grid(False)
                    rtick_locs = Radii
                    rtick_labels = ['%.f' % r for r in rtick_locs]
                    ax3.set_rgrids(rtick_locs, rtick_labels)
                    ax3.set_thetalim(thetamin=-10, thetamax=190)
                    cmap = ax3.pcolormesh(Pitch, Radii, Z,**ESAmovie.pcolormesh_style)
                    ax3.set_thetagrids([theta * 10 - 10 for theta in range(210 // 10)])
                    plt.grid(True)
                else:
                    X= Pitch
                    Y = Energy[ESAmovie.EnergyStart:ESAmovie.EnergyEnd]
                    Z = np.array([ESAData[1][location][ptch][ESAmovie.EnergyStart:ESAmovie.EnergyEnd] for ptch in range(len(Pitch))]).transpose()  # Get a slice in pitch data

                    ax3.set_yscale('log')
                    ax3.set_xlabel('Pitch [deg]',**ESAmovie.axes_label_ESA_style)
                    ax3.set_ylabel('Energy [eV]',**ESAmovie.axes_label_ESA_style)
                    ax3.tick_params(**ESAmovie.majorTicks_style)
                    ax3.tick_params(**ESAmovie.minorTicks_style)
                    ax3.xaxis.set_major_locator(MultipleLocator(ESAmovie.spacing_of_majorticks_ESA))
                    ax3.xaxis.set_major_formatter('{x:.0f}')
                    cmap = ax3.pcolormesh(X, Y, Z,**ESAmovie.pcolormesh_style)

                ##################
                # --- COLORBAR ---
                ##################
                cbar = plt.colorbar(cmap, **ESAmovie.colorbar_style, cax=ax_colorbar)
                cbar.minorticks_on()
                cbar.ax.tick_params(**ESAmovie.majorTicks_style)

                if plotDistFunc:
                    cbar.set_label('Distribution Function [$m^{-6}s^{3}$]', **ESAmovie.cbar_label_style)
                else:
                    cbar.set_label('Counts', **ESAmovie.cbar_label_style)

                fig.savefig(f'D:\Data\ACESII\\trajectories\\trajectory_plots\\photos_for_movie\\{location}.png')
                plt.close(fig)


            Done(start_time)

        #########################
        # --- BUILD THE MOVIE ---
        #########################
        prgMsg('Building the Movie')
        from moviepy.editor import ImageSequenceClip

        # Get the images
        pathToImages = 'D:\Data\ACESII\\trajectories\\trajectory_plots\\photos_for_movie\\'
        images = glob(f'{pathToImages}*.png*')

        # sort the images
        images_sorted = []
        for i, txt in enumerate(images):
            images_sorted.append(int(txt.replace(pathToImages,'').replace('.png','')))

        images_sorted.sort()

        images = [ pathToImages + str(num) + '.png' for num in images_sorted]


        # creating a Image sequence clip with fps = 1
        animation = ImageSequenceClip(images, fps=ESAmovie.fps)

        # showing  clip
        animation.write_videofile(f'D:\Data\ACESII\\trajectories\\trajectory_plots\\ACESII_{wInstr_movie}.mp4', fps=20)

        Done(start_time)




# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
ACESIIplotting()
