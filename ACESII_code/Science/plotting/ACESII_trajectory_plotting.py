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

# --- Alt vs Lat/geoLat ---
plotALTvsLAT = False

# --- Lat vs Long ---
plotLATvsLONG = False

# --- B-Field Projection ---
plotILATvsILONG = False
useKM = False
wInstr_projection = 'eepaa'

# --- ESA overlay ---
plotESAoverlay = False
wRocket_overlay = 4
wPitches_overlay = [1, 10, 19] # 0deg, 90deg, 180deg
wInstr_overlay = 'eepaa'

# --- ESA Movie ---
plotESAmovie = True
frame_skips = 2 # 1 is no skip i.e. all frames used, 5 is use only every 5th frame, 10 etc...
wInstr_movie = 'eepaa'
plotSpecificLocations = False
specific_locations = [6000 + i for i in range(0,100)]

# -- parameters involved in plotting the interpolated distribution function --
wDataPlot = [1] # 0 - Counts, 1 - DiffEFlux, 2 - Dist. Func
N = 200 # number of points for NxN interpolation grid
fillvalue = 0 # interpolated plots fillvalue
normalizeToThermal = True # normalize ESA plots to the electron thermal velocity, Te = 0.1

# Use connecting lines in all the relevant plots?
connectingLines = True


# --- Plot the all sky data ---
plotAllSkyData = True

# --- --- --- ---
# --- import ---
# --- --- --- ---
import numpy as np
from matplotlib import pyplot as plt, animation
from matplotlib.patches import Wedge
from ACESII_code.class_var_func import m_e, q0
from copy import deepcopy
from tqdm import tqdm
from glob import glob
from scipy.interpolate import LinearNDInterpolator
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from ACESII_code.data_paths import fliers, ACES_data_folder
from ACESII_code.class_var_func import color, Rotation3D, loadDictFromFile
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
                try:
                    data_dict['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_esa'][0][i]) for i in range(len(data_dict['Epoch_esa'][0]))])
                except:
                    data_dict['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][i]) for i in range(len(data_dict['Epoch'][0]))])

                parentDict[pkey].append(data_dict)
            except:
                parentDict[pkey].append(data_dict)

    return parentDict

def ACESIIplotting():

    rocketAttrs,b,c = ACES_mission_dicts()

    if plotAllSkyData == False:
        # Trajectory Data
        trajFolderPath = f'{ACES_data_folder}trajectories\\'
        dataPath_traj = [glob(trajFolderPath + rf'{fliers[0]}\\*_ILat_ILong*'),
                         glob(trajFolderPath + rf'{fliers[1]}\\\\*_ILat_ILong*')]

        # L1 ESA data
        FolderPath = f'{ACES_data_folder}L1\\'
        dataPath_L1ESA = {'eepaa':[glob(FolderPath + rf'{fliers[0]}\\*eepaa_*'), glob(FolderPath + rf'{fliers[1]}\\\\*eepaa_*')],
                          'iepaa':[glob(FolderPath + rf'{fliers[0]}\\*iepaa*'),glob(FolderPath + rf'{fliers[1]}\\\\*iepaa*')],
                          'leesa':[glob(FolderPath + rf'{fliers[0]}\\*leesa*'),glob(FolderPath + rf'{fliers[1]}\\\\*leesa*')]}
        # L2 ESA data
        FolderPath = f'{ACES_data_folder}L2\\'
        dataPath_L2ESA = {'eepaa': [glob(FolderPath + rf'{fliers[0]}\\*eepaa_fullCal*'), glob(FolderPath + rf'{fliers[1]}\\\\*eepaa_fullCal*')],
                          'iepaa': [glob(FolderPath + rf'{fliers[0]}\\*iepaa*'), glob(FolderPath + rf'{fliers[1]}\\\\*iepaa*')],
                          'leesa': [glob(FolderPath + rf'{fliers[0]}\\*leesa*'), glob(FolderPath + rf'{fliers[1]}\\\\*leesa*')]}

        # Distribution Function Data
        FolderPath = f'{ACES_data_folder}\science\DistFunc\\'
        dataPath_Dist = {'eepaa':[glob(FolderPath + rf'{fliers[0]}\\*eepaa_fullCal*'), glob(FolderPath + rf'{fliers[1]}\\\\*eepaa_fullCal*')],
                          'iepaa':[glob(FolderPath + rf'{fliers[0]}\\*iepaa*'),glob(FolderPath + rf'{fliers[1]}\\\\*iepaa*')],
                          'leesa':[glob(FolderPath + rf'{fliers[0]}\\*leesa*'),glob(FolderPath + rf'{fliers[1]}\\\\*leesa*')]}

        prgMsg('Collecting Trajectory data')

        # --- Create TRAJECTORY Data Dicts ---
        data_dicts_traj = []

        for i in range(len(dataPath_traj)):
            data_dict_traj = {}
            data_dict_traj = loadDictFromFile(dataPath_traj[i][0],data_dict_traj)
            data_dict_traj['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch_esa'][0][i]) for i in range(len(data_dict_traj['Epoch_esa'][0]))])
            data_dicts_traj.append(data_dict_traj)

        Done(start_time)

        # --- Create ESA L1 Data Dicts ---
        prgMsg('Collecting ESA L1 data')
        data_dicts_template = {'eepaa':[],'iepaa':[],'leesa':[]}
        data_dicts_ESA_L1 = storeData(data_dicts_template,dataPath_L1ESA)
        Done(start_time)

        # --- Create ESA L2 Data Dicts ---
        prgMsg('Collecting ESA L2 data')
        data_dicts_template = {'eepaa': [], 'iepaa': [], 'leesa': []}
        data_dicts_ESA_L2 = storeData(data_dicts_template, dataPath_L2ESA)
        Done(start_time)

        # --- Create Distribution_Func Data Dicts ---
        prgMsg('Collecting Distribution Function data')
        data_dicts_template = {'eepaa': [], 'iepaa': [], 'leesa': []}
        data_dicts_dist = storeData(data_dicts_template,dataPath_Dist)

        ####################################################################
        # --- find the points 100s, 200s ... 600s in the Trajectory data ---
        ####################################################################
        prgMsg('Finding Target Times')
        timeTargets = [[100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600], [100, 150, 200, 250, 300, 350, 400]]

        timeTargetsIndices = [[], []]
        timeTargetsEpoch = []
        for i in range(len(timeTargets)):
            for timeTarg in timeTargets[i]:
                timeTargetsIndices[i].append(np.abs(data_dicts_traj[i]['Epoch_esa'][0] - (rocketAttrs.Launch_Times[i] + timeTarg * (10 ** (9)))).argmin())

            target_in_Epoch = [pycdf.lib.tt2000_to_datetime(data_dicts_traj[i]['Epoch_esa'][0][index]) for index in timeTargetsIndices[i]]
            timeTargetsEpoch.append([targetTime.strftime("%H:%M:%S") for targetTime in target_in_Epoch])

        timeTargetsData = {
            'Epoch': [[pycdf.lib.tt2000_to_datetime(data_dicts_traj[0]['Epoch_esa'][0][index]).time().strftime("%H:%M:%S") for index in timeTargetsIndices[0]],
                      [pycdf.lib.tt2000_to_datetime(data_dicts_traj[1]['Epoch_esa'][0][index]).time().strftime("%H:%M:%S") for index in timeTargetsIndices[1]]],
            'geoAlt': [[data_dicts_traj[0]['geoAlt'][0][index] for index in timeTargetsIndices[0]],
                       [data_dicts_traj[1]['geoAlt'][0][index] for index in timeTargetsIndices[1]]],
            'geoLat': [[data_dicts_traj[0]['geoLat'][0][index] for index in timeTargetsIndices[0]],
                       [data_dicts_traj[1]['geoLat'][0][index] for index in timeTargetsIndices[1]]],
            'geoLong': [[data_dicts_traj[0]['geoLong'][0][index] for index in timeTargetsIndices[0]],
                        [data_dicts_traj[1]['geoLong'][0][index] for index in timeTargetsIndices[1]]],
            'geomagLat': [[data_dicts_traj[0]['geomagLat'][0][index] for index in timeTargetsIndices[0]],
                          [data_dicts_traj[1]['geomagLat'][0][index] for index in timeTargetsIndices[1]]],
            'geomagLong': [[data_dicts_traj[0]['geomagLong'][0][index] for index in timeTargetsIndices[0]],
                           [data_dicts_traj[1]['geomagLong'][0][index] for index in timeTargetsIndices[1]]],
            'geoLat_km': [[data_dicts_traj[0]['geoLat_km'][0][index] for index in timeTargetsIndices[0]],
                          [data_dicts_traj[1]['geoLat_km'][0][index] for index in timeTargetsIndices[1]]],
            'geoLong_km': [[data_dicts_traj[0]['geoLong_km'][0][index] for index in timeTargetsIndices[0]],
                           [data_dicts_traj[1]['geoLong_km'][0][index] for index in timeTargetsIndices[1]]],
            'geoILat': [[data_dicts_traj[0]['geoILat'][0][index] for index in timeTargetsIndices[0]],
                        [data_dicts_traj[1]['geoILat'][0][index] for index in timeTargetsIndices[1]]],
            'geoILong': [[data_dicts_traj[0]['geoILong'][0][index] for index in timeTargetsIndices[0]],
                         [data_dicts_traj[1]['geoILong'][0][index] for index in timeTargetsIndices[1]]],
            'geoILat_km': [[data_dicts_traj[0]['geoILat_km'][0][index] for index in timeTargetsIndices[0]],
                           [data_dicts_traj[1]['geoILat_km'][0][index] for index in timeTargetsIndices[1]]],
            'geoILong_km': [[data_dicts_traj[0]['geoILong_km'][0][index] for index in timeTargetsIndices[0]],
                            [data_dicts_traj[1]['geoILong_km'][0][index] for index in timeTargetsIndices[1]]]}
        Done(start_time)
    else:
        # --- Create Attitude Data Dict ---
        prgMsg('Collecting Attitude data')
        FolderPath = f'{ACES_data_folder}\\attitude\\'
        dataPath_attitude = {'attitude': [glob(FolderPath + rf'{fliers[0]}\\*Attitude_Solution*'), glob(FolderPath + rf'{fliers[1]}\\\\*Attitude_Solution*')]}
        data_dicts_template = {'attitude':[]}
        data_dicts_attitude = storeData(data_dicts_template,dataPath_attitude)

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

        data_dicts_ESA = data_dicts_ESA_L1[wInstr_overlay]

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

            # PCOLORMESH Plot
            cmap = ax[i].pcolormesh(X, Y, Z, **ESAoverlay.pcolormesh_style)
            cbar = plt.colorbar(cmap, ax = ax[i])
            cbar.ax.tick_params(**ESAoverlay.majorTicks_style_cbar)
            cbar.set_label('Counts', **ESAoverlay.cbar_label_style)
            ax[i].set_ylabel('Energy [eV]', **ESAoverlay.axes_label_style)

            ax[i].minorticks_on()
            ax[i].yaxis.set_tick_params(**ESAoverlay.tick_style)
            ax[i].xaxis.set_tick_params(**ESAoverlay.tick_style)

            if i == 2:
                ax[i].set_xlabel('Geodetic Lat [deg]', **ESAoverlay.axes_label_style)

            ax[i].set_ylim(**ESAoverlay.yaxis_lim)

            # title
            ax[i].set_title(f'{wInstr_overlay}, Pitch = {Pitch[wPitches_overlay[i]]}$^\circ$', **ESAoverlay.title_style)

            ##################
            # --- Alt AXIS ---
            ##################
            ax1 = ax[i].twinx()
            ax1.plot(data_dicts_traj[wRocket_overlay - 4]['geoLat'][0][start:end], data_dicts_traj[wRocket_overlay - 4]['geoAlt'][0][start:end], **ESAoverlay.trajectory_style_twin[wRocket_overlay - 4])
            ax1.set_ylabel('Alt [km]', **ESAoverlay.axes_label_style_twin)
            ax1.minorticks_on()
            ax1.yaxis.set_tick_params(**ESAoverlay.tick_style)
            ax1.xaxis.set_tick_params(**ESAoverlay.tick_style)
            plt.tight_layout()

        fig.savefig(rf'D:\Data\ACESII\trajectories\trajectory_plots\{fliers[wRocket_overlay - 4]}\{wInstr_overlay}_overlay_ALTvsLAT.png')
        Done(start_time)

    if plotESAmovie:
        wInstr_movie = 'eepaa'

        prgMsg('Creating ESA movie (This may take awhile) \n')

        if wDataPlot[0] == 2: # Distribution Function
            data_dicts_movie = data_dicts_dist[wInstr_movie] # DONT CHANGE THE ORDER OF THIS LINE AND THE NEXT LINE
            wInstr_movie = 'Distribution_Function'
        elif wDataPlot[0] == 1:
            data_dicts_movie = data_dicts_ESA_L2[wInstr_movie]  # DONT CHANGE THE ORDER OF THIS LINE AND THE NEXT LINE
            wInstr_movie = 'Differential_Energy_Flux'
        else:
            wInstr_movie = wInstr_movie
            data_dicts_movie = data_dicts_ESA_L1[wInstr_movie]

        from plottingStyles import ESAmovie

        # Get the data
        Energy = data_dicts_movie[0]['Energy'][0]
        Pitch = data_dicts_movie[0]['Pitch_Angle'][0]
        geoAlt = [data_dicts_traj[0]['geoAlt'][0],data_dicts_traj[1]['geoAlt'][0]]
        geoLat = [data_dicts_traj[0]['geoLat'][0],data_dicts_traj[1]['geoLat'][0]]
        geoLong = [data_dicts_traj[0]['geoLong'][0],data_dicts_traj[1]['geoLong'][0]]
        EpochData = [np.array(data_dicts_movie[0]['Epoch_esa'][0]),np.array(data_dicts_movie[1]['Epoch_esa'][0])]
        EpochData_dates = [np.array([pycdf.lib.tt2000_to_datetime(EpochData[0][i]).strftime("%H:%M:%S:%f") for i in range(len(EpochData[0]))]), np.array([pycdf.lib.tt2000_to_datetime(EpochData[1][i]).strftime("%H:%M:%S:%f") for i in range(len(EpochData[1]))]),]
        ESAData = [data_dicts_movie[0][wInstr_movie][0],data_dicts_movie[1][wInstr_movie][0]]

        ##################################################
        # --- RESIZE LOW FLYER DATA TO MATCH HIGHFLYER ---
        ##################################################
        highFlyerSize, lowFlyerSize = len(geoAlt[0]),len(geoAlt[1])

        # --- Append start Values ---
        no_of_points_start = int((EpochData[1][0] - EpochData[0][0]) / (rocketAttrs.MinorFrameTime))
        newAlt = [geoAlt[1][0] for i in range(no_of_points_start)]
        newLat = [geoLat[1][0] for i in range(no_of_points_start)]
        newLong = [geoLong[1][0] for i in range(no_of_points_start)]
        ESAfillVal = [[0 for engy in range(len(Energy))] for ptch in range(len(Pitch)) ]
        newESA = [ESAfillVal for i in range(no_of_points_start)]

        for i in range(len(geoAlt[1])):
            newAlt.append(geoAlt[1][i])
            newLat.append(geoLat[1][i])
            newESA.append(ESAData[1][i])
            newLong.append(geoLong[1][i])

        # --- Append the ending values ---
        remainingIndicies = highFlyerSize - (lowFlyerSize + no_of_points_start)

        for i in range(remainingIndicies):
            newAlt.append(geoAlt[1][-1])
            newLat.append(geoLat[1][-1])
            newLong.append(geoLong[1][-1])
            newESA.append(ESAfillVal)

        geoAlt[1], geoLat[1], geoLong[1] = np.array(newAlt), np.array(newLat), np.array(newLong)
        ESAData[1] = np.array(newESA)

        ##################################
        # --- FORMAT DATA FOR PLOTTING ---
        ##################################

        if wDataPlot[0] in [1,2]: # Handles Distribution Function and DiffEFlux case

            # --- Velocity Space coordinate transformation ---
            modelData = ESAData[0]
            scienceVar = [[], []]

            # The vperp and vparallel coordinates never change, so pre-calculate them here for use later
            Vpars = np.array([np.cos(np.radians(Pitch[ptch])) * np.sqrt(2 * q0 * Energy[engy] / m_e) for engy in
                        range(ESAmovie.EnergyStart, ESAmovie.EnergyEnd) for ptch in range(len(modelData[0]))])
            Vperps = np.array([np.sin(np.radians(Pitch[ptch])) * np.sqrt(2 * q0 * Energy[engy] / m_e) for engy in
                         range(ESAmovie.EnergyStart, ESAmovie.EnergyEnd) for ptch in range(len(modelData[0]))])

            for i in range(2):

                if plotSpecificLocations:
                    iterateThis = specific_locations
                else:
                    iterateThis = range(len(modelData))

                for tme in iterateThis:
                    scienceVar[i].append([ESAData[i][tme][ptch][engy] for engy in range(ESAmovie.EnergyStart, ESAmovie.EnergyEnd) for ptch in range(len(Pitch))])

            scienceVar = np.array(scienceVar)

            # --- Interpolate and store the Data for each flyer ----

            # Note: I'm not interested in interpolating over the entire Vspace range. Much of the interesting stuff happens <2000eV
            # Reduce the X,Y linspace to only relvant data

            if normalizeToThermal:
                if wDataPlot[0] == 2:
                    xlimitBot = ESAmovie.normalizedXLimits_dist[0] * ESAmovie.vth_e
                    xlimitTop = ESAmovie.normalizedXLimits_dist[1] * ESAmovie.vth_e
                    ylimitBot = ESAmovie.normalizedYLimits_dist[0] * ESAmovie.vth_e
                    ylimitTop = ESAmovie.normalizedYLimits_dist[1] * ESAmovie.vth_e
                elif wDataPlot[0] == 1:
                    xlimitBot = ESAmovie.normalizedXLimits_diff[0] * ESAmovie.vth_e
                    xlimitTop = ESAmovie.normalizedXLimits_diff[1] * ESAmovie.vth_e
                    ylimitBot = ESAmovie.normalizedYLimits_diff[0] * ESAmovie.vth_e
                    ylimitTop = ESAmovie.normalizedYLimits_diff[1] * ESAmovie.vth_e
            else:
                if wDataPlot[0] == 2:
                    ylimitBot = ESAmovie.YLimits_dist[0]
                    ylimitTop = ESAmovie.YLimits_dist[1]
                    xlimitBot = ESAmovie.XLimits_dist[0]
                    xlimitTop = ESAmovie.XLimits_dist[1]
                elif wDataPlot[0] == 1:
                    ylimitBot = ESAmovie.YLimits_diff[0]
                    ylimitTop = ESAmovie.YLimits_diff[1]
                    xlimitBot = ESAmovie.XLimits_diff[0]
                    xlimitTop = ESAmovie.XLimits_diff[1]

            X = np.linspace(xlimitBot, xlimitTop, int(N/2))
            Y = np.linspace(ylimitBot, ylimitTop, int(N/2))

            X, Y = np.meshgrid(X, Y, indexing='xy')
            all_Z_data =[[],[]]

            for i in range(2):
                prgMsg(f'Interpolating for {rocketAttrs.rocketID[i]}\n')

                if plotSpecificLocations:
                    iterateThis = range(len(scienceVar[i]))
                else:
                    iterateThis = range(len(ESAData[i]))

                for tme in tqdm(iterateThis):
                    interp = LinearNDInterpolator(list(zip(Vpars, Vperps)), scienceVar[i][tme], fill_value=fillvalue)
                    all_Z_data[i].append(interp(X, Y))

                Done(start_time)

        elif wDataPlot[0] == 0:
            all_X_data = Pitch
            all_Y_data = Energy[ESAmovie.EnergyStart:ESAmovie.EnergyEnd]
            all_Z_data = [[], []]
            for i in range(2):
                for tme in range(1, len(ESAData[i])):
                    all_Z_data[i].append(np.array([ESAData[i][tme][ptch][ESAmovie.EnergyStart:ESAmovie.EnergyEnd] for ptch in range(len(Pitch))]).transpose())

        prgMsg('Initializing Figure')

        ##########################
        # --- SETUP THE FIGURE ---
        ##########################
        fig = plt.figure()
        fig.set_figwidth(ESAmovie.figure_width)
        fig.set_figheight(ESAmovie.figure_height)

        # The overall plot shape
        gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1, 5])

        if wDataPlot[0] in [1,2]:
            # Shape of altitude and latvslong plots
            gsAlt = gs[0].subgridspec(nrows=1, ncols=3, width_ratios=[60, 60,1])
            ax1 = fig.add_subplot(gsAlt[0, 0])
            ax_lat = fig.add_subplot(gsAlt[0, 1])

            # Shape of ESA plots
            gsESA = gs[1].subgridspec(nrows=1, ncols=3, width_ratios=[60, 60, 1])
            ax2 = fig.add_subplot(gsESA[0])
            ax3 = fig.add_subplot(gsESA[1])

            # Give colorbar its own axis
            ax_colorbar = fig.add_subplot(gsESA[2])
        else:
            # Shape of altitude and latvslong plots
            gsAlt = gs[0].subgridspec(nrows=1, ncols=2, width_ratios=[1, 1])
            ax1 = fig.add_subplot(gsAlt[0, 0])
            ax_lat = fig.add_subplot(gsAlt[0, 1])

            # Shape of ESA plots
            gsESA = gs[1].subgridspec(nrows=2, ncols=2, hspace=0.3, height_ratios=[1, 1], width_ratios=[90, 1])

            # ESA plots
            ax2 = fig.add_subplot(gsESA[0, 0])
            ax3 = fig.add_subplot(gsESA[1, 0], sharex=ax2,sharey=ax2)

            # ColorBars
            ax_colorbar_high = fig.add_subplot(gsESA[0, 1])
            ax_colorbar_low = fig.add_subplot(gsESA[1, 1])

        ##############################
        # --- Initialize the plots ---
        ##############################

        title = fig.suptitle(f'{wInstr_movie.upper()}\n{EpochData_dates[0][0]}', **ESAmovie.title_style)

        # --- INITIALIZE TRAJECTORY DATA ---
        latAltMarker_high = ax1.scatter(geoLat[0][0], geoAlt[1][0], **ESAmovie.marker_style[0])
        latAltMarker_low = ax1.scatter(geoLat[0][0], geoAlt[1][0], **ESAmovie.marker_style[1])

        # Plot the altitude text
        latAltText_high = ax1.text(geoLat[0][0] - ESAmovie.textOffset_alt, geoAlt[0][0], f'{round(geoAlt[0][0], ESAmovie.altrounding)} km', ha=ESAmovie.text_alignment[1], **ESAmovie.textUTC_style_alt[0])
        latAltText_low = ax1.text(geoLat[1][0] - ESAmovie.textOffset_alt, geoAlt[1][0], f'{round(geoAlt[1][0], ESAmovie.altrounding)} km', ha=ESAmovie.text_alignment[1], **ESAmovie.textUTC_style_alt[1])

        latLongMarker_high = ax_lat.scatter(geoLat[0][0], geoLong[0][0], **ESAmovie.marker_style[0])
        latLongMarker_low = ax_lat.scatter(geoLat[1][0], geoLong[1][0], **ESAmovie.marker_style[1])

        latlongText_high = ax_lat.text(geoLat[0][0] - ESAmovie.textOffset_latlong[0][0], geoLong[0][0] - ESAmovie.textOffset_latlong[0][1], f'({round(geoLat[0][0], ESAmovie.altrounding)}' + '$^{\circ}$' + f', {round(geoLong[0][0], ESAmovie.altrounding)}' + '$^{\circ}$)', ha=ESAmovie.text_alignment[0], **ESAmovie.textUTC_style_lat[0])
        latlongText_low = ax_lat.text(geoLat[1][0] - ESAmovie.textOffset_latlong[1][0], geoLong[1][0] - ESAmovie.textOffset_latlong[1][1], f'({round(geoLat[1][0], ESAmovie.altrounding)}' + '$^{\circ}$' + f', {round(geoLong[1][0], ESAmovie.altrounding)}' + '$^{\circ}$)', ha=ESAmovie.text_alignment[1], **ESAmovie.textUTC_style_lat[1])

        for i in range(2):
            ax1.plot(geoLat[i], geoAlt[i], **ESAmovie.trajectory_style[i])  # Alttitude plot trajectory lines
            ax_lat.plot(geoLat[i], geoLong[i],**ESAmovie.trajectory_style[i]) # Lat vs Long plot trajectory lines

        # --- SET TRAJECTORY PLOTS LABEL PARAMETERS ---
        ax1.set_xlabel('Geographic Latitude', **ESAmovie.axes_label_traj_style)
        ax1.set_ylabel('Altitude [km]', **ESAmovie.axes_label_traj_style)
        ax_lat.set_xlabel('Geographic Latitude', **ESAmovie.axes_label_traj_style)
        ax_lat.set_ylabel('Geographic Longitude', **ESAmovie.axes_label_traj_style)
        trajAxes = [ax1,ax_lat]

        for i in range(2):
            trajAxes[i].tick_params(**ESAmovie.majorTicks_style)
            trajAxes[i].tick_params(**ESAmovie.minorTicks_style)

        # --- INITIALIZE ESA DATA ---
        subtitles = ['ACESII 36359', 'ACESII 36364']

        if wDataPlot[0] in [1, 2]:

            axes = [ax2, ax3]

            # Plot background color first before
            if normalizeToThermal:
                blackRadius = (2.5 * np.sin(np.radians(90)) * np.sqrt(2 * q0 * Energy[12] / m_e)) / ESAmovie.vth_e
            else:
                blackRadius = (2.5 * np.sin(np.radians(90)) * np.sqrt(2 * q0 * Energy[12] / m_e))

            center = (0, 0)

            for i in range(2):
                background_black_circle = Wedge(center, blackRadius, 340, 200, fc=(0.18995, 0.07176, 0.23217), edgecolor=None)
                axes[i].add_artist(background_black_circle)

            # --- INITIALIZE INTERPOLATED PLOTS ---
            if normalizeToThermal:
                if wDataPlot[0] == 2:
                    xlimits = ESAmovie.normalizedXLimits_dist
                    ylimits = ESAmovie.normalizedYLimits_dist

                elif wDataPlot[0] == 1:
                    xlimits = ESAmovie.normalizedXLimits_diff
                    ylimits = ESAmovie.normalizedYLimits_diff

                X = X / ESAmovie.vth_e
                Y = Y / ESAmovie.vth_e
            else:
                if wDataPlot[0] == 2:
                    xlimits = ESAmovie.XLimits_dist
                    ylimits = ESAmovie.YLimits_dist
                elif wDataPlot[0] == 1:
                    xlimits = ESAmovie.XLimits_diff
                    ylimits = ESAmovie.YLimits_diff

            xlimBot, xlimTop = xlimits[0], xlimits[1]  # values chosen based on Te = 1000eV
            ylimBot, ylimTop = ylimits[0], ylimits[1]

            if wDataPlot[0] == 2:
                # PLOT THE DISTRIBUTION DATA
                cmap1 = axes[0].pcolormesh(X, Y, all_Z_data[0][0], **ESAmovie.dist_cmap_style)
                cmap2 = axes[1].pcolormesh(X, Y, all_Z_data[1][0], **ESAmovie.dist_cmap_style)
                cmap_style = dict(cmap='turbo', vmin=ESAmovie.Vmin_dist, vmax=ESAmovie.Vmax_dist, norm='log')

            elif wDataPlot[0] == 1:
                # PLOT THE DIFFERENTIAL FLUX DATA
                cmap1 = axes[0].pcolormesh(X, Y, all_Z_data[0][0], **ESAmovie.diff_cmap_style)
                cmap2 = axes[1].pcolormesh(X, Y, all_Z_data[1][0], **ESAmovie.diff_cmap_style)
                cmap_style = dict(cmap='turbo', vmin=ESAmovie.Vmin_diff, vmax=ESAmovie.Vmax_diff, norm='log')

            # --- ADJUST PLOT LABELS ---
            for i in range(2):
                axes[i].set_xlim(xlimBot, xlimTop)
                axes[i].set_ylim(ylimBot, ylimTop)

                if normalizeToThermal:
                    axes[i].set_xlabel('$v_{\parallel}/v_{th,e}$',**ESAmovie.axes_label_ESA_style)
                    axes[i].set_ylabel('$v_{\perp}/v_{th,e}$',**ESAmovie.axes_label_ESA_style)
                else:
                    axes[i].set_xlabel('$v_{\parallel}$', **ESAmovie.axes_label_ESA_style)
                    axes[i].set_ylabel('$v_{\perp}$', **ESAmovie.axes_label_ESA_style)

                axes[i].set_title(subtitles[i], **ESAmovie.subtitle_style)
                axes[i].tick_params(**ESAmovie.majorTicks_style)
                axes[i].tick_params(**ESAmovie.minorTicks_style)

                # Get the "order of magnitude" text object
                if not normalizeToThermal:
                    text = axes[i].yaxis.get_offset_text()
                    text.set_size(ESAmovie.size_of_axis_orderOfMagnitude)
                    text = axes[i].xaxis.get_offset_text()
                    text.set_size(ESAmovie.size_of_axis_orderOfMagnitude)

            # --- INITIALIZE COLORBAR ---
            cbar = plt.colorbar(mappable=cmap1,cax=ax_colorbar, **ESAmovie.colorbar_style)

            if wDataPlot[0] == 2:
                cbar.set_label('Distribution Function [$m^{-6}s^{3}$]', **ESAmovie.cbar_label_style)
            elif wDataPlot[0] == 1:
                cbar.set_label('Differential Energy Flux [$cm^{-2}str^{-1}s^{-1} eV/eV$]', **ESAmovie.cbar_label_style)

            # --- COLORBAR LABEL PARAMS ---
            ax_colorbar.tick_params(**ESAmovie.majorTicks_colorbar_style)
            ax_colorbar.tick_params(**ESAmovie.minorTicks_colorbar_style)

            # --- INITIALIZED BLOCKING SHAPES ---
            center = (0, 0)
            theta1, theta2 = 0, 360
            radius = ESAmovie.radius_modifier * np.sin(np.radians(90)) * np.sqrt(2 * q0 * Energy[ESAmovie.EnergyEnd - 1] / m_e)

            if normalizeToThermal:
                radius = radius/ESAmovie.vth_e

            for i in range(2):

                # add sector lines
                for j in range(21):
                    sector_circle = Wedge(center, blackRadius, 10*j + 345, 10*(j+1) + 345, fc=(1,1,1,0), edgecolor='black')
                    axes[i].add_artist(sector_circle)

                semicircle = Wedge(center, radius, theta1, theta2, **ESAmovie.shapes_style)
                leftTriangle = plt.Polygon([[0.0, -0.0 * radius], [0.0, ylimBot], [-1 * ylimBot / np.tan(np.radians(15.5)), ylimBot]],**ESAmovie.shapes_style)
                rightTriangle = plt.Polygon([[0.0, -0.0 * radius], [0.0, ylimBot], [ylimBot / np.tan(np.radians(15.5)), ylimBot]],**ESAmovie.shapes_style)
                axes[i].add_artist(semicircle)
                axes[i].add_artist(leftTriangle)
                axes[i].add_artist(rightTriangle)

        else:
            vExtremes = [[1E-17, 1E-13], [1E-17, 1E-13]]
            cmap1 = ax2.pcolormesh(all_X_data, all_Y_data, all_Z_data[0][0],vmin=vExtremes[0][0], vmax=vExtremes[0][1], **ESAmovie.pcolormesh_style)
            cmap2 = ax3.pcolormesh(all_X_data, all_Y_data, all_Z_data[1][0],vmin=vExtremes[1][0], vmax=vExtremes[1][1], **ESAmovie.pcolormesh_style)

            axes = [ax2, ax3]

            # INITIALIZE COLORBARS
            axes_colorbars = [ax_colorbar_high, ax_colorbar_low]
            for i in range(2):
                cbar = plt.colorbar(**ESAmovie.colorbar_style, cax=axes_colorbars[i])
                cbar.minorticks_on()
                cbar.ax.tick_params(**ESAmovie.majorTicks_style)
                cbar.set_label('Counts', **ESAmovie.cbar_label_style)

            # --- ADJUST PLOT LABELS ---
            for i in range(2):
                axes[i].set_title(subtitles[i],**ESAmovie.subtitle_style)
                axes[i].tick_params(**ESAmovie.majorTicks_style)
                axes[i].tick_params(**ESAmovie.minorTicks_style)
                axes[i].set_ylabel('Energy [eV]', **ESAmovie.axes_label_ESA_style)
                axes[i].set_yscale('log')
                axes[i].xaxis.set_major_locator(MultipleLocator(ESAmovie.spacing_of_majorticks_ESA))
                axes[i].xaxis.set_major_formatter('{x:.0f}')

        plt.tight_layout()
        Done(start_time)
        prgMsg('Compiling the Movie')
        def animatePlot(i):

            # update Epoch title
            title.set_text(f'{wInstr_movie.upper()}\n{EpochData_dates[0][i]}')

            # update ESA DATA
            cmap1.set_array(all_Z_data[0][i])
            cmap2.set_array(all_Z_data[1][i])

            ##### UPDATE ATTITUDE DATA #####
            # --- high ---
            latAltMarker_high.set_offsets([geoLat[0][i], geoAlt[0][i]])
            latAltText_high.set_x(geoLat[0][i])
            latAltText_high.set_y(geoAlt[0][i])
            latAltText_high.set_text(f'{round(geoAlt[0][i], ESAmovie.altrounding)} km')

            latLongMarker_high.set_offsets([geoLat[0][i], geoLong[0][i]])
            latlongText_high.set_x(geoLat[0][i])
            latlongText_high.set_y(geoLong[0][i])
            latlongText_high.set_text(f'({round(geoLat[0][i], ESAmovie.altrounding)}' + '$^{\circ}$' + f', {round(geoLong[0][i], ESAmovie.altrounding)}' + '$^{\circ}$)')

            # --- low ---
            latAltMarker_low.set_offsets([geoLat[1][i],geoAlt[1][i]])
            latAltText_low.set_x(geoLat[1][i])
            latAltText_low.set_y(geoAlt[1][i])
            latAltText_low.set_text(f'{round(geoAlt[1][i], ESAmovie.altrounding)} km')

            latLongMarker_low.set_offsets([geoLat[1][i],geoLong[1][i]])
            latlongText_low.set_x(geoLat[1][i])
            latlongText_low.set_y(geoLong[1][i])
            latlongText_low.set_text(f'({round(geoLat[1][i], ESAmovie.altrounding)}' + '$^{\circ}$' + f', {round(geoLong[1][i], ESAmovie.altrounding)}' + '$^{\circ}$)')

        if plotSpecificLocations:
            locations = [i for i in range(len(specific_locations))]
        else:
            locations = [i for i in range(0, int(len(EpochData_dates[0])), frame_skips)]  # NEEDS TO BE THE HIGH FLYER LENGTH

        anim = animation.FuncAnimation(fig=fig, func=animatePlot, interval= 1000/ESAmovie.fps, frames=locations)
        anim.save(f'C:\Data\ACESII\\trajectories\\trajectory_plots\\ACESII_{wInstr_movie}.mp4')

        Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if not plotAllSkyData:
    ACESIIplotting()
else:

    rocketAttrs, data_dicts_attitude = ACESIIplotting()

    from plottingStyles import AttitudeMovie


    def sphere2cart(r, theta, phi):
        return [
            r * np.sin(np.radians(theta)) * np.cos(np.radians(phi)),
            r * np.sin(np.radians(theta)) * np.sin(np.radians(phi)),
            r * np.cos(np.radians(theta))
        ]


    # --- --- --- --- --- --- --- --- -
    # --- PREPARE DATA FOR PLOTTING ---
    # --- --- --- --- --- --- --- --- -

    # convert all Epoch data to datetimes to display them
    Epoch = data_dicts_attitude["attitude"][wFlyer - 4]["Epoch"][0]
    timeOffset = [pycdf.lib.datetime_to_tt2000(datetime.datetime(2022, 11, 20, 17, 20, 00, 107800)) - Epoch[0],
                  pycdf.lib.datetime_to_tt2000(datetime.datetime(2022, 11, 20, 17, 21, 40, 115700)) - Epoch[0]]  # number of nanoseconds from 17:20:00 each was launched
    datetimesEpoch = np.array([pycdf.lib.tt2000_to_datetime(Epoch[i] + timeOffset[wFlyer - 4]).strftime("%H:%M:%S:%f") for i in range(len(Epoch))])

    YawI =data_dicts_attitude['attitude'][wFlyer - 4]['YawI'][0]
    PitchI = data_dicts_attitude['attitude'][wFlyer - 4]['PitchI'][0]
    RollI = data_dicts_attitude['attitude'][wFlyer - 4]['RollI'][0]
    T0position = sphere2cart(1, AttitudeMovie.launcherSettings[wFlyer-4][1], 90 - AttitudeMovie.launcherSettings[wFlyer-4][2] )

    X_Az = data_dicts_attitude['attitude'][wFlyer - 4]['X_Az'][0]
    X_El = data_dicts_attitude['attitude'][wFlyer - 4]['X_El'][0]
    Y_Az = data_dicts_attitude['attitude'][wFlyer - 4]['Y_Az'][0]
    Y_El = data_dicts_attitude['attitude'][wFlyer - 4]['Y_El'][0]
    Z_Az = data_dicts_attitude['attitude'][wFlyer - 4]['Z_Az'][0]
    Z_El = data_dicts_attitude['attitude'][wFlyer - 4]['Z_El'][0]

    # Convert all the axes data into cartesian vectors
    xAxisData = np.array(
        [sphere2cart(1, 90 - X_El[i], X_Az[i]) for i in range(len(X_Az))]
    )
    yAxisData = np.array(
        [sphere2cart(1, 270 - Y_El[i], Y_Az[i]) for i in range(len(Y_Az))]
    )
    zAxisData = np.array(
        [sphere2cart(1, 90 - Z_El[i], Z_Az[i]) for i in range(len(Z_Az))]
    )

    # plot the attitude data at one time
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    ax.set_zlim([-1.5, 1.5])

    def get_arrow(i, data):
        x, y, z = 0, 0, 0
        u = data[i][0]
        v = data[i][1]
        w = data[i][2]
        return x, y, z, u, v, w


    # --- INITIALIZE ANIMATION ---

    testval = 0

    # initial direction vectors
    xQ = ax.quiver(*get_arrow(testval, xAxisData), color='red')
    yQ = ax.quiver(*get_arrow(testval, yAxisData), color='green')
    zQ = ax.quiver(*get_arrow(testval, zAxisData), color='blue')



    title = ax.set_title(f'ACESII {rocketAttrs.rocketID[wFlyer - 4]} \n {datetimesEpoch[0]}')

    ax.set_xlabel('North')
    ax.set_ylabel('West')
    ax.set_zlabel('Up')

    plt.show()


    # --- --- --- --- --- --- --
    # --- ANIMATION FUNCTION ---
    # --- --- --- --- --- --- --
    def animatePlot(i):

        # update quiver
        global xQ,yQ,zQ
        xQ.remove()
        yQ.remove()
        zQ.remove()

        xQ = ax.quiver(*get_arrow(i, xAxisData),color='red')
        yQ = ax.quiver(*get_arrow(i, yAxisData),color='green')
        zQ = ax.quiver(*get_arrow(i, zAxisData),color='blue')

        # update Epoch title
        title.set_text(f'ACESII {rocketAttrs.rocketID[wFlyer - 4]} \n {datetimesEpoch[i]}')

    # --- --- --- --- --- --- -
    # --- PERFORM ANIMATION ---
    # --- --- --- --- --- --- -
    prgMsg('Creating Movie')
    ax.view_init(40,-40)
    locations = [i for i in range(len(X_Az))]
    # locations = [i for i in range(0,10000)]
    anim = animation.FuncAnimation(fig=fig, func=animatePlot, interval=1000 / AttitudeMovie.fps, frames=locations)
    anim.save(f'D:\Data\ACESII\\trajectories\\trajectory_plots\\ACESII_{rocketAttrs.rocketID[wFlyer - 4]}_Trajectory.mp4')
    Done(start_time)
