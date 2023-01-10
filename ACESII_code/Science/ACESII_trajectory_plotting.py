# --- csv_to_cdf_attitude.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Turn the .cdf files of the TRICE attitude data into cdf files


# --- --- --- --- ---
import time
from ACESII_code.class_var_func import Done, setupPYCDF,prgMsg
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

##### PLOT 1: Altitude vs Geodetic lat vs Geomagnetic lat ###
Plot1 = False
Plot2 = True

# Andoya_magnetic_inclination = 78.1300 deg
figure_size = (13,10)

# high/low flyer plot
plot_style11 = dict(color = 'red',linewidth = 1)
plot_style12 = dict(color = 'blue',linewidth = 1)
plot_style_MAG = dict(color = 'blue',linewidth = 1, alpha = 0.0)
plot_label_style1 = dict(size=18,labelpad=10)
plot_label_style2 = dict(size=18,labelpad=10)

#ticks
spacing_of_majorticks = 1
tick_params_major1 = dict(labelsize=15,which='major',size=15,pad=5)
tick_params_major2 = dict(labelsize=15,which='major',size=15,pad=5)
tick_params_minor1 = dict(labelsize=15,which='minor',size=10)
tick_params_minor2 = dict(labelsize=15,which='minor',size=10)

# x-mark lines
scatter_plot_Marker_style1 = dict(marker = 'x',color = 'red',s = 150)
scatter_plot_Marker_style2 = dict(marker = 'x',color = 'blue',s = 150)
scatter_text_style1 = dict(size=14, color='red')
scatter_text_style2 = dict(size=14, color='blue')
targets = [[100,150,200,250,300,350,400,450,500,550,600],[100,150,200,250,300,350,400]]# these values are "time since launch" numbers. They're chosen just so the numbers fit on the graph
scatter_text_alignment_high = ['right','right','right','right','right','left','left','left','left','left','left']
scatter_text_alignment_high_values = [-0.1,-0.1,-0.1,-0.1,-0.1,0.1,0.1,0.1,0.1,0.1,0.1] # geographic latitude shift
scatter_text_alignment_low = ['right','right','right','left','left','left','left']
scatter_text_alignment_low_values = [-0.1,-0.1,-0.1,0.1,0.1,0.1,0.1] #geomagnetic

# magnetic lines
background_style = dict(color = 'black',linewidth = 0.75, linestyle='-.')


##### PLOT 2: Lat vs Long ###
scatter_text_alignment_high_plot2 = ['left','left','left','left','left','left','left','left','left','left','left']
scatter_text_alignment_high_values_plot2 = [0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15] # geographic latitude shift
scatter_text_alignment_low_plot2 = ['right','right','right','right','right','right','right']
scatter_text_alignment_low_values_plot2 = [-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,] #geomagnetic








# --- --- --- ---
# --- import ---
# --- --- --- ---
import numpy as np
import matplotlib.pyplot as plt
from os import path
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from ACESII_code.data_paths import fliers, ACES_data_folder
from ACESII_code.data_paths import ACES_cdf_trajectories
from ACESII_code.class_var_func import color
from ACESII_code.missionAttributes import ACES_mission_dicts

setupPYCDF()
from spacepy import pycdf
from spacepy import coordinates as coord

coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field

print(color.BOLD + color.CYAN + 'csv_to_cdf.py' + color.END + color.END)



def plot_Range_vs_Altitude(inputFiles, dataFolderPath,missionDicts):

    trajectFolderPath = [dataFolderPath + rf'trajectories\\{fliers[0]}\\',dataFolderPath + rf'trajectories\\{fliers[1]}\\']

    # --- Create Data Dicts ---
    data_dicts = []
    ECEFCvals = []
    ECEFMAGCvals = []
    geodeticCvals = []
    geodeticMAGCvals = []
    geoTargetsLat = []
    geoTargetsAlt = []
    geoTargetsLong = []
    target_indices = []
    geomagnetic_lat = []


    #Create the data dictonary and store the coordinates
    for i in range(2):
        prgMsg(f'Loading data from cdf GPS files {fliers[i]}')
        data_dict = {}
        with pycdf.CDF(inputFiles[i]) as cdfDataFile:
            for key, val in cdfDataFile.items():
                data_dict = {**data_dict, **{key: [cdfDataFile[key][...], {key: val for key, val in cdfDataFile[key].attrs.items()}]}}

        ECEFcoords = np.array([data_dict['ECEFXPOS'][0],
                               data_dict['ECEFYPOS'][0],
                               data_dict['ECEFZPOS'][0]]).transpose()
        geodetic = np.array([data_dict['Alt'][0],
                             data_dict['Lat'][0],
                             data_dict['Long'][0]]).transpose() #THESE COORDINATES WONT work. Requires a reference elliposid of earth so we can add the earth radii to the altitude.

        data_dicts.append(data_dict)

        Done(start_time)

        #############################
        # --- Convert Coordinates ---
        #############################

        prgMsg(f'Converting Coordiantes for {fliers[i]}')

        # Get the times that the mission was launched in ISO datetime. Needed for geomagnetic coordinates
        ISOtime = [data_dicts[i]['Epoch'][0][j].isoformat() for j in range(len(data_dict['Epoch'][0]))]

        # --- Convert to geoMAG Coordinates ---
        cvals_GEO = coord.Coords(ECEFcoords, 'GEO', 'sph')  # define ECEF coordinates in cartesian or spherical
        cvals_GEO.ticks = Ticktock(ISOtime, 'ISO')
        cvals_GEO_MAG = cvals_GEO.convert('MAG', 'sph')

        # --- Geodetic Coordinates ---
        cvals_GDZ = coord.Coords(geodetic, 'GDZ', 'sph')
        cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
        cvals_GDZ_MAG = cvals_GDZ.convert('MAG', 'sph')

        # store data
        ECEFCvals.append(cvals_GEO)
        ECEFMAGCvals.append(cvals_GEO_MAG)

        geodeticCvals.append(cvals_GEO)
        geodeticMAGCvals.append(cvals_GDZ_MAG)

        geomagnetic_lat.append(geodeticMAGCvals[i].lati)

        # --- find the points 100s, 200s ... 600s in the data ---
        target_indices.append([ np.abs(np.array(data_dicts[i]['FlightTime'][0]) - targets[i][j] ).argmin() for j in range(len(targets[i]))])
        geoTargetsLat.append([data_dicts[i]['Lat'][0][k] for k in target_indices[i]])
        geoTargetsAlt.append([data_dicts[i]['Alt'][0][k] for k in target_indices[i]])
        geoTargetsLong.append([data_dicts[i]['Long'][0][k] for k in target_indices[i]])

        Done(start_time)


    if Plot1:

        ##################
        # --- PLOTTING ---
        ##################

        # --- Plot the Data ---

        prgMsg('Plotting Data')
        fig = plt.figure(figsize=figure_size)

        #########################
        # --- geographic data ---
        #########################
        ax1 = fig.add_subplot(111)
        ax1.plot(data_dicts[0]['Lat'][0], data_dicts[0]['Alt'][0], **plot_style11) # High Flyer
        ax1.plot(data_dicts[1]['Lat'][0], data_dicts[1]['Alt'][0], **plot_style12)  # High Flyer
        ax1.scatter(geoTargetsLat[0], geoTargetsAlt[0], **scatter_plot_Marker_style1) # High Flyer text labels
        ax1.scatter(geoTargetsLat[1], geoTargetsAlt[1], **scatter_plot_Marker_style2)  # Low Flyer text labels

        # ticks major/minor
        ax1.minorticks_on()
        ax1.set_ylabel('Altitude (km)', **plot_label_style1)
        ax1.set_xlabel('Geographic Latitude', **plot_label_style1)
        ax1.xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
        ax1.xaxis.set_major_formatter('{x:.0f}')
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(**tick_params_major1)
        ax1.tick_params(**tick_params_minor1)

        for i, txt in enumerate(targets[0]):
            # convert to UTC time
            epoch = data_dicts[0]['Epoch'][0][target_indices[0][i]]
            txtUTC = epoch.time()
            ax1.text(geoTargetsLat[0][i]+scatter_text_alignment_high_values[i],geoTargetsAlt[0][i],str(txtUTC),ha=scatter_text_alignment_high[i] ,**scatter_text_style1)
        for i, txt in enumerate(targets[1]):
            # convert to UTC time
            epoch = data_dicts[1]['Epoch'][0][target_indices[1][i]]
            txtUTC = epoch.time()


            ax1.text(geoTargetsLat[1][i] + scatter_text_alignment_low_values[i] ,geoTargetsAlt[1][i],str(txtUTC) ,ha=scatter_text_alignment_low[i] ,**scatter_text_style2)

        ##########################
        # --- geomagnetic data ---
        ##########################
        geomagnetic_lat = [geodeticMAGCvals[0].lati,geodeticMAGCvals[1].lati]
        ax2 = ax1.twiny()
        ax2.plot(geomagnetic_lat[0], data_dicts[0]['Alt'][0], **plot_style_MAG) # High Flyer
        ax2.plot(geomagnetic_lat[1], data_dicts[1]['Alt'][0], **plot_style_MAG)  # Low Flyer

        # ticks major/minor
        ax2.minorticks_on()
        ax2.xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
        ax2.xaxis.set_major_formatter('{x:.0f}')
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(**tick_params_major2)
        ax2.tick_params(**tick_params_minor2)

        ax2.set_xlabel('Geomagnetic Latitude',**plot_label_style2)

        ax3 = ax1.twinx()
        ax3.get_yaxis().set_visible(False)
        slope = -2.378842439667985558915 # corresponds to line with -78.13deg inclination
        for i in range(11):
            ax3.axline(xy1=(69+i*0.5,-1*slope), xy2=(69 + (0.5*(i+1)), 0),**background_style)

        fig.savefig(r'D:\Data\ACESII\trajectories\magnetic_latitude_trajectories.png')
        Done(start_time)

    if Plot2:

        prgMsg('Plotting Lattitude vs Longitude')

        #######################
        # --- Lat/Long data ---
        #######################
        fig = plt.figure(figsize=figure_size)

        ax1 = fig.add_subplot(111)
        ax1.grid(True)

        # lat vs long
        ax1.plot(data_dicts[0]['Lat'][0], data_dicts[0]['Long'][0], **plot_style11)  # High Flyer
        ax1.plot(data_dicts[1]['Lat'][0], data_dicts[1]['Long'][0], **plot_style12)  # Low Flyer

        # ticks
        ax1.minorticks_on()
        ax1.xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
        ax1.xaxis.set_major_formatter('{x:.0f}')
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.set_xlabel('Geographic Latitude', **plot_label_style1)
        ax1.set_ylabel('Geographic Longitude', **plot_label_style1)
        ax1.legend(['High Flyer','Low Flyer'])
        ax1.tick_params(**tick_params_major1)
        ax1.tick_params(**tick_params_minor1)

        # ax 2
        ax2 = ax1.twiny()
        ax2.plot(geomagnetic_lat[0], data_dicts[0]['Long'][0], **plot_style_MAG)  # High Flyer
        ax2.plot(geomagnetic_lat[1], data_dicts[1]['Long'][0], **plot_style_MAG)  # Low Flyer
        ax2.set_xlabel('Geomagnetic Latitude',**plot_label_style2)

        # ticks major/minor
        ax2.minorticks_on()
        ax2.xaxis.set_major_locator(MultipleLocator(spacing_of_majorticks))
        ax2.xaxis.set_major_formatter('{x:.0f}')
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(**tick_params_major2)
        ax2.tick_params(**tick_params_minor2)

        # Scatter
        ax1.scatter(geoTargetsLat[0], geoTargetsLong[0], **scatter_plot_Marker_style1) # High Flyer text labels
        ax1.scatter(geoTargetsLat[1], geoTargetsLong[1], **scatter_plot_Marker_style2)  # Low Flyer text labels


        # UTC text
        for i, txt in enumerate(targets[0]):
            epoch = data_dicts[0]['Epoch'][0][target_indices[0][i]] # convert to UTC time
            txtUTC = epoch.time()
            ax1.text(geoTargetsLat[0][i]+scatter_text_alignment_high_values_plot2[i],geoTargetsLong[0][i],str(txtUTC),ha=scatter_text_alignment_high_plot2[i] ,**scatter_text_style1)
        for i, txt in enumerate(targets[1]):
            # convert to UTC time
            epoch = data_dicts[1]['Epoch'][0][target_indices[1][i]]
            txtUTC = epoch.time()
            ax1.text(geoTargetsLat[1][i] + scatter_text_alignment_low_values_plot2[i],geoTargetsLong[1][i],str(txtUTC) ,ha=scatter_text_alignment_low_plot2[i] ,**scatter_text_style2)

        fig.savefig(r'D:\Data\ACESII\trajectories\lat_long_trajectories.png')


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
inputFiles = [ACES_cdf_trajectories[0][0],ACES_cdf_trajectories[1][0]]

if path.exists(inputFiles[0]) and path.exists(inputFiles[1]):
    missionDicts = ACES_mission_dicts()
    dataFolderPath = ACES_data_folder
    plot_Range_vs_Altitude(inputFiles, dataFolderPath,missionDicts)

else:
    print(color.RED + f'There are no .cdf files in the directory' + color.END)
