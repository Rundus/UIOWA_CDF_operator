# --- plottingStyles.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Location to store the settings for the plots in ACESII_trajectory_plotting



###############################################################
# --- PLOT 1: BOTH FLYERS, ALT VS LAT (w/ CONNECTING LINES) ---
###############################################################

class ALTvsLAT:

    figure_size = (13, 10)

    # high/low flyer plot
    trajectory_style = [dict(color='red', linewidth=1),dict(color='blue', linewidth=1)]
    plot_style_MAG = dict(color='blue', linewidth=1, alpha=0.0)
    axes_label_style = [dict(size=18, labelpad=10),dict(size=18, labelpad=10)]

    # ticks
    spacing_of_majorticks = 1
    majorTicks_style = [dict(labelsize=15, which='major', size=15, pad=5),dict(labelsize=15, which='major', size=15, pad=5)]
    minorTicks_style = [dict(labelsize=15, which='minor', size=10),dict(labelsize=15, which='minor', size=10)]

    # x-mark lines
    marker_style = [dict(marker='x', color='red', s=150), dict(marker='x', color='blue', s=150)]
    textUTC_style = [dict(size=14, color='red'),dict(size=14, color='blue')]

    # text alignement
    scatter_text_alignment = [['right', 'right', 'right', 'right', 'right', 'left', 'left', 'left', 'left','left', 'left'],['right', 'right', 'right', 'left', 'left', 'left', 'left']]
    scatter_text_alignment_offsets = [[-0.1, -0.1, -0.1, -0.1, -0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],[-0.1, -0.1, -0.1, 0.1, 0.1, 0.1, 0.1]]

    # magnetic lines
    mag_background_style = dict(color='black', linewidth=2, linestyle='-.', alpha=0.2)

    # Connecting lines
    connecting_lines_style = dict(color='green', linestyle="--", alpha=0.7)



################################################################
# --- PLOT 2: BOTH FLYERS, LAT VS LONG (w/ CONNECTING LINES) ---
################################################################
class LATvsLONG:

    figure_size = (13, 10)

    # Trajectory Plot
    trajectory_style = [dict(color='red', linewidth=1),dict(color='blue', linewidth=1)]
    axes_label_style = [dict(size=18, labelpad=10),dict(size=18, labelpad=10)]

    # invisible geomag plots
    plot_style_MAG = dict(color='blue', linewidth=1, alpha=0.0)

    # ticks
    spacing_of_majorticks = 1
    majorTicks_style = [dict(labelsize=15, which='major', size=15, pad=5),
                        dict(labelsize=15, which='major', size=15, pad=5)]
    minorTicks_style = [dict(labelsize=15, which='minor', size=10), dict(labelsize=15, which='minor', size=10)]

    # x-mark lines
    marker_style = [dict(marker='x', color='red', s=150), dict(marker='x', color='blue', s=150)]
    textUTC_style = [dict(size=14, color='red'), dict(size=14, color='blue')]

    # text alignement
    scatter_text_alignment = [['left', 'left', 'left', 'left', 'left', 'left', 'left', 'left', 'left', 'left', 'left'],
                              ['right', 'right', 'right', 'right', 'right', 'right', 'right']]
    scatter_text_alignment_offsets = [[0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15],
                                      [-0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, ]]


###########################################################################
# --- PLOT 3: BOTH FLYERS, ILAT,LAT VS ILONG,LONG (w/ CONNECTING LINES) ---
###########################################################################
class ILATvsILONG:

    targetProjectionAltitude = 100
    refLong = 16.020833
    refLat = 69.294167

    figureSize = (9, 9)

    trajectory_style = [dict(color='red', linewidth=1,label='High Flyer Trajectory'), dict(color='blue', linewidth=1,label='Low Flyer Trajectory')]
    intersections_style = [dict(color='salmon', linewidth=1,label='High Flyer - Geograhpic Projection'),
                           dict(color='deepskyblue', linewidth=1,label='Low Flyer - Geographic Projection')]
    axes_label_style = dict(size=18, labelpad=10)

    # invisible geomag plots
    plot_style_MAG = dict(color='blue', linewidth=1, alpha=0.0)

    # Grid
    gridMajor_style = dict(which='Major', alpha=0.5)
    gridMinor_style = dict(which='Minor', alpha=0.2)

    # ticks
    spacing_of_majorticks = 1
    spacing_of_majorticks_km = 50
    majorTicks_style = [dict(labelsize=15, which='major', size=15, pad=5), dict(labelsize=15, which='major', size=15, pad=5)]
    minorTicks_style = [dict(labelsize=15, which='minor', size=10), dict(labelsize=15, which='minor', size=10)]

    # x-mark lines
    marker_style = [dict(marker='x', color='red', s=150), dict(marker='x', color='blue', s=150)]
    marker_intersect_style = dict(marker='x', color='salmon', s=150)
    ha = ['left', 'right']
    textUTC_style = [dict(size=10, color='red',ha = ha[0]), dict(size=10, color='blue',ha=ha[1])]
    textUTC_intersect_style = dict(size=10, color='salmon', ha=ha[0])
    adjusts_km = [[1, 0], [-1 * 2, -1 * 8]]  # km adjustments to the timestamp labels
    adjusts = [[0.05, 0.05], [-0.05, -0.05]] # deg adjustments to the timestamp labels

    # Connecting lines
    connecting_lines_style = dict(color='green', linestyle="--", alpha=0.7)




############################################
# --- PLOT 4: ALT vs EPOCH (ESA OVERLAY) ---
############################################

class ESAoverlay:

    figure_height = 15
    figure_width = 23

    # "Zoom" into data
    lowFlyer_adjust = 1000

    # given in geodetic Latitude
    start = [71,70.5]
    end = [73.4,73.5]

    # pcolormesh
    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.1, 0.5, 0.5),
                     (0.2, 0.0, 0.0),
                     (0.4, 0.2, 0.2),
                     (0.6, 0.0, 0.0),
                     (0.8, 1.0, 1.0),
                     (1.0, 1.0, 1.0)),'green': ((0.0, 0.0, 0.0),
                       (0.1, 0.0, 0.0),
                       (0.2, 0.0, 0.0),
                       (0.4, 1.0, 1.0),
                       (0.6, 1.0, 1.0),
                       (0.8, 1.0, 1.0),
                       (1.0, 0.0, 0.0)),'blue': ((0.0, 0.0, 0.0),
                      (0.1, 0.5, 0.5),
                      (0.2, 1.0, 1.0),
                      (0.4, 1.0, 1.0),
                      (0.6, 0.0, 0.0),
                      (0.8, 0.0, 0.0),
                      (1.0, 0.0, 0.0))}
    from matplotlib import colors
    my_cmap = colors.LinearSegmentedColormap('my_colormap', cdict, 256)
    pcolormesh_style = dict(vmin=2,vmax=15,cmap=my_cmap,alpha=0.8)
    colorbar_style = dict(pad=0.04)
    cbar_label_style = dict(size=10)

    # high/low flyer plot
    trajectory_style = [dict(color='red', linewidth=1), dict(color='blue', linewidth=1)]
    plot_style_MAG = dict(color='blue', linewidth=1, alpha=0.0)
    axes_label_style = dict(size=13, labelpad=5)
    yaxis_lim = dict(bottom = 20, top = 1300)
    title_style = dict(size=10)
    margins_style = dict(tight=True)


    # ticks
    spacing_of_majorticks = 1
    majorTicks_style = dict(labelsize=12, which='major', size=10, pad=5)
    minorTicks_style = dict(labelsize=8, which='minor', size=8)


class ESAmovie:

    # fps of the video
    fps = 20

    # figure size
    figure_height = 15
    figure_width = 25

    ################
    # TRAJECT PLOT #
    ################

    # high/low flyer plot
    trajectory_style = [dict(color='red', linewidth=2), dict(color='blue', linewidth=2)]
    axes_label_traj_style = dict(size=20, labelpad=5)
    title_style = dict(size=25)

    # ticks
    spacing_of_majorticks_Traj = 0.5
    spacing_of_majorticks_ESA = 10
    majorTicks_style = dict(labelsize=20, which='major', size=15, pad=2)
    minorTicks_style = dict(labelsize=15, which='minor', size=10)

    # Altitude Text
    textOffset = 0.15
    text_alignment = 'right'
    textUTC_style = [dict(size=22, color='red'), dict(size=22, color='blue')]
    altrounding = 1 # how many decimals to round to

    # x-mark lines
    marker_style = [dict(marker='o', color='red', s=175), dict(marker='o', color='blue', s=175)]



    ##################
    # ESA DATA PLOTS #
    ##################

    # if plotDistFunc:
    vmin = [1E-18, 1E-18]
    vmax = [1E-12, 1E-12]
    # else:
    #     vmin = [ESAData[0].min(), ESAData[1].min()]
    #     vmax = [ESAData[0].max(), ESAData[1].max()]

    # pcolormesh
    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.1, 0.5, 0.5),
                     (0.2, 0.0, 0.0),
                     (0.4, 0.2, 0.2),
                     (0.6, 0.0, 0.0),
                     (0.8, 1.0, 1.0),
                     (1.0, 1.0, 1.0)), 'green': ((0.0, 0.0, 0.0),
                                                 (0.1, 0.0, 0.0),
                                                 (0.2, 0.0, 0.0),
                                                 (0.4, 1.0, 1.0),
                                                 (0.6, 1.0, 1.0),
                                                 (0.8, 1.0, 1.0),
                                                 (1.0, 0.0, 0.0)), 'blue': ((0.0, 0.0, 0.0),
                                                                            (0.1, 0.5, 0.5),
                                                                            (0.2, 1.0, 1.0),
                                                                            (0.4, 1.0, 1.0),
                                                                            (0.6, 0.0, 0.0),
                                                                            (0.8, 0.0, 0.0),
                                                                            (1.0, 0.0, 0.0))}
    from matplotlib import colors
    my_cmap = colors.LinearSegmentedColormap('my_colormap', cdict, 256)
    pcolormesh_style = dict(cmap=my_cmap, alpha=1)
    colorbar_style = dict(pad=-4)
    cbar_label_style = dict(size=25, labelpad=50,rotation=270)

    # axes label spectrogram
    axes_label_ESA_style = dict(size=25, labelpad=5)

    # "Zoom" into energy range
    EnergyStart = 12
    EnergyEnd = 42





