# --- csv_to_cdf_attitude.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Turn the .cdf files of the TRICE attitude data into cdf files


# --- --- --- --- ---
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
wRocket = 4

# --- --- --- ---
# --- import ---
# --- --- --- ---
import datetime
from matplotlib import pyplot as plt, animation
from ACESII_code.class_var_func import color

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


    # --- Create Attitude Data Dict ---
    prgMsg('Collecting Attitude data')
    FolderPath = f'{ACES_data_folder}\\attitude\\'
    dataPath_attitude = {'attitude': [glob(FolderPath + rf'{fliers[0]}\\*Attitude_Solution*'), glob(FolderPath + rf'{fliers[1]}\\\\*Attitude_Solution*')]}
    data_dicts_template = {'attitude':[]}
    data_dicts_attitude = storeData(data_dicts_template,dataPath_attitude)
    Done(start_time)

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
    Epoch = data_dicts_attitude["attitude"][wRocket - 4]["Epoch"][0]
    timeOffset = [pycdf.lib.datetime_to_tt2000(datetime.datetime(2022, 11, 20, 17, 20, 00, 107800)) - Epoch[0],
                  pycdf.lib.datetime_to_tt2000(datetime.datetime(2022, 11, 20, 17, 21, 40, 115700)) - Epoch[0]]  # number of nanoseconds from 17:20:00 each was launched
    datetimesEpoch = np.array([pycdf.lib.tt2000_to_datetime(Epoch[i] + timeOffset[wRocket - 4]).strftime("%H:%M:%S:%f") for i in range(len(Epoch))])

    YawI =data_dicts_attitude['attitude'][wRocket - 4]['YawI'][0]
    PitchI = data_dicts_attitude['attitude'][wRocket - 4]['PitchI'][0]
    RollI = data_dicts_attitude['attitude'][wRocket - 4]['RollI'][0]
    T0position = sphere2cart(1, AttitudeMovie.launcherSettings[wRocket-4][1], 90 - AttitudeMovie.launcherSettings[wRocket-4][2] )

    X_Az = data_dicts_attitude['attitude'][wRocket - 4]['X_Az'][0]
    X_El = data_dicts_attitude['attitude'][wRocket - 4]['X_El'][0]
    Y_Az = data_dicts_attitude['attitude'][wRocket - 4]['Y_Az'][0]
    Y_El = data_dicts_attitude['attitude'][wRocket - 4]['Y_El'][0]
    Z_Az = data_dicts_attitude['attitude'][wRocket - 4]['Z_Az'][0]
    Z_El = data_dicts_attitude['attitude'][wRocket - 4]['Z_El'][0]

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



    title = ax.set_title(f'ACESII {rocketAttrs.rocketID[wRocket - 4]} \n {datetimesEpoch[0]}')

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
        title.set_text(f'ACESII {rocketAttrs.rocketID[wRocket - 4]} \n {datetimesEpoch[i]}')

    # --- --- --- --- --- --- -
    # --- PERFORM ANIMATION ---
    # --- --- --- --- --- --- -
    prgMsg('Creating Movie')
    ax.view_init(40,-40)
    locations = [i for i in range(len(X_Az))]
    # locations = [i for i in range(0,10000)]
    anim = animation.FuncAnimation(fig=fig, func=animatePlot, interval=1000 / AttitudeMovie.fps, frames=locations)
    anim.save(f'D:\Data\ACESII\\trajectories\\trajectory_plots\\ACESII_{rocketAttrs.rocketID[wRocket - 4]}_Trajectory.mp4')
    Done(start_time)
