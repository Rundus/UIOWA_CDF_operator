# --- csv_to_cdf_LPcal.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: convert the LPcalibration files that were performed at iowa into .cdf files
# with all the data into one file


# --- --- --- --- ---
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

justPrintFileNames = False


wRocket = 5 # 4 or 5 for high and low, respectively

inputPath_modifier = "calibration\LP_calibration\sweptCals"
outputPath_modifier = 'calibration\LP_calibration'
modifier = ''

# information to read the .csv files
wRow_names = 0 ; wRow_startData = 1  # row index corresponding to the name of the variable # which row does the data begin
getVars = [1, 3, 5, 7] # column index corresponding to the desired variables to collect from each file. These should be 2 - DAC, 4 - P1, 6 - ni_stepped, 8- ne_stepped

# plot each individual curve
plotIndividualCurveData = False

# plot the overlay of all the individual curves
plotCalCurveOverlay = False

# Size of analog window when determining the error in current derived from V_applied/R_resistor
analogWindowSize = 39 # 39 is the smallest window that has >=2 datapoints in each window bin
plotErrorsInCurrent= False

# Average the calibration curves to get a smoothed result
AverageCalCurves = False # The point of averaging is to remove the assumed 60Hz noise that was present on the output of the LP input point. We then fit the averaged dataset with a linear line
plotAverageCalCurves = True

# plot the final calibration curve with the fits
plotCalibrationCurveFitted = True

# --- Andoya Calibration toggles ---
useAndoyaData = False
useSingleSweeptData = True
slope = [0.004296482412060302, 0.004267676767676767]
intercept = [- 4.72, -4.68]

# --- Iowa fit toggles ---
# low flyer analog range [-4065,3713]
# high flyer analog range [-4061,2824]

# LINEAR - MIDDLE
fitRangeYMid = [-1, 1]
fitRangeXMid = [-4200, 2000]

# EXPO - UPPER
fitRangeXUp = [4000, 5000] # nominally [1400,3000]


# EXPO - LOWER
fitRangeXDown = [-5000, -4500] # nominally [-4070,-3000]

unitConv = 1E9


# output the data
outputData = False



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from os import remove, path
from more_itertools import sort_together
from csv import reader
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.

print(color.BOLD + color.CYAN + 'csv_to_cdf_LPcal.py' + color.END + color.END)

def csv_to_cdf_LPcal(wRocket, rocketFolderPath, justPrintFileNames):

    if wRocket == 4:
        wflyer = 0
    elif wRocket == 5:
        wflyer = 1

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'science'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}{modifier}\*.csv')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}{modifier}\\', '') for ifile in inputFiles]

    fileoutName = 'ACESII_LP_Cals_swept_Iowa.cdf'

    if justPrintFileNames:
        for i, file in enumerate(input_names):
            print('[{:.0f}] {:20s}{:5.1f} MB'.format(i, input_names[i], round(getsize(inputFiles[i]) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data' + color.END)

        # --- Create Data Dict ---
        data_dict = {}
        rocketAttrs, missionDict, data_dict_temp = ACES_mission_dicts()
        exampleAttrs = {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                        'FORMAT': 'I5', 'UNITS': None, 'LABLAXIS': None, 'VALIDMIN': None, 'VALIDMAX': None,
                        'VAR_TYPE': 'data', 'SCALETYP': 'linear'}

        #############################
        # --- STORE THE .CSV DATA ---
        #############################
        Rnames = []
        calCurveData_Analog = []
        calCurveData_Current = []
        for i, inputfile in enumerate(inputFiles):

            # get the filename identifier
            filenameOnly = inputfile.replace(f'{rocketFolderPath}{inputPath_modifier}{modifier}\\', '').replace('.csv','').replace('LP_','')
            Rnames.append(filenameOnly)
            prgMsg(f'Converting {filenameOnly}')

            # determine the resistance value of the file
            resitanceNum = filenameOnly[:filenameOnly.index("p")]
            resitanceDeci =filenameOnly[filenameOnly.index("p") + 1 : filenameOnly.index("k")]
            resistanceVal = (float(resitanceNum) + float(resitanceDeci)/1000)*1000

            # --- collect the csv data ---
            with open(inputfile) as csvFile:
                csvAllData = [row for row in reader(csvFile)]

            csvData = np.array(csvAllData[wRow_startData:], dtype='float64').transpose()

            # --- get v_stepped ---
            v_stepped = np.array(csvData[getVars[0]][:], dtype='float64')
            R = 4.997/(1.21 + 4.997)
            v_stepped = np.array([4*(R*vStep + (-5)*(1 - R)) for vStep in v_stepped]) # apply function to convert DAC output to v_stepped
            v_stepped_current = np.array(v_stepped / resistanceVal)

            # --- determine what the ADC would report. The ADC was a 12-bit, 5V converter ---
            ni_stepped = np.array(csvData[getVars[2]][:], dtype='float64')
            ne_stepped = np.array(csvData[getVars[3]][:], dtype='float64')

            # firstly, there was a diode that killed ANY voltage that was <0, so:
            ni_stepped = np.array([val if val > 0 else 0 for val in ni_stepped])
            ne_stepped = np.array([val if val > 0 else 0 for val in ne_stepped])

            ni_stepped_digital = np.array([(4096/5) * (ni_stepped[i]) for i in range(len(ni_stepped))])
            ne_stepped_digital = np.array([(4096/5) * (ne_stepped[i]) for i in range(len(ne_stepped))])
            sweep = np.array(ne_stepped_digital - ni_stepped_digital)

            if plotIndividualCurveData:
                if i in range(len(inputFiles)):
                    fig, ax = plt.subplots(2)
                    ax[0].scatter(v_stepped, sweep)
                    ax[0].set_xlabel('Voltage applied to Resistor')
                    ax[0].set_ylabel('Current seen by ADC [Analog Units]')

                    ax[1].plot([j for j in range(len(sweep))], ni_stepped_digital, [j for j in range(len(sweep))], ne_stepped_digital)
                    plt.suptitle(filenameOnly)
                    ax[1].set_xlabel('Datapoint No.')
                    ax[1].set_ylabel('Voltage')
                    plt.show()

            # --- store data ---
            calCurveData_Analog.append(sweep)
            calCurveData_Current.append(v_stepped_current)
            Done(start_time)

        ###########################
        ### GET ANDOYA CAL DATA ###
        ###########################
        if useAndoyaData:

            prgMsg('Collecting Swept Calibration Data')
            LangmuirSweptCalFiles = glob(f'{rocketFolderPath}\science\LP_calibration\{fliers[wRocket - 4]}\*_345deg_*')

            # Collect the LP data except deltaNdivN into a data dict
            data_dict_cal = {}
            for file in LangmuirSweptCalFiles:
                if 'deltaNdivN' not in file:
                    with pycdf.CDF(file) as LangmuirSweptCalFiles:
                        for key, val in LangmuirSweptCalFiles.items():
                            if key not in data_dict_cal:
                                data_dict_cal = {**data_dict_cal, **{key: [LangmuirSweptCalFiles[key][...], {key: val for key, val in LangmuirSweptCalFiles[key].attrs.items()}]}}

            Done(start_time)

            if useSingleSweeptData:
                sweptCalRanges = rocketAttrs.LPswept_cal_epoch_ranges_single_sweep[wflyer]
            else:
                sweptCalRanges = rocketAttrs.LPswept_cal_epoch_ranges[wflyer]

            # --- --- --- --- --- --- --- --- ---
            # --- COLLECT/CALC CAL RANGE DATA ---
            # --- --- --- --- --- --- --- --- ---
            calEpoch = np.array([pycdf.lib.datetime_to_tt2000(data_dict_cal['Epoch_ne_swept'][0][i]) for i in range(len(data_dict_cal['Epoch_ne_swept'][0]))])

            sweptCal_voltage = []  # voltage of the step
            sweptCal_calCurrent = []  # corresponding current based on voltage
            sweptCal_analog_current = []  # current determined from the analag calibration data
            sweptCal_Epoch = []  # epoch values of the start/end of this cal range
            for i in range(len(sweptCalRanges)):
                # find the indicies of the calibration range
                start = np.abs(calEpoch - pycdf.lib.datetime_to_tt2000(sweptCalRanges[i][0])).argmin()
                end = np.abs(calEpoch - pycdf.lib.datetime_to_tt2000(sweptCalRanges[i][1])).argmin()

                # calculate the analog current, voltage and calibration current
                sweptCal_analog_current.append(np.array(data_dict_cal['ne_swept'][0][start:end]) - np.array(data_dict_cal['ni_swept'][0][start:end]))
                sweptCal_Epoch.append(calEpoch[start:end])
                resistance = rocketAttrs.LPswept_cal_resistances[i]
                sweptCal_step = data_dict_cal['step'][0][start:end]

                def step_to_Voltage(analog_voltage):
                    return slope[wRocket - 4]*analog_voltage + intercept[wRocket - 4]

                sweptCal_voltage.append(np.array([step_to_Voltage(sweptCal_step[i]) for i in range(len(sweptCal_step))]))
                sweptCal_calCurrent.append(sweptCal_voltage[-1] / resistance)

            # --- --- --- --- --- --- --- --- ---
            # --- Remove outliers from data ---
            # --- --- --- --- --- --- --- --- ---

            # surprisingly, only the high flyer cal needed needs to be filtered for outliers
            if wRocket == 4:
                removeOutliers = True
            elif wRocket == 5:
                removeOutliers = False

            if removeOutliers:
                prgMsg('Removing Outliers')
                threshold_percent = 0.05  # if the next point is >= 50% the the distance between max/min of the dataset, it is an outlier
                repeat_process = 20  # how many times to repeat outlier search

                for process in range(repeat_process):
                    for i in range(len(sweptCalRanges)):
                        maxV = sweptCal_analog_current[i].max()
                        minV = sweptCal_analog_current[i].min()
                        threshDistance = np.abs(maxV - minV) * threshold_percent

                        # check every value and record the indices
                        remove_indices = []

                        for j in range(len(sweptCal_analog_current[i]) - 1):
                            # take an average of data to see where I am:
                            no_of_points = 10
                            AVG = sum([sweptCal_analog_current[i][j - k] for k in range(no_of_points)]) / no_of_points

                            if AVG < -3400:  # apply a stricter threshold for things close to when the current is low
                                if np.abs(sweptCal_analog_current[i][j + 1] - sweptCal_analog_current[i][j]) >= threshDistance:
                                    remove_indices.append(j + 1)
                            elif np.abs(sweptCal_analog_current[i][j + 1]) >= 4095:  # remove any points at absolute maximum
                                remove_indices.append(j + 1)
                            elif sweptCal_Epoch[i][j + 1] < 3178576.8709884263:
                                remove_indices.append(j + 1)

                        sweptCal_calCurrent[i] = np.delete(sweptCal_calCurrent[i], remove_indices, axis=0)
                        sweptCal_voltage[i] = np.delete(sweptCal_voltage[i], remove_indices, axis=0)
                        sweptCal_analog_current[i] = np.delete(sweptCal_analog_current[i], remove_indices, axis=0)
                        sweptCal_Epoch[i] = np.delete(sweptCal_Epoch[i], remove_indices, axis=0)

                Done(start_time)

            # put data into plotting datastructures
            for i in range(len(sweptCal_analog_current)):
                name = str(rocketAttrs.LPswept_cal_resistances[i]) + '_Andoya'
                Rnames.append(name)
                calCurveData_Analog.append(sweptCal_analog_current[i])
                calCurveData_Current.append(sweptCal_calCurrent[i])

        # --- --- --- --- --- --- --- --- --- --- -
        # --- FILTER THE CALIBRATION CURVE DATA ---
        # --- --- --- --- --- --- --- --- --- --- -

        if plotCalCurveOverlay:
            legend = []
            # for resistance in [Rnames[0],Rnames[1],Rnames[2],Rnames[3]]:
            for i in range(len(calCurveData_Analog)):
                stepped_current = calCurveData_Current[i]
                sweep = calCurveData_Analog[i]
                plt.scatter(sweep, stepped_current)
                plt.xlabel('Analog Value ADC Would Record [Analog]')
                plt.ylabel('Current derived from V_applied/R [A]')
                legend.append(Rnames[i])

            plt.legend(legend)
            plt.show()

        # Flatten all the data and convert it to nano-amps
        calCurveData_Analog = np.array([item for sublist in calCurveData_Analog for item in sublist])
        calCurveData_Current = np.array([item*unitConv for sublist in calCurveData_Current for item in sublist])

        # --- COLLECT LP ERROR DATA ---
        prgMsg('Collecting Error Data')
        analogBins = np.array([[num, num+analogWindowSize] for num in range(int(calCurveData_Analog.min()) - 1, int(calCurveData_Analog.max())+1, analogWindowSize)])
        stdDevErrors = []

        for bin in tqdm(analogBins):

            errors = [calCurveData_Current[i] for i in range(len(calCurveData_Analog)) if ((calCurveData_Analog[i] >= bin[0]) and (calCurveData_Analog[i] < bin[1]))]

            stdDevErrors.append(np.std(errors))

        analogBins = np.array(analogBins)
        stdDevErrors = np.array(stdDevErrors)

        if plotErrorsInCurrent:
            xData = [bin[0] for bin in analogBins]
            yData = stdDevErrors*unitConv

            fig, ax = plt.subplots()
            ax.scatter(xData, yData)
            ax.set_ylim(-0.025, 0.18)
            ax.set_title(f'$\Delta I$ vs Analog Value\n Analog Bin Size: {analogWindowSize}')
            ax.set_xlabel('Lower Bound of Analog Bin [Analog]')
            ax.set_ylabel('Std Dev Error [nA]')
            plt.show()

        Done(start_time)

        # --- Sort the data ---

        # assign the data and also sort it based on the analog values
        calCurveData_Current,calCurveData_Analog = sort_together([calCurveData_Current,calCurveData_Analog])

        # --- AVERAGE LP CAL CURVES ---
        # The point of averaging is to remove the assumed 60Hz noise that was present on the output of the LP input point. We then
        # fit the averaged dataset with a linear line

        if AverageCalCurves:
            prgMsg('Averaging Calibration Curves')

            from scipy.signal import savgol_filter
            AveragingWindow = 100
            polyOrder = 2

            # calCurveData_Current_filtered = moving_average(calCurveData_Current, AveragingWindow)
            calCurveData_Current_filtered = savgol_filter(calCurveData_Current, window_length=AveragingWindow, polyorder=polyOrder)

            if plotAverageCalCurves:
                fig, ax = plt.subplots()
                ax.scatter(calCurveData_Analog, calCurveData_Current)
                ax.plot(calCurveData_Analog, calCurveData_Current_filtered, color='red')
                ax.set_ylabel('Current [nA]')
                ax.set_xlabel('ADC Analog Value')
                ax.set_title('SG filtering Cal data')
                plt.show()

            calCurveData_Current = calCurveData_Current_filtered
            Done(start_time)

        # --- --- --- --- --- --- --- -
        # --- ORGANIZE THE CAL DATA ---
        # --- --- --- --- --- --- --- -
        yData = calCurveData_Current
        xData = calCurveData_Analog

        # break the data up into 3 fitted curves

        yData_fit_Up = []; xData_fit_Up = []
        yData_fit_Mid = []; xData_fit_Mid = []
        yData_fit_Down = []; xData_fit_Down = []

        for i in range(len(yData)):

            # Middle
            if xData[i] <= fitRangeXMid[1] and xData[i] > fitRangeXMid[0]:
                xData_fit_Mid.append(xData[i])
                yData_fit_Mid.append(yData[i])
            # upper
            elif xData[i] <= fitRangeXUp[1] and xData[i] > fitRangeXUp[0]:
                xData_fit_Up.append(xData[i])
                yData_fit_Up.append(yData[i])
            # lower
            elif xData[i] <= fitRangeXDown[1] and xData[i] > fitRangeXDown[0]:
                xData_fit_Down.append(xData[i])
                yData_fit_Down.append(yData[i])

        # convert to numpy
        yData_fit_Mid = np.array(yData_fit_Mid)
        xData_fit_Mid = np.array(xData_fit_Mid)
        yData_fit_Down = np.array(yData_fit_Down)
        xData_fit_Down = np.array(xData_fit_Down)
        yData_fit_Up = np.array(yData_fit_Up)
        xData_fit_Up = np.array(xData_fit_Up)

        #########################
        ### FIT THE CAL CURVE ###
        #########################
        def linear(x, A, B):
            y = A * x + B
            return y
        def expo(x,A,B,C,D):
            y = A*np.exp((x-B)/C) + D
            return y

        def poly(x,A,B,C,D,E):
            y = A*((B)**(C*(x-D))) - A*((B)**(-1*C*(x-D))) + E
            return y


        xData_fit_Mid = xData_fit_Mid
        yData_fit_Mid = yData_fit_Mid
        paramsMid, cov = curve_fit(poly, xData_fit_Mid, yData_fit_Mid, maxfev=1000000)


        if plotCalibrationCurveFitted:

            fig, ax = plt.subplots()

            # Middle Down
            xData = np.linspace(fitRangeXMid[0], fitRangeXMid[1], 1000)
            yData = np.array(poly(xData, *paramsMid))
            ax.plot(xData, yData, color='red')

            plt.hlines(fitRangeYMid[0], xmin=fitRangeXMid[0], xmax=fitRangeXMid[1], color='green', label='Data used for cal fit')
            plt.hlines(fitRangeYMid[1], xmin=fitRangeXMid[0], xmax=fitRangeXMid[1], color='green')
            plt.vlines(fitRangeXMid[0], ymin=fitRangeYMid[0], ymax=fitRangeYMid[1], color='green')
            plt.vlines(fitRangeXMid[1], ymin=fitRangeYMid[0], ymax=fitRangeYMid[1], color='green')

            # extra
            ax.scatter(calCurveData_Analog,calCurveData_Current)
            ax.set_ylabel('Current [nA]')
            ax.set_xlabel('Analog Val')
            ax.set_xlim(-4500, 4000)
            plt.legend([f'Func: A*B^(C*(x-D)) - A*B^(-C*(x-D)) + E\n'
                        f'A: {paramsMid[0]}\n'
                        f'B: {paramsMid[1]}\n'
                        f'C: {paramsMid[2]}\n'
                        f'D: {paramsMid[3]}\n'
                        f'E: {paramsMid[4]}\n'])
            plt.suptitle('V/R resistor Current vs Analog Circuit Response\n'
                         f'Analog Fit Range: [{fitRangeXMid[0]},{fitRangeXMid[1]}] ')
            plt.show()


        ###################################
        ### CALULCATE ERRORS IN CURRENT ###
        ###################################

        # The analog bins for the error in the current need to be converted into current themselves
        # in order to be useful
        #
        analogBins_current = [] # bins for the error in the current

        for i in range(len(analogBins)):

            # Lower
            if analogBins[i][1] <= fitRangeXMid[0]:
                analogBins_current.append([expo(bin[0],paramsDown[0],paramsDown[1],paramsDown[2],paramsDown[3]),expo(bin[1],paramsDown[0],paramsDown[1],paramsDown[2],paramsDown[3])])


        analogBins_current= np.array(analogBins_current)


        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        if outputData:
            prgMsg('Creating output file')

            # store the fit parameters
            fit_params = np.array([paramsDown,[paramsMid[0],paramsMid[1],0,0],paramsUp])
            data_dict = {**data_dict,
                         **{'fit_params': [fit_params, {'LABLAXIS': 'fit_params',
                                                              'DEPEND_0': None,
                                                              'DEPEND_1': None,
                                                              'DEPEND_2': None,
                                                              'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                              'UNITS': 'nA',
                                                              'VALIDMIN': fit_params.min(),
                                                              'VALIDMAX': fit_params.max(),
                                                              'VAR_TYPE': 'data',
                                                              'SCALETYP': 'linear'}]}}


            # store the error in the current and the respective bins
            data_dict = {**data_dict,
                         **{'errorInCurrent': [stdDevErrors, {'LABLAXIS': 'errorInCaldCurrent',
                                                                        'DEPEND_0': None,
                                                                        'DEPEND_1': None,
                                                                        'DEPEND_2': None,
                                                                        'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                                        'UNITS': 'nA',
                                                                        'VALIDMIN': stdDevErrors.min(),
                                                                        'VALIDMAX': stdDevErrors.max(),
                                                                        'VAR_TYPE': 'data',
                                                                        'SCALETYP': 'linear'}]}}
            data_dict = {**data_dict,
                         **{'errorBins': [analogBins, {'LABLAXIS': 'errorBins',
                                                              'DEPEND_0': None,
                                                              'DEPEND_1': None,
                                                              'DEPEND_2': None,
                                                              'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                              'UNITS': 'analogValue',
                                                              'VALIDMIN': analogBins.min(),
                                                              'VALIDMAX': analogBins.max(),
                                                              'VAR_TYPE': 'support_data',
                                                              'SCALETYP': 'linear'}]}}

            data_dict = {**data_dict,
                         **{'errorBins_current': [analogBins_current, {'LABLAXIS': 'errorBins_current',
                                                       'DEPEND_0': None,
                                                       'DEPEND_1': None,
                                                       'DEPEND_2': None,
                                                       'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                       'UNITS': 'nA',
                                                       'VALIDMIN': analogBins.min(),
                                                       'VALIDMAX': analogBins.max(),
                                                       'VAR_TYPE': 'support_data',
                                                       'SCALETYP': 'linear'}]}}

            # store the regions where the different fits apply
            fitRegions = np.array(fitRangeXMid)
            data_dict = {**data_dict,
                         **{'fitRegions': [fitRegions, {'LABLAXIS': 'fitRegions',
                                                                       'DEPEND_0': None,
                                                                       'DEPEND_1': None,
                                                                       'DEPEND_2': None,
                                                                       'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                       'FORMAT': 'E12.2',
                                                                       'UNITS': 'nA',
                                                                       'VALIDMIN': fitRegions.min(),
                                                                       'VALIDMAX': fitRegions.max(),
                                                                       'VAR_TYPE': 'support_data',
                                                                       'SCALETYP': 'linear'}]}}



            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'
            globalAttrsMod['Descriptor'] = "Langmuir_Probe"
            outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod,'LangmuirProbe')

            Done(start_time)

# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II Integration High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Integration Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\*.csv')) == 0:
    print(color.RED + 'There are no .csv files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        csv_to_cdf_LPcal(wRocket, rocketFolderPath, justPrintFileNames)
    else:
        csv_to_cdf_LPcal(wRocket, rocketFolderPath, justPrintFileNames)


