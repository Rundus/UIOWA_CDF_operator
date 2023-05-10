# --- csv_to_cdf_LPcal.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: convert the LPcalibration files that were performed at iowa into .cdf files
# with all the data into one file


# --- --- --- --- ---
import time
from ACESII_code.class_var_func import Done, setupPYCDF,prgMsg
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

justPrintFileNames = False


wRocket = 5 # 4 or 5 for high and low, respectively

inputPath_modifier = "science\LP_calibration\sweptCals"
outputPath_modifier = 'science\LP_calibration'
modifier = ''

# information to read the .csv files
wRow_names = 0 ; wRow_startData = 1  # row index corresponding to the name of the variable # which row does the data begin
getVars = [1,3,5,7] # column index corresponding to the desired variables to collect from each file. These should be 2 - DAC, 4 - P1, 6 - ni_stepped, 8- ne_stepped

# plot each individual curve
plotIndividualCurveData = False

# plot the overlay of all the individual curves
plotCalCurveOverlay = True

# plot the final calibration curve with the fits
plotCalibrationCurveFitted = True

# --- Andoya Calibration toggles ---
useAndoyaData = False
useSingleSweeptData = True
slope = [0.004296482412060302, 0.004267676767676767]
intercept = [- 4.72, -4.68]

# --- Iowa fit toggles ---
AnalogRange = [-2900, 3720] # choose the range of DAC values to accept in the final cal-curve dataset
fitRangeY = [-2900, 1500]
fitRangeX = [-0.15, 0.02]
unitConv = 1E9

# output the data
outputData = False



# --- --- --- ---
# --- import ---
# --- --- --- ---
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.optimize import curve_fit
from os import remove, path
from ACESII_code.class_var_func import outputCDFdata
from os.path import getsize
from csv import reader
from ACESII_code.data_paths import fliers, ACES_data_folder
from ACESII_code.class_var_func import color, L2_TRICE_Quick
from ACESII_code.missionAttributes import ACES_mission_dicts
from copy import deepcopy
setupPYCDF()
from spacepy import pycdf
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
                        'VAR_TYPE': 'data',
                        'SCALETYP': 'linear'}

        #############################
        # --- STORE THE .CSV DATA ---
        #############################
        Rnames = []
        calCurveData_Analog = []
        calCurveData_v_stepped = []
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
            R = 4.99/(1.25 + 4.99)
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
                    ax[1].plot([j for j in range(len(sweep))], ni_stepped_digital,[j for j in range(len(sweep))], ne_stepped_digital)
                    plt.suptitle(filenameOnly)
                    plt.show()

            # --- store data ---
            calCurveData_Analog.append(sweep)
            calCurveData_v_stepped.append(v_stepped_current)

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
                                if np.abs(sweptCal_analog_current[i][j + 1] - sweptCal_analog_current[i][
                                    j]) >= threshDistance:
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
                calCurveData_v_stepped.append(sweptCal_calCurrent[i])


        # --- --- --- --- --- --- --- --- --- --- -
        # --- FILTER THE CALIBRATION CURVE DATA ---
        # --- --- --- --- --- --- --- --- --- --- -

        if plotCalCurveOverlay:
            legend = []
            print(Rnames)
            # for resistance in [Rnames[0],Rnames[1],Rnames[2],Rnames[3]]:
            for i in range(len(calCurveData_Analog)):
                stepped_current = calCurveData_v_stepped[i]
                sweep = calCurveData_Analog[i]
                plt.scatter(stepped_current, sweep)
                legend.append(Rnames[i])

            plt.legend(legend)
            plt.show()


        calCurveData_Analog = np.array([item for sublist in calCurveData_Analog for item in sublist])
        calCurveData_v_stepped = np.array([item for sublist in calCurveData_v_stepped for item in sublist])

        # Apply a filter to remove unwanted values
        analogData_temp = []
        v_steppedData_temp = []
        for i in range(len(calCurveData_Analog)):
            if (calCurveData_Analog[i] < AnalogRange[1]) and (calCurveData_Analog[i] > AnalogRange[0]):
                analogData_temp.append(calCurveData_Analog[i])
                v_steppedData_temp.append(calCurveData_v_stepped[i])

        # assign the data and also sort it so it's ascending
        calCurveData_Analog = np.array([x for _, x in sorted(zip(v_steppedData_temp, analogData_temp))])
        calCurveData_v_stepped = np.array(sorted(v_steppedData_temp))

        # --- --- --- --- --- --- --- -
        # --- ORGANIZE THE CAL DATA ---
        # --- --- --- --- --- --- --- -

        yData = calCurveData_Analog
        xData = calCurveData_v_stepped * unitConv

        yData_fit = [];
        xData_fit = []

        for i in range(len(yData)):
            if yData[i] <= fitRangeY[1] and yData[i] > fitRangeY[0]:
                if xData[i] <= fitRangeX[1] and xData[i] > fitRangeX[0]:
                    xData_fit.append(xData[i])
                    yData_fit.append(yData[i])

        # sort the data based on the xdata
        yData_fit = np.array(yData_fit);
        xData_fit = np.array(xData_fit)


        #########################
        ### FIT THE CAL CURVE ###
        #########################
        def fitFunc(x, A, B):
            return A * x + B

        # Sort the yData based on the xData
        yData_fit = np.array([x for _, x in sorted(zip(xData_fit, yData_fit))])
        xData_fit = np.array(sorted(xData_fit))

        parameters, covariance = curve_fit(fitFunc, xData_fit, yData_fit, maxfev=100000)

        if plotCalibrationCurveFitted:
            xData_fit = np.linspace(-1, 1, 1000)
            yData_fit = np.array([fitFunc(x, parameters[0], parameters[1]) for x in xData_fit])

            plt.scatter(xData, yData)
            plt.plot(xData_fit, yData_fit, color='red')
            plt.xlabel('Current [nA]')
            plt.ylabel('Analog Val')
            plt.show()


        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        if outputData:
            prgMsg('Creating output file')

            # store the fit parameters
            data_dict = {**data_dict, **{'fit_params': [parameters, exampleAttrs]}}

            # store the error in the

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


