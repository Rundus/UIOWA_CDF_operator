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

wRocket = 4 # 4 or 5 for high and low, respectively

inputPath_modifier = "science\LP_calibration\sweptCals"
outputPath_modifier = 'science\LP_calibration'
modifier = ''

wRow_names = 0 ; wRow_startData = 1  # row index corresponding to the name of the variable # which row does the data begin
getVars = [1,3,5,7] # column index corresponding to the desired variables to collect from each file. These should be 2 - DAC, 4 - P1, 6 - ni_stepped, 8- ne_stepped

# plot each individual curve
plotIndividualCurveData = False

# plot the overlay of all the individual curves
plotCalCurveOverlay = False

# toggles for the final dataset curve
AnalogRange = [-4294,3720] # choose the range of DAC values to accept in the final cal-curve dataset


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
from ACESII_code.data_paths import fliers, ACES_data_folder,ACES_csv_trajectories,TRICE_data_folder,ACES_L0_files,TRICE_L0_files
from ACESII_code.class_var_func import color, L2_TRICE_Quick
from ACESII_code.missionAttributes import ACES_mission_dicts,TRICE_mission_dicts
from copy import deepcopy
setupPYCDF()
from spacepy import pycdf
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field

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
                    fig,ax = plt.subplots(2)
                    ax[0].scatter(v_stepped, sweep)
                    ax[1].plot([j for j in range(len(sweep))], ni_stepped_digital,[j for j in range(len(sweep))], ne_stepped_digital)
                    plt.suptitle(filenameOnly)
                    plt.show()

            # --- store data ---
            calCurveData_Analog.append(sweep)
            calCurveData_v_stepped.append(v_stepped_current)

            Done(start_time)

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

        # --- --- --- --- --- --- --- --- --- --- -
        # --- FILTER THE CALIBRATION CURVE DATA ---
        # --- --- --- --- --- --- --- --- --- --- -

        # flatten the data
        calCurveData_Analog = np.array([item for sublist in calCurveData_Analog for item in sublist])
        calCurveData_v_stepped = np.array([item for sublist in calCurveData_v_stepped for item in sublist])

        # Apply a filter to remove unwanted values
        analogData_temp = []
        v_steppedData_temp = []
        for i in range(len(calCurveData_Analog)):
            if (calCurveData_Analog[i] < AnalogRange[1]) and (calCurveData_Analog[i] > AnalogRange[0]) :
                analogData_temp.append(calCurveData_Analog[i])
                v_steppedData_temp.append(calCurveData_v_stepped[i])

        # assign the data and also sort it so it's ascending
        calCurveData_Analog = np.array([x for _, x in sorted(zip(v_steppedData_temp,analogData_temp))])
        calCurveData_v_stepped = np.array(sorted(v_steppedData_temp))


        # --- --- --- --- --- --- ---
        # --- FIT THE CAL CURVE ---
        # --- --- --- --- --- --- ---
        unit_convert = 1E9 # for easier fitting, switch to nA

        yData = calCurveData_Analog
        xData = calCurveData_v_stepped * unit_convert

        # break data above and below a certain point
        adjust = 10
        breakpoint = 100
        breakIndex = np.abs(yData - breakpoint).argmin()
        rangeLow = np.abs(yData - AnalogRange[0]).argmin() ; rangeHigh = np.abs(yData - AnalogRange[1]).argmin()
        yData_D = yData[rangeLow:breakIndex]; xData_D = xData[rangeLow:breakIndex]
        yData_U = yData[breakIndex:rangeHigh]; xData_U = xData[breakIndex:rangeHigh]

        # UP
        def upFunc(x, A, B, C):
            return A*np.log(1*(x + B)) + C

        parameters, covariance = curve_fit(upFunc, xData_U, yData_U, maxfev=100000)

        xData_U_fit = np.linspace(xData_D.max(), xData[rangeLow:(breakIndex+adjust)],100000)
        yData_U_fit = np.array([upFunc(x, parameters[0], parameters[1], parameters[2]) for x in xData_U_fit])

        # DOWN
        def downFunc(x, A, B, C):
            return -1*A*np.log(-1*(x - B)) + C

        parameters, covariance = curve_fit(downFunc, xData_D, yData_D, maxfev=10000)

        xData_D_fit = np.linspace(xData_D.min(), xData[breakIndex:(rangeHigh - adjust)], 100000)
        yData_D_fit = np.array([downFunc(x, parameters[0], parameters[1], parameters[2]) for x in xData_D_fit])

        # plot the data
        fig, ax = plt.subplots(2)
        ax[0].scatter(xData_U, yData_U)
        ax[1].scatter(xData_D, yData_D)
        ax[0].plot(xData_U_fit, yData_U_fit, color='red')
        ax[1].plot(xData_D_fit, yData_D_fit,color='red')

        # plot the base data
        plt.suptitle('I_analog vs I_stepped')
        plt.show()

        # overplot the two fits to check overlap
        fig, ax = plt.subplots(1)
        ax.scatter(xData, yData)
        ax.plot(xData_U_fit, yData_U_fit,color='purple')
        ax.plot(xData_D_fit, yData_D_fit, color='green')
        plt.show()

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        if outputData:
            prgMsg('Creating output file')

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


