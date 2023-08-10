# --- L2_Langmuir_to_SciLangmuir.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Takes in the L2 Langmuir Data and performs the Chi Square
# analysis to get the Temperature and Density from the characteristic curves


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
justPrintFileNames = False

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'L2'  # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\Langmuir'  # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
errorPath_modifier = 'calibration\LP_calibration'

# --- Fitting Toggles ---
unitConv = 1E9  # converts from A to nA
wSweeps = []  # [] --> all sweeps, [#1,#2,...] specific sweeps

IgnoreTheseSweeps = [[22, 23], []] # some particular sweeps in the data are bad, ignore them

###################
### FIXED PROBE ###
###################
calculateFixedni = True

###################
### SWEPT PROBE ###
###################
# bad data: High Flyer 36
# auroral case: 260
# calm case:
# Exponential case: 40, 100, 70,90, 150
# Linear Cases: 130


# plot the sweeps in wSweeps
plotSpecificSweeps = True

# --- TRADITIONAL FIT TOGGLES ---
traditionalFit = False

showTeFittingMethod = False # shows how Te is derived. Using the lowest Chisquare
fitsStartHere = 0.7  # in volts. Start searching for the linear region here
fitsEndHere = 1.8  # in volts. Finish searching for the linear region here

showVspFittingMethod = False


alternateData = False
wSet = 1 # nominally == 0 or 1, determines keeping every odd sweep or even sweep


#^^^ I noticed a pattern in the data: every other sweep follows its own pattern,
# almost like two datasets are being sampled. Splitting the data in 2 may reveal something

# --- GRID SEARCH Vsp FIT TOGGLES ---
gridSearchFit = False # the toggle parameters are in LP_gridSearch_toggles
showGridSearchFitting = False # shows the best gridsearch fit

# --- OutputData ---
outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---

from LP_gridSearch_toggles import searchStartsHere, transFitParameters,satFitParameters,transGuess,transBounds,a0_range,a1_range,a2_range,satNumPoints


##############################
# --- FITTED CHI FUNCTIONS ---
##############################
rocketAttrs, b, c = ACES_mission_dicts()
e_coefficient = (((q0 ** 3) * (rocketAttrs.LP_probe_areas[0][0] ** (2))) / (2 * np.pi * m_e)) ** (1 / 2)

def transitionFunc(x, a0, a1, a2):
    y = (e_coefficient) * a0 * np.exp((x - a1) / a2)
    return y
def saturationFunc(x, a0, a1, a2):
    y = a0*(e_coefficient)*(x/a2 + (1- a1/a2))
    return y


#######################
# --- MAIN FUNCTION ---
#######################
def main(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer,wSweeps):
    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'LangmuirData'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*langmuir*')
    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*Temp&Density*')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace(inputPath_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace(outputPath_modifier.lower() + '_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\',
                                              '')
    fileoutName = f'ACESII_{rocketID}_langmuir_Temp&Density.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made L2: {:3s} '.format(i, input_names_searchable[i], round(os.path.getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(os.path.getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the L2 file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict = {}
        data_dict = loadDictFromFile(inputFiles[wFile],data_dict)
        data_dict['Epoch_swept_Current'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_swept_Current'][0][i]) for i in (range(len(data_dict['Epoch_swept_Current'][0])))])
        Done(start_time)

        # --- get error data from calibration ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict_errors = {}
        errorsPath = glob(f'{ACES_data_folder}{errorPath_modifier}\{fliers[wflyer]}\*Iowa*')
        data_dict_errors = loadDictFromFile(errorsPath[0],data_dict_errors)
        Done(start_time)


        if calculateFixedni:
            ##################################################################
            # --- Calculate the plasma density from Ion Saturation Current ---
            ##################################################################

            ni = data_dict['fixed_Current'][0]
            unit_conversion = 10 ** 9  # 1 for Amps 10**9 for nano amps, etc

            # determining n_i from Ion saturation
            Ti_assumed = 0.1  # assuming an ion temperature of 0.1eV

            # using the fixed LP data (now calibrated), determine n_i from the basic ion saturation current equation
            vth_i = np.sqrt((8 * Ti_assumed * q0) / (IonMasses[0] * np.pi))
            ni_density = np.array([(4 * (current)) / ((cm_to_m ** 3) * unit_conversion * q0 * vth_i * rocketAttrs.LP_probe_areas[0][0]) for current in ni])

            # quality assurance check
            for i, val in enumerate(ni_density):
                if np.abs(val) > 1E20:
                    ni_density[i] = rocketAttrs.epoch_fillVal

            data_dict = {**data_dict, **{'fixed_ni_density': [ni_density, {'LABLAXIS': 'plasma density',
                                                                           'DEPEND_0': 'fixed_Epoch',
                                                                           'DEPEND_1': None,
                                                                           'DEPEND_2': None,
                                                                           'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                           'FORMAT': 'E12.2',
                                                                           'UNITS': 'cm!U-3!',
                                                                           'VALIDMIN': ni_density.min(),
                                                                           'VALIDMAX': ni_density.max(),
                                                                           'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
            import scipy
            ni_density_filtered = scipy.signal.savgol_filter(ni_density, window_length=10000, polyorder=2, mode='nearest')
            data_dict = {**data_dict,
                         **{'fixed_ni_density_smoothed': [ni_density_filtered, {'LABLAXIS': 'plasma density_smoothed',
                                                                                'DEPEND_0': 'fixed_Epoch',
                                                                                'DEPEND_1': None,
                                                                                'DEPEND_2': None,
                                                                                'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                                'FORMAT': 'E12.2',
                                                                                'UNITS': ''
                                                                                         '',
                                                                                'VALIDMIN': ni_density.min(),
                                                                                'VALIDMAX': ni_density.max(),
                                                                                'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}

            Done(start_time)


        ######################################
        # --- FIND THE SWEPT CURVE INDICES ---
        ######################################

        prgMsg('Finding the swept curve indicies')

        # Find the start and end point for each sweep
        indvEpochThresh = 250000000  # Value, in tt2000, that determines the time diff needed between epoch points to identify particular sweeps
        counter = 0
        sweepIndices = []

        for i in (range(len(data_dict['Epoch_swept_Current'][0]) - 1)):

            if np.abs((data_dict['Epoch_swept_Current'][0][i + 1] / 100000 - data_dict['Epoch_swept_Current'][0][i] / 100000)) >= indvEpochThresh / 100000:

                if counter == 0:
                    sweepIndices.append([0, i])
                    counter += 1
                else:
                    sweepIndices.append([sweepIndices[-1][1] + 1, i])

        # include the last sweep
        sweepIndices.append([sweepIndices[-1][1] + 1, len(data_dict['Epoch_swept_Current'][0])])

        # Add a QA set where we remove anyplace where sweepIndicies[i][0] == sweepIndicies[i][1]
        removeThese = []
        for i in range(len(sweepIndices)):
            if sweepIndices[i][0] == sweepIndices[i][1]:
                removeThese.append(i)

        for thing in removeThese:
            del sweepIndices[thing]

        # Find the index corresponding to the break in upleg/downleg voltage for the sweeps
        breakIndices = []
        qualityCounter = 0
        for i in range(len(sweepIndices)):
            start = sweepIndices[i][0]
            end = sweepIndices[i][1]
            sweeps = np.array(data_dict['swept_Current'][0][start:end])
            sweepMAX, sweepMIN = sweeps.max(), sweeps.min()
            threshold = (sweepMAX - sweepMIN) * 0.75

            for j in range(len(sweeps) - 1):
                if np.abs(sweeps[j + 1] - sweeps[j]) >= threshold:
                    breakIndices.append(j + start)

                    if start < (j + start) < end:
                        qualityCounter += 1

        Done(start_time)

        print(len(sweepIndices), qualityCounter)

        sweepIndices, breakIndices = np.array(sweepIndices), np.array(breakIndices)

        # --- store the sweep data into a single data variables filled with individual sweeps---
        prgMsg('Reorganizing Data')
        sweepsCurrent, sweepsVoltage, sweepsCurrent_epoch = [], [], []

        for sweep in range(len(sweepIndices)):
            start = sweepIndices[sweep][0]
            end = sweepIndices[sweep][1]
            breakPoint = breakIndices[sweep]

            # Get the data and sort it
            xData1 = np.array(data_dict['swept_Voltage'][0][start:breakPoint])
            yData1 = np.array(data_dict['swept_Current'][0][start:breakPoint])
            yData1 = np.array([x for _, x in sorted(zip(xData1, yData1))])
            xData1 = np.array(sorted(xData1))

            xData2 = data_dict['swept_Voltage'][0][breakPoint:end]
            yData2 = np.array(data_dict['swept_Current'][0][breakPoint:end])

            yData2 = np.array([x for _, x in sorted(zip(xData2, yData2))])
            xData2 = np.array(sorted(xData2))

            # append data to lists
            sweepsVoltage.append(xData1)
            sweepsVoltage.append(xData2)
            sweepsCurrent_epoch.append(data_dict['Epoch_swept_Current'][0][start])
            sweepsCurrent.append(yData1)
            sweepsCurrent.append(yData2)
            sweepsCurrent_epoch.append(data_dict['Epoch_swept_Current'][0][breakPoint])

        sweepsVoltage = np.array(sweepsVoltage, dtype='object')
        sweepsCurrent = np.array(sweepsCurrent, dtype='object')
        sweepsCurrent_epoch = np.array(sweepsCurrent_epoch, dtype='object')
        Done(start_time)

        ###############################################
        # --- FIT THE DATA TO GET PLASMA PARAMETERS ---
        ###############################################


        # If wSweeps = [], do all sweeps otherwise do the specified ones
        if wSweeps == []:
            wSweeps = [i for i in range(len(sweepsVoltage))]
        else:
            wSweeps = wSweeps


        if plotSpecificSweeps:

            # Make a plot of the sweeps in sweeps
            for sweep in wSweeps:
                xData, yData = sweepsVoltage[sweep], np.array(sweepsCurrent[sweep])
                EpochStamp = pycdf.lib.tt2000_to_datetime(sweepsCurrent_epoch[sweep])

                fig, ax = plt.subplots()
                ax.scatter(xData, yData)
                ax.set_xlabel('Probe Voltage')
                ax.set_ylabel('Current [nA]')
                ax.set_title(f'ACESII {rocketAttrs.rocketID[wRocket - 4]} \n'
                             f'{EpochStamp.strftime("%H:%M:%S")}\n'
                             f'Sweep No: {sweep}')
                plt.show()


        # use a grid-search technique and chiSquare fitting to get the plasma parameters
        if gridSearchFit:

            ###########################
            # --- APPLY GRID SEARCH ---
            ###########################

            prgMsg('Performing Grid Search')
            print('\n')

            # data that will be outputed
            sweeps_Te = []
            sweeps_Vsp = []
            sweeps_n0 = []
            sweeps_Epoch = []

            # analyze the specific sweeps in wSweeps or do all of them
            theseSweeps = [i for i in range(len(sweepsCurrent))] if wSweeps == [] else wSweeps

            # load in the error data that will be used later
            errorBins_current = [data_dict_errors['errorBins_current'][0][i][0] for i in range(len(data_dict_errors['errorBins_current'][0]))]  # get just the left boundary of the bins so we can use np.digitize
            errorInCurrent = np.array(data_dict_errors['errorInCurrent'][0])

            # loop through all the sweeps
            for sweep in theseSweeps:

                # [1] DETERMINE ESTIMATE OF ION SATURATION CURRENT
                # description: find the most negative value of current in the data and upshift data by it
                xData, yData = np.array(sweepsVoltage[sweep]), np.array(sweepsCurrent[sweep]/unitConv)
                EpochStamp = pycdf.lib.tt2000_to_datetime(sweepsCurrent_epoch[sweep])

                print('xData',xData)
                print('yData',yData)

                Isat_Est = yData.min()

                yData_upshift = np.array([y - Isat_Est for y in yData])

                # [2] GRID SEARCH OVER FIT PARAMETERS TO FIND BEST CHI
                # description: The parameters Vsp, nTe^1/2, Te show up in both the transition and saturation region
                # perform a grid search over these variables by choosing a range of values of each of them,
                # then loop over all possiblities until a best ChiSquare is found. Store these values in "fitresults"
                # Furthermore, the electron saturation region is fit using 3,4,5 points in each curve to determine
                # the best chiSquare fit

                fitresults = [1E8,1,1,1,1,1,1,1,1] # format: [Chitotal, chiTrans, XTrans, YTrans, ChiSat, XSat, YSat, [nTe, Vsp, Te],numPointsIe0]]


                # loop over grid search paramters: Vsp, Te, nTe^1/2

                gridSearchRanges = [a1_range,a2_range,a0_range,satNumPoints]

                for Vsp,Te,nTe,numPoints in tqdm(itertools.product(*gridSearchRanges)):

                    # [3] FIT TRANSITION REGION, CALC CHI SQUARE
                    # description: using the exponential fit function, fit the transition region
                    # Everything is done in Amps, not nanoamps to derive the parameters easier

                    # get the transition region data
                    breakIndex_Vsp = np.abs(xData - Vsp).argmin() # maximum x-val for transition region
                    breakIndex_Vmin = np.abs(xData - searchStartsHere).argmin() # minimum x-val for the transition region
                    xDataTrans = np.array(xData[breakIndex_Vmin:breakIndex_Vsp])
                    yDataTrans = np.array(yData_upshift[breakIndex_Vmin:breakIndex_Vsp])

                    # convert Ii0_Est to amps
                    Isat_Est = Isat_Est

                    # Assign the errors
                    try:
                        wBins = np.digitize(yDataTrans, errorBins_current)
                        wErrorsTrans = np.array([errorInCurrent[y] for y in wBins])/unitConv
                    except:
                        wErrorsTrans = [10E-20 for y in wBins]

                    # calc ChiSquare
                    nu = 1/(len(yDataTrans) - 3)
                    ChiSquareTrans = (1/nu) * sum([(yDataTrans[i] - (transitionFunc(xDataTrans[i], nTe, Vsp, Te) + Isat_Est)  )**2 / (wErrorsTrans[i]**2) for i in range(len(yDataTrans))])

                    # fit the saturation curve
                    xDataSat = np.array(xData[-1*numPoints:])
                    yDataSat = np.array(yData[-1*numPoints:])

                    # Assign the errors
                    wBins = np.digitize(yDataSat, errorBins_current)
                    wErrorsSat = np.array([errorInCurrent[y] for y in wBins])

                    # calculate ChiSquare
                    nu = 1/(len(yDataSat) - 2)
                    ChiSquareSat = (1/nu) * sum([((yDataSat[i] - (saturationFunc(xDataSat[i],nTe,Vsp,Te) - np.abs(Isat_Est)))**2)/(wErrorsSat[i]**2) for i in range(len(xDataSat))])

                    # store the ChiSquares if they're good
                    ChiTotal = ChiSquareTrans + ChiSquareSat

                    if ChiTotal <fitresults[0]:
                        # [Chitotal, chiTrans, XTrans, YTrans, ChiSat, XSat, YSat, [nTe, Vsp, Te],numPointsIe0]]
                        fitresults[0] = ChiTotal
                        fitresults[1] = ChiSquareTrans
                        fitresults[2] = xDataTrans
                        fitresults[3] = yDataTrans
                        fitresults[4] = ChiSquareSat
                        fitresults[5] = xDataSat
                        fitresults[6] = yDataSat
                        fitresults[7] = [nTe, Vsp, Te]
                        fitresults[8] = numPoints

                if showGridSearchFitting:
                    xDataTrans = fitresults[2]
                    nTe = fitresults[7][0]
                    Vsp = fitresults[7][1]
                    Te = fitresults[7][2]
                    xDataSat = fitresults[5]
                    yDataSat = fitresults[6]
                    numPoints = fitresults[8]
                    ChiSquareTrans = fitresults[1]
                    ChiSquareSat = fitresults[4]

                    fig,ax = plt.subplots()

                    # plot data
                    ax.scatter(xData, yData)
                    ax.set_ylabel('Current [A]')
                    ax.set_xlabel('Probe Voltage [V]')

                    # plot the transition fit
                    xTransLine = np.linspace(min(xDataTrans), max(xDataTrans))
                    yTransLine = np.array([transitionFunc(x, nTe, Vsp, Te) - np.abs(Isat_Est) for x in xTransLine])
                    ax.plot(xTransLine, yTransLine, color='red')

                    # plot saturation fit
                    xSatLine = np.linspace(max(xDataTrans), max(xDataSat))
                    ySatLine = np.array([saturationFunc(x, nTe, Vsp, Te) for x in xSatLine])
                    ax.scatter(xDataSat, yDataSat, color='purple')
                    ax.plot(xSatLine, ySatLine, color='black')

                    plt.suptitle(f'ACESII {rocketAttrs.rocketID[wRocket-4]}\n'
                                 f'{EpochStamp.strftime("%H:%M:%S")}\n'
                                 f'Ie0 Samples: {numPoints}\n'
                                 f'Sweep No: {sweep}')

                    plt.legend(['$nT_e^{1/2}$:'+f' {nTe}\n'
                                f'Vsp: {Vsp}V\n'
                                f'Te: {Te}eV\n'
                                f'n: {  (nTe/(np.sqrt(Te)))/(100*100*100) } $cm^-3$\n'
                                +'$\chi^{2}_{trans}$' + f'{ChiSquareTrans}\n'
                                '$\chi^{2}_{sat}$' + f'{ChiSquareSat}\n'])
                    plt.show()


        if outputData:
            # --- --- --- --- --- --- ---
            # --- WRITE OUT THE DATA ---
            # --- --- --- --- --- --- ---
            prgMsg('Creating output file')

            # remove unwanted data from output file
            sweeps_Epoch = np.array(sweeps_Epoch)
            sweeps_Vsp = np.array(sweeps_Vsp)
            sweeps_n0 = np.array(sweeps_n0)
            sweeps_Te = np.array(sweeps_Te)


            del data_dict['swept_Current'],data_dict['swept_Voltage'],data_dict['errorInProbeVoltage'],\
                data_dict['fixed_current'],data_dict['step_Voltage'], \
                data_dict['Epoch_step'],data_dict['Epoch_swept_Current']

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'
            outputCDFdata(outputPath,data_dict,outputModelData,globalAttrsMod,'Langmuir Probe')

            Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 0:  # ACES II Integration High
    rocketFolderPath = Integration_data_folder
    wflyer = 0
elif wRocket == 1:  # ACES II Integration Low
    rocketFolderPath = Integration_data_folder
    wflyer = 1
elif wRocket == 2:  # TRICE II High
    rocketFolderPath = TRICE_data_folder
    wflyer = 0
elif wRocket == 3:  # TRICE II Low
    rocketFolderPath = TRICE_data_folder
    wflyer = 1
elif wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5:  # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        main(wRocket, 0, rocketFolderPath, justPrintFileNames, wflyer,wSweeps)
    else:
        main(wRocket, 0, rocketFolderPath, justPrintFileNames, wflyer,wSweeps)
