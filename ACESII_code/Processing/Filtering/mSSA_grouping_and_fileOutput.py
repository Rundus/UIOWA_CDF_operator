# --- mSSA_grouping_and_fileOutput.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:
# The output are SSA component files found in \L3\SSAcomponents_B
# or deltaB/deltaE files.

# it this code must be able to do the following SIMPLY
# [3] Plot the component grouping step
# [4] produce the deltaB/deltaE datafiles through only a couple toggles

# --- bookkeeping ---
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# --- Select the DataSet ---
inputPath_modifier = 'L3\SSAcomponents_S'


#################
# --- TOGGLES ---
#################

# --- Plotting ---
plot_GroupingSSA = False
plot_wSSAComponetFile = 2
plot_wAxesSSA = 1
# --- OUTPUT DATA ---
outputDataIntoOneMasterFile = True
# -------------------

def mSSA_grouping_and_fileOutput(wRocket, rocketFolderPath, justPrintFileNames):

    # --- ACES-II Attributes Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    inputFiles, input_names, input_names_searchable = getInputFiles(rocketFolderPath=rocketFolderPath,wRocket=wRocket,inputPath_modifier=inputPath_modifier)

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:95s}{:5.1f} MB'.format(i, file.rsplit('\\')[-1], round(getsize(file) / (10 ** 6), 1)))
        return

    # --- Get the Input Data for one File to read its contents---
    data_dict, globalAttrs = loadDictFromFile(inputFiles[plot_wSSAComponetFile], getGlobalAttrs=True)

    # --- Determine which coordinate system the data uses ---
    compNames, coordSys, coordSet = getCoordinateKeys(data_dict)
    varName = compNames[0].replace(coordSet[0], '')

    wInstr = ''
    for vari in ['RingCore', 'PoyntingFlux', 'mBeEr', 'BrEe', 'E_Field']:
        if vari in inputFiles[plot_wSSAComponetFile]:
            wInstr = vari
            break

    # DETERMINE SSA Window Length
    firstIndex = input_names_searchable[plot_wSSAComponetFile].find('_WL')
    lastIndex = input_names_searchable[plot_wSSAComponetFile].find('_subset')
    SSA_window_Size = int(input_names_searchable[plot_wSSAComponetFile][firstIndex+3:lastIndex])

    # DETERMINE if data is mirror
    MirrorData = True if 'mirrored' in input_names_searchable[plot_wSSAComponetFile] else False


    #############################
    # --- DEFINE THE GROUPING ---
    #############################
    groupings_dict = {'RingCore_high': [[i for i in range(16)] +
                                        [20, 21, 22, 23, 24, 25, 26, 28, 29, 30, 31, 32, 37, 38, 39, 40, 41, 43, 44, 45,
                                         51, 52, 53, 54, 56, 57, 58, 59] +
                                        [60 + i for i in range(10)] +
                                        [71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 85, 86, 88, 89] +
                                        [90, 91, 92, 93, 94, 95, 96, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107,
                                         108, 109, 110] +
                                        [],  # identified harmonics (10,11,48,51,81,82 maybe?)
                                        [[i] for i in range(14, 18)],  # indicies to investigation
                                        200],
                      # limit of the noise, between this value and 3*SSA_window size is discarded as noise
                      'E_Field_low': [[0, 1, 2, 3, 4, 5, 8, 9],
                                      [[i] for i in range(10, 20)],
                                      3 * SSA_window_Size],
                      'RingCore_low': [[0, 1, 2, 3, 4, 5],
                                       [[i] for i in range(0, 10)],
                                       3 * SSA_window_Size],
                      'PoyntingFlux_low': [[0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 14],
                                           [[i] for i in range(5, 11)],
                                           3 * SSA_window_Size]
                      }


    # generate the groupings variables
    from ACESII_code.class_var_func import generateGrouping
    groupings = generateGrouping(SSA_window_Size = SSA_window_Size,
                                 badCompIndicies=groupings_dict[wInstr+f'_{fliers[wRocket-4]}'][0],
                                 InvestigateIndicies = groupings_dict[wInstr+f'_{fliers[wRocket-4]}'][1],
                                 noiseLimit=groupings_dict[wInstr+f'_{fliers[wRocket-4]}'][2])


    ###############################
    # --- PLOTTING THE GROUPING ---
    ###############################
    if plot_GroupingSSA:

        # plot the grouping
        from ACESII_code.class_var_func import mSSA_grouping_Plotting
        B_SSA_vec = [data_dict[compNames[i]][0] for i in range(len(compNames))]
        mSSA_grouping_Plotting(B_SSA=B_SSA_vec,
                               groupings=groupings,
                               compNames=compNames,
                               SSA_window_Size=SSA_window_Size,
                               wAxesSSA = plot_wAxesSSA,
                               wRocket = wRocket,
                               Epoch=data_dict['Epoch'][0],
                               mirrored = MirrorData)

    ##############################
    # --- OUTPUT THE DATA DICT ---
    ##############################

    if outputDataIntoOneMasterFile:
        prgMsg('Outputting Data')

        # --- load all the data ---
        allData = []
        for file in inputFiles:
            dict =deepcopy(loadDictFromFile(file))
            allData.append(dict)

        # --- apply the mSSA grouping to the data ---
        newData = [[] for i in range(len(compNames))] # should look like [[],[],[]]

        # loop through all the dictionaries and apple the groupings
        for i, SSAdict in enumerate(allData):

            for k in range(len(compNames)):  # loop through all the components, but only take the last grouping, which is the physical signal

                data = np.array(SSAdict[compNames[k]][0])

                # combine the groupings of the last group set only
                formattedData = np.array([0 for i in range(len(data[:, 0]))], dtype='float64')
                for j, componentIndex in enumerate(groupings[-1]):
                    formattedData += np.array(data[:, componentIndex])

                newData[k].append(formattedData)

        # --- eliminate the 20% padded data on the sides ---
        if MirrorData:

            # loop through all the subset dictionaries
            for i, SSAdict in enumerate(allData):

                # for each dictionary, determine what the start/end point of the real data is
                dataTrueLength = len(SSAdict['Epoch'][0])
                startPoint = int(0.2*dataTrueLength)
                endPoint = startPoint + dataTrueLength

                # apply the reduction to the subsets
                for k in range(len(compNames)):
                    newData[k][i] = newData[k][i][startPoint:endPoint]

        # --- combine all the subsets into one dataset ---
        stitchedData = np.array([[item for sublist in dataset for item in sublist] for dataset in newData])
        stitchedEpoch = np.array([item for dict in allData for item in dict['Epoch'][0]])

        data_dict_output = {'Epoch': [stitchedEpoch, allData[0]['Epoch'][1]],
                            compNames[0]:[np.array(stitchedData[0]), allData[0][compNames[0]][1]],
                            compNames[1]:[np.array(stitchedData[1]), allData[0][compNames[1]][1]],
                            compNames[2]:[np.array(stitchedData[2]), allData[0][compNames[2]][1]]}

        # --- output the data ---
        outputPath = f'{rocketFolderPath}\\L3\\delta{varName}\\{fliers[wRocket-4]}\\ACESII_{rocketID}_{wInstr}_{coordSys}_WL{SSA_window_Size}.cdf'
        outputCDFdata(outputPath, data_dict_output)
        Done(start_time)









# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    mSSA_grouping_and_fileOutput(wRocket, rocketFolderPath, justPrintFileNames)
