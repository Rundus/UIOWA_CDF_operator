# --- L2_mag_to_mSSA.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Takes as input the despun B-Field data or E-Field data with/without bandpass filters. Then mSSAs the data.
# The output are SSA component files found in \L3\SSAcomponents_B
# or deltaB/deltaE files.

# it this code must be able to do the following SIMPLY
# [0] handle Electric or magnetic data
# [1] perform mSSA on an input data_dict
# [2] be able to output components for the whole time series, a subset or the whole series broken into subsets
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
# 0 -> ACES II Electric Field Data
# 1 -> ACES II Magnetic Field Data
wData = 0
searchMod = 'E_Field' if wData == 0 else 'RingCore'

# --- Select the specific DataFile ---
wFile = 2


#################
# --- TOGGLES ---
#################

# --- Mirror Data ---
# for calculating the components, the data is cut in half
MirrorData = True

# --- Window Size ---
SSA_window_Size = 601

# ---- SSA Components ----
computeSSAcomponents = False
breakIntoThisManyFiles = 5
generate_only_these_files = [] # if I only want to re-generate specific component subsets
# ------------------------

# ---- SSA Grouping and delta ----
plotGroupingSSA = True
wAxesSSA = 1
groupings_dict = {'RingCore_high': [[i for i in range(16)]+
                                    [20,21,22,23,24,25,26,28,29,30,31, 32,37,38,39,40,41,43,44,45,51,52, 53,54,56,57,58,59] +
                                    [60+i for i in range(10)] +
                                    [71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 85, 86, 88, 89]+
                                    [90, 91, 92, 93, 94, 95, 96, 98,99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110]+
                                    [], # identified harmonics (10,11,48,51,81,82 maybe?)
                                   [[i] for i in range(14,18)], # indicies to investigation
                                   200], # limit of the noise, between this value and 3*SSA_window size is discarded as noise
                  'E_Field_low': [[0,1,2,3,4,5,8,9],
                                   [[i] for i in range(10,20)],
                                   3*SSA_window_Size],
                  'RingCore_low': [[0,1,2,3,4,5],
                                   [[i] for i in range(0,10)],
                                   3*SSA_window_Size]
                  }
# ------------------------

# --- OUTPUT DATA ---
outputData = False
# -------------------


def mSSA_filtering(wRocket, wData, inputDataFiles, rocketFolderPath, justPrintFileNames, computeSSAcomponents):

    # --- ACES-II Attributes Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    globalAttrsMod = rocketAttrs.globalAttributes[wRocket-4]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L3'
    outputModelData = L0_ACES_Quick(wRocket-4)

    if justPrintFileNames:
        for i, file in enumerate(inputDataFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, file.rsplit('\\')[-1], round(getsize(file) / (10 ** 6), 1)))

        return

    # --- Get the Input Data ---
    prgMsg('Loading Data for ' + inputDataFiles[wFile].rsplit("\\")[-1])
    data_dict = loadDictFromFile(inputDataFiles[wFile])
    Done(start_time)

    # --- Determine which coordinate system the data uses ---
    from ACESII_code.class_var_func import determineCoordCompNames
    coordSys, compNames, fileWL = determineCoordCompNames(inputDataFiles[wFile])

    ##############################
    # --- CALC mSSA Components ---
    ##############################
    if computeSSAcomponents:

        # --- break data into segments based on "BreakIntoThisManyFiles" variable ---
        # find the inidcies to break up the dataset into "breakIntoThisManyPieces" parts
        rangeLimits = [[ar[0], ar[-1]+1] for ar in np.array_split([i for i in range(len(data_dict['Epoch'][0]))], breakIntoThisManyFiles)]

        # modify rangeLimits to only include those in "generate_these_files"
        if generate_only_these_files != []:
            rangeLimits = [rangeLimits[num] for num in generate_only_these_files]


        # --- DO THE mSSA component Calculation ---
        from ACESII_code.class_var_func import mSSA_components

        for i, subset in enumerate(rangeLimits):

            prgMsg(f'Calculating SSA components for Indicies {subset[0]} to {subset[1]}')

            # get the data_dict subset
            temp_dict = deepcopy(data_dict)
            data_dict_chunk = {key: [val[0][subset[0]:subset[1]], val[1]] for key, val in temp_dict.items()}
            data_dict_output = mSSA_components(data_dict_input=data_dict_chunk, compNames=compNames,
                                               SSA_window_Size=SSA_window_Size, mirrorData=MirrorData,
                                               mirrorPercent=0.2)

            # --- --- --- --- ---
            # --- OUTPUT DATA ---
            # --- --- --- --- ---
            pathModifier = 'SSAcomponents_E' if wData == 0 else 'SSAcomponents_B'
            wInstr = 'E_Field' if wData == 0 else 'RingCore'
            outputFilePath = rf'{rocketFolderPath}L3\\{pathModifier}\\{fliers[wRocket-4]}\\ACESII_{rocketID}_{wInstr}_SSAComponents_{coordSys}_WL{SSA_window_Size}_subset_{i}'
            if MirrorData:
                outputFilePath = outputFilePath + '_mirrored'
            outputCDFdata(outputFilePath+".cdf", data_dict_output, outputModelData, globalAttrsMod, pathModifier)
            Done(start_time)


    #############################
    # --- DEFINE THE GROUPING ---
    #############################


    # generate the groupings variables
    from ACESII_code.class_var_func import generateGrouping
    groupings = generateGrouping(SSA_window_Size = fileWL,
                                 badCompIndicies=groupings_dict[searchMod+f'_{fliers[wRocket-4]}'][0],
                                 InvestigateIndicies = groupings_dict[searchMod+f'_{fliers[wRocket-4]}'][1],
                                 noiseLimit=groupings_dict[searchMod+f'_{fliers[wRocket-4]}'][2])


    ###############################
    # --- PLOTTING THE GROUPING ---
    ###############################
    if plotGroupingSSA:

        # plot the grouping
        from ACESII_code.class_var_func import mSSA_grouping_Plotting
        B_SSA_vec = [data_dict[compNames[i]][0] for i in range(len(compNames))]
        mSSA_grouping_Plotting(B_SSA=B_SSA_vec,
                               groupings=groupings,
                               compNames=compNames,
                               SSA_window_Size=fileWL,
                               wAxesSSA = wAxesSSA,
                               wRocket = wRocket,
                               Epoch=data_dict['Epoch'][0],
                               mirrored = MirrorData)


    ##############################
    # --- OUTPUT THE DATA DICT ---
    ##############################

    if outputData:
        prgMsg('Outputting Data')

        # --- load all the data ---
        allData = []
        for file in inputDataFiles:
            dict =deepcopy(loadDictFromFile(file))
            allData.append(dict)

        # --- apply the mSSA grouping to the data ---
        newComps = ['d' + comp for comp in compNames]
        newData = [[] for i in range(len(newComps))] # should look like [[],[],[]]

        # loop through all the dictionaries and apple the groupings
        for i, SSAdict in enumerate(allData):

            for k in range(len(newComps)):  # loop through all the components, but only take the last grouping, which is the physical signal

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
                for k in range(len(newComps)):
                    newData[k][i] = newData[k][i][startPoint:endPoint]

        # --- combine all the subsets into one dataset ---
        stitchedData = np.array([[item for sublist in dataset for item in sublist] for dataset in newData])
        stitchedEpoch = np.array([item for dict in allData for item in dict['Epoch'][0]])

        data_dict_output = {'Epoch':[stitchedEpoch, allData[0]['Epoch'][1]],
                            newComps[0]:[np.array(stitchedData[0]), allData[0][compNames[0]][1]],
                            newComps[1]:[np.array(stitchedData[1]), allData[0][compNames[1]][1]],
                            newComps[2]:[np.array(stitchedData[2]), allData[0][compNames[2]][1]]}

        # --- output the data ---
        filePaths = ['deltaE','deltaB']
        outputPath = f'{rocketFolderPath}\\L3\\{filePaths[wData]}\\{fliers[wRocket-4]}\\ACESII_{rocketID}_{searchMod}_{coordSys}_WL{fileWL}_stitchedFlight.cdf'
        outputCDFdata(outputPath, data_dict_output,outputModelData,globalAttrsMod,searchMod)
        Done(start_time)









# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
if computeSSAcomponents:
    inputDataFiles = glob(f'{rocketFolderPath}\\l2\\{fliers[wRocket-4]}\\*{searchMod}*')
elif plotGroupingSSA or outputData:
    inputPath_modifier = 'SSAcomponents_E' if wData == 0 else 'SSAcomponents_B'
    inputDataFiles = glob(f'{rocketFolderPath}\\l3\\{inputPath_modifier}\\{fliers[wRocket - 4]}\\*{searchMod}*')
else:
    raise Exception('Need to set computeSSAcomponents == True or plotGroupingSSA == True')

if len(inputDataFiles) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    mSSA_filtering(wRocket, wData, inputDataFiles, rocketFolderPath, justPrintFileNames, computeSSAcomponents)