# --- L2_mag_to_mSSA.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Takes as input the despun B-Field data or E-Field data,bandpass filters, then mSSAs the data.
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

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

# --- Select the DataSet ---
# 0 -> ACES II Electric Field Data
# 1 -> ACES II Magnetic FIeld Data
wData = 0

# --- Select the specific DataFile ---
wFile = 0

inputPath_modifier = 'l2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder

# Just print the names of files
justPrintFileNames = False

#################
# --- TOGGLES ---
#################
MirrorData = True
SSA_window_Size = 1201

# ---- SSA Components ----
computeSSAcomponents = True
breakIntoThisManyFiles = 0
# ------------------------

# ---- SSA Grouping and delta ----
plotGroupingSSA = True
# ------------------------

# --- OUTPUT DATA ---
outputData = True
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
            print('[{:.0f}] {:40s}{:5.1f} MB'.format(i, file.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\\', ''), round(getsize(file) / (10 ** 6), 1)))

        return

    # --- Get the Input Data ---
    prgMsg('Loading Data')
    data_dict = loadDictFromFile(inputDataFiles[wFile])
    Done(start_time)

    # --- determine the names of the components to mSSA ---
    keyNames = [['B_East', 'B_North', 'B_Up'], ['B_e', 'B_p', 'B_r'], ['Bx', 'By', 'Bz']]
    for keySet in keyNames: # check if set of keys exist in dictonary
        if data_dict.keys() >= {key for key in keySet}:
            compNames = keySet
            break
        raise Exception('No keyNames found in input Dictionary. Check component names in dictionary or consider defining new set of keys')



    ##############################
    # --- CALC mSSA Components ---
    ##############################
    if computeSSAcomponents:
        prgMsg('Calculating mSSA components')

        # --- break data into segments based on "BreakIntoThisManyFiles" variable


        from ACESII_code.class_var_func import mSSA_components
        data_dict_output = mSSA_components(data_dict_input=data_dict, compNames=compNames, SSA_window_Size=SSA_window_Size, mirrorData=MirrorData)

        Done(start_time)




    ##############################
    # --- OUTPUT THE DATA DICT ---
    ##############################
    if computeSSAcomponents:
        modifier = 'SSAcomponents_E' if wData == 0 else 'SSAcomponents_B'
    else:
        modifier = 'deltaE' if wData == 0 else 'deltaE'

    outputCDFdata(outputPathSSA, data_dict_SSAcomps, outputModelData, globalAttrsMod, 'RingCore')









# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
searchMod = 'E_Field_ENU' if wData == 0 else 'RingCore_ENU'
inputDataFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*{searchMod}*')

if len(inputDataFiles) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    mSSA_filtering(wRocket, wData, inputDataFiles, rocketFolderPath, justPrintFileNames, computeSSAcomponents)