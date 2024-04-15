# --- PoyntingFlux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Determine the PoyntingFLux of the data using E-Field and B-Field Measurements.
# ONLY works if E and B data is ALREADY the same length


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
printMagFiles = True
printElecFiles = True

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

modifier = ''
inputPath_modifier_elec = 'L2'
wMagFile = 2
Bscale = 1E-9 # what to multiply B-Field data to get into SI units

inputPath_modifier_mag = 'L2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
wEFIFile = 2
Escale = 1
outputPath_modifier = 'science/PoyntingFlux' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder


# --- --- --- PLOT --- --- ---
plotSPoynting = False
# --- --- --- OUTPUT --- --- ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---

from ACESII_code.class_var_func import u0,coordinatesNames,coordinatesSets

def PoyntingFlux(wRocket, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]

    inputFiles_elec = glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wRocket-4]}{modifier}\*E_Field*')
    inputFiles_mag = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wRocket-4]}{modifier}\*RingCore*')

    if justPrintFileNames:
        if printMagFiles:
            print('--- B-FIELD FILES ---')
            for i, file in enumerate(inputFiles_mag):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_mag[i], round(getsize(file) / (10 ** 6), 1)))
            print('\n')

        if printElecFiles:
            print('--- E-FIELD FILES ---')
            for i, file in enumerate(inputFiles_elec):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_elec[i], round(getsize(file) / (10 ** 6), 1)))
            print('\n')

        return

    print('\n')
    print(color.UNDERLINE + f'Calculating Poynting flux for ACESII {rocketID}' + color.END)

    # --- get the data from the mag file ---
    prgMsg(f'Loading data from mag Files')
    data_dict_mag = loadDictFromFile(inputFiles_mag[wMagFile])
    compNames_B,cordSysB,a = getCoordinateKeys(data_dict_mag)
    Done(start_time)

    # --- get the data from the electric file ---
    prgMsg(f'Loading data from Electric Field Files')
    data_dict_elec = loadDictFromFile(inputFiles_elec[wEFIFile])
    compNames_E,cordSysE,a = getCoordinateKeys(data_dict_elec)
    Done(start_time)

    #################################
    # --- CALCULATE POYNTING FLUX ---
    #################################
    prgMsg('Calculating Poynting Flux')

    # get the electric field and convert it to V/m, get the magnetic field and convert it to T
    B_Field = Bscale*np.array([data_dict_mag[compNames_B[0]][0],data_dict_mag[compNames_B[1]][0],data_dict_mag[compNames_B[2]][0]]).T
    E_Field = Escale*np.array([data_dict_elec[compNames_E[0]][0],data_dict_elec[compNames_E[1]][0],data_dict_elec[compNames_E[2]][0]]).T
    S = (1/u0)*np.array(np.cross(E_Field, B_Field))

    for i,vec in enumerate(S):
        if np.abs(sum(vec)) > 1E10:
            S[i] = np.array([rocketAttrs.epoch_fillVal,rocketAttrs.epoch_fillVal,rocketAttrs.epoch_fillVal])






    Done(start_time)
    if plotSPoynting:
        Epoch = data_dict_mag['Epoch'][0]
        fig, ax = plt.subplots(3)
        fig.suptitle('Poynting Flux')
        ax[0].plot(Epoch, S[:, 0])
        ax[0].set_ylabel('S_East')
        ax[1].plot(Epoch, S[:, 1])
        ax[1].set_ylabel('S_up')
        ax[2].plot(Epoch, S[:, 2])
        ax[2].set_ylabel('S_North')
        for i in range(3):
            ax[i].set_ylim(-3E-1,3E-1)
        plt.show()


    # --- prepare data for output ---
    prgMsg('Preparing Data')
    data_for_output = np.array([[S[i][0], S[i][1], S[i][2], np.linalg.norm(S[i])] for i in range(len(S))])
    compNamesS = [thing.replace('E','S') for thing in compNames_E]
    newComps = compNamesS + ['Smag']
    Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        prgMsg('creating output file')

        fileoutName = f'ACESII_{rocketID}_PoyntingFlux_{cordSysB}.cdf'

        data_dict_output = {}

        for i in range(len(newComps)):
            data = data_for_output[:, i]
            varAttrs = {'LABLAXIS': newComps[i], 'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None,
                        'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'W/m!A2!N',
                        'VALIDMIN': data.min(), 'VALIDMAX': data.max(),
                        'VAR_TYPE': 'data', 'SCALETYP': 'linear'}

            data_dict_output = {**data_dict_output, **{newComps[i]:[data, deepcopy(varAttrs)]}}

        Epoch_output = deepcopy(data_dict_mag['Epoch'])
        Epoch_output[1]['VAR_TYPE'] = 'support_data'
        data_dict_output = {**data_dict_output, **{'Epoch': Epoch_output}}

        # add in the attitude data
        keys = ['Alt', 'Lat', 'Long', 'Alt_geom', 'Lat_geom', 'Long_geom','ILat','ILong']
        for key in keys:
            data_dict_output = {**data_dict_output, **{key:data_dict_mag[key]}}

        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wRocket-4]}\\{fileoutName}'

        outputCDFdata(outputPath, data_dict_output,instrNam='PoyntingFlux')

        Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder


if len(glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wRocket-4]}\*.cdf')) == 0 and wRocket ==5:
    print(color.RED + 'There are no electric field .cdf files in the specified directory' + color.END)
elif len(glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no B-field .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        PoyntingFlux(wRocket, rocketFolderPath, justPrintFileNames)
    else:
        PoyntingFlux(wRocket, rocketFolderPath, justPrintFileNames)
