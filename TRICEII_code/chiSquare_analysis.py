# --- chiSquare_corrections.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: organize the chiSquare data from chiSquare_padPair.py and apply the
# correction factor obtained from chiSquare_analysis.py to calibrate TRICEII data post-flight
import math


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
wRocket = []
# [] --> Both Fliers
# [0] --> High Flier
# [1] --> Low Flier


# big reference plot
Prin_vs_UnCal_size_plot = True

# individual plots
Prin_vs_UnCal_chi_plots = True


# --- --- --- --- ---
import time
from class_var_func import Done, setupPYCDF, prgMsg

start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- import ---
# --- --- --- ---
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from copy import deepcopy
from tqdm import tqdm
from itertools import product
from os import remove, path
from data_paths import fliers, TRICE_data_folder, COUNTS_p1_files, CHISQUARE_padPairs_p4_files
from class_var_func import color, L1_TRICE
from missionAttributes import TRICE_mission_dicts
from scipy.special import gammainc

setupPYCDF()
from spacepy import pycdf

print(color.BOLD + color.CYAN + 'chiSquare_analysis.py' + color.END + color.END)


def chiSquare_corrections(wRocket):
    print(color.UNDERLINE + f'{fliers[wRocket]} flier' + color.END)
    prgMsg('Acquiring data')

    # --- rocket attribute data ---
    rocketAttrs,missionAttrs, data_dict_temp = TRICE_mission_dicts()
    pycdf.lib.set_backward(False)

    # --- get the masked eepaa data ---
    data_dict_padPairs = {}
    inputFile = CHISQUARE_padPairs_p4_files[wRocket][0]
    with pycdf.CDF(inputFile) as cdfFile:
        for key, val in cdfFile.items():
            data_dict_padPairs = {**data_dict_padPairs, **{key: [np.array(cdfFile[key][...]), {key: val for key, val in cdfFile[key].attrs.items()}]}}

    padPairs = data_dict_padPairs['chiSquare_padPairs'][0]
    Pitch_Angles = np.array([10*(i - 1) for i in range(21)])

    Done(start_time)

    # --- --- --- --- --- --- --- --- ---
    # --- Organize calibration pairs ---
    # --- --- --- --- --- --- --- --- ---

    pairData = [  [ [] for j in range(21)] for i in range(21)]

    # append format:
    # [0] pad_principle,
    # [1] counts_principle,
    # [2] pad_uncal,
    # [3] counts_uncal,
    # [4] time_index
    # [5] energy_index
    for pair in tqdm(padPairs):
        pairData[int(pair[0]) ][int(pair[2])].append( [int(pair[1]), int(pair[3])] ) # pairData[ShouldBe_ptch].append([FoundIn, ShouldBe])

    # --- --- --- --- --- ---
    # --- Some analysis ---
    # --- --- --- --- --- ---

    if Prin_vs_UnCal_size_plot:
        pairSizeMatrix = np.zeros(shape=(21, 21))

        for pr in range(len(pairData)):
            for uc in range(len(pairData[0])):
                pairSizeMatrix[uc][pr] = int(len(pairData[uc][pr]))

        xData = Pitch_Angles
        yData = Pitch_Angles
        figure = plt.figure()
        plt.title('Principle Pad vs Uncal Pad # of padPairs ')
        if wRocket == 0:
            vmax = 30000
        elif wRocket == 1:
            vmax = 5000

        plt.imshow(pairSizeMatrix, interpolation='none', cmap=plt.cm.nipy_spectral,vmax=vmax,vmin=0,origin='lower')
        plt.xticks([i for i in range(21)])
        plt.yticks([i for i in range(21)])
        plt.ylabel('Principal Pad')
        plt.xlabel('Uncalibrated Pad')
        plt.colorbar().set_label('# of contacts', rotation=270)
        plt.grid(data = [[i -0.5 for i in range(22 + 1)],[i -0.5 for i in range(22 + 1)]])
        plt.savefig(rf'C:\Users\cfeltman\PycharmProjects\UIOWA_CDF_operator\TRICEII_code\extra\padPair_{fliers[wRocket]}.png', bbox_inches='tight', pad_inches=0)

    if Prin_vs_UnCal_chi_plots:

        # --- input data ---
        print('\n')
        prgMsg('Plotting Chi Square pair pads')
        for wPrinPad in tqdm(range(21)):
            for wUncalPad in range(21):
                padData = pairData[wPrinPad][wUncalPad]

                if padData != [] and len(padData) >= 10:
                    xData = [data[0] for data in padData]
                    yData = [data[1] for data in padData]

                    # --- analysis ---
                    roundTo = 5 # how many decimals to round to


                    # METHOD 1: Polyfit
                    dof = 2  # degrees of freedom in model
                    degree = 1
                    polyData, covariance = np.polyfit(xData, yData, degree,cov = True)
                    xData_m1_list = np.linspace(min(xData), max(xData),len(xData) + 1)
                    yData_m1 = [polyData[0]*xData_m1_list[i] for i in range(len(xData_m1_list))]
                    rChiSquare_m1 = round( (1/(len(xData) - dof)) * sum([(yData[i] - polyData[0]*xData[i])**2 / (math.sqrt(xData[i] + yData[i])) for i in range(len(xData))]), roundTo)
                    Qgof_poly = round(1 -  gammainc( (len(xData) - dof)/2 , rChiSquare_m1/2) , roundTo+10)

                    # METHOD 2: Calculus
                    dof = 1  # degrees of freedom in model
                    Sxy = sum([(xData[i] * yData[i]) / (xData[i] + yData[i]) for i in range(len(xData))])
                    Sxx = sum([(xData[i] * xData[i]) / (xData[i] + yData[i]) for i in range(len(xData))])
                    Sx  = sum([(xData[i] ) / (xData[i] + yData[i]) for i in range(len(xData))])
                    Sy  = sum([(yData[i]) / (xData[i] + yData[i]) for i in range(len(xData))])
                    S   = sum([ 1 / (xData[i] + yData[i]) for i in range(len(xData))])

                    alpha = round(Sxy / Sxx,5)
                    alphaError = round(1/math.sqrt(Sxx),roundTo)

                    xData_m2_list = np.linspace(min(xData),max(xData),len(xData) + 1)
                    yData_m2 = [alpha*xData_m1_list[i] for i in range(len(xData_m2_list))]
                    rChiSquare_m2 = round( (1 / (len(xData) - dof)) * sum([(yData[i] - alpha * xData[i]) ** 2 / (math.sqrt(xData[i] + yData[i])) for i in range(len(xData))]),roundTo)
                    Qgof_alpha = round(1 -  gammainc( (len(xData) - dof)/2 , rChiSquare_m2/2),roundTo+10)

                    # --- plot ---
                    fig, ax = plt.subplots()
                    fig.set_size_inches(10,10)
                    ax.plot(xData_m1_list, yData_m1, label= r'$\chi^{2}_{poly}$ ' + str(rChiSquare_m1))
                    ax.plot(xData_m2_list, yData_m2, label= r'$\chi^{2}_{\alpha}$ ' + str(rChiSquare_m2))
                    ax.scatter(xData,yData,color='black',s=3)
                    ax.set_title(f'Prin {Pitch_Angles[wPrinPad]}'+ '$^{\circ}$' + f'vs UnCal data {Pitch_Angles[wUncalPad]}' + '$^{\circ}$',fontsize=15)
                    ax.set_xlabel('UnCal Pad [Counts]',fontsize = 15)
                    ax.set_ylabel('Principle Pad [Counts]',fontsize=15)
                    ax.legend(loc='upper left')
                    ax.annotate(f'N = {len(xData)} \n' +
                                r'$Q_{poly}$ =' + f'{Qgof_poly}\n' +
                                r'$\alpha$ = ' + f'{alpha}$\pm${alphaError}\n' +
                                r'$Q_{\alpha}$ =' + f'{Qgof_alpha}',
                                xy=(0.01, 0.8),
                                xycoords = 'axes fraction',
                                fontsize = 10,
                                weight='bold'
                                )
                    fig.savefig(rf'C:\Users\cfeltman\PycharmProjects\UIOWA_CDF_operator\TRICEII_code\extra\padPair_Chi\padPair_Prin{Pitch_Angles[wPrinPad]}deg_UnCal{Pitch_Angles[wUncalPad]}deg.png',bbox_inches='tight' ,pad_inches=0.3)
                    plt.close()

        Done(start_time)



    # --- --- --- --- --- --- --- --- ---
    # --- Apply ChiSquare calibration ---
    # --- --- --- --- --- --- --- --- ---


    # --- --- --- --- --- --- --- ---
    # --- Prepare Data for Output ---
    # --- --- --- --- --- --- --- ---
    data_dict = {}
    data_dict = {**data_dict, **{'chiSquare_correctionFactors': [np.array([0]), deepcopy(data_dict_temp['data'][1])]}}


    # --- --- --- --- ---
    # --- OUTPUT Data ---
    # --- --- --- --- ---

    prgMsg('Writing out Data')
    fileoutName = f'TRICE_{rocketAttrs.rocketID[wRocket]}_correctionFactors'
    outputPath = f'{TRICE_data_folder}\\science\\{fliers[wRocket]}\\{fileoutName}_p5.cdf'

    # --- delete output file if it already exists ---
    if path.exists(outputPath):
        remove(outputPath)
    pycdf.lib.set_backward(False)

    # --- open the output file ---
    with pycdf.CDF(outputPath, '') as cdfFile:
        cdfFile.readonly(False)

        # --- WRITE OUT GLOBAL ATTRS ---
        globalAttrsMod = {'nothing': None}
        ModelData = L1_TRICE(wRocket)[0]
        inputGlobDic = ModelData.cdfFile.globalattsget()

        for key, val in inputGlobDic.items():
            if key in globalAttrsMod:
                cdfFile.attrs[key] = globalAttrsMod[key]
            else:
                cdfFile.attrs[key] = val

        # --- WRITE OUT DATA ---
        for varKey, varVal in data_dict.items():
            if len(varVal[0]) != 0:
                if varKey == 'Epoch':
                    cdfFile.new(varKey, data=varVal[0], type=33)
                else:
                    cdfFile.new(varKey, data=varVal[0])

                # --- Write out the attributes and variable info ---
                for attrKey, attrVal in data_dict[varKey][1].items():

                    if attrKey == 'VALIDMIN':
                        cdfFile[varKey].attrs[attrKey] = varVal[0].min()
                    elif attrKey == 'VALIDMAX':
                        cdfFile[varKey].attrs[attrKey] = varVal[0].max()
                    elif attrVal != None:
                        cdfFile[varKey].attrs[attrKey] = attrVal
    Done(start_time)
    print('\n')


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

if wRocket == []:
    for i in range(2):
        if path.exists(CHISQUARE_padPairs_p4_files[i][0]):
            chiSquare_corrections(i)
        else:
            print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if path.exists(CHISQUARE_padPairs_p4_files[wRocket[0]][0]):
        chiSquare_corrections(wRocket[0])
    else:
        print(color.RED + 'There are no .cdf files in the specified directory' + color.END)