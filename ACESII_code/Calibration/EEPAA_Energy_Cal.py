# --- EEPAA_Energy_Cal.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Using the 401 chamber calibration data, fit a linear curve to
# between the data to determine a calibration curve for EEPAA steps/energy
import matplotlib.pyplot as plt
import numpy as np

# --- --- --- --- ---
from ACESII_code.myImports import *

start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
inputFile = r'C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\Calibration\EEPAA_response_comparison.xlsx'

wRow = 0 # row corresponding to where all the names of the variables are


plotCalData_vs_Theoretical = False

# Thereotical EEPAA Energies
# 12750.00, 10905.65, 9328.10, 7978.75, 6824.58, 5837.38, 4992.97, 4270.72,
#                           3652.94,  3124.52, 2672.55, 2285.95, 1955.28, 1672.44, 1430.51, 1223.58,
#                           1046.58,   895.19,  765.70,  654.94,  560.20,  479.16,  409.85,  350.56,
#                            299.85,   256.48,  219.38,  187.64,  160.50,  137.28,  117.42,  100.44,
#                             85.91,    73.48,   62.85,   53.76,   45.98,   39.33,   33.64,   28.78,
#                             24.61,    21.05,   18.01,   15.40,   13.17,   11.27,    9.64,    8.24, 7.05


# --- --- --- ---
# --- import ---
# --- --- --- ---

import pandas as pd


def EEPAA_Energy_Cal(inputFile):

    # collect the csv data
    df = pd.read_excel(inputFile)
    steps = np.array(df.STEP)
    Energy_Theory = np.array(df['Energy (Theoretical) [eV]'])
    Energy_cal = np.array(df['Calibration Estimate [eV]'])

    dataForCal = np.array([list(pair) for pair in zip(Energy_Theory,Energy_cal) if not (np.isnan(pair[0]) or np.isnan(pair[1]))])

    xData = dataForCal[:,0]
    yData = dataForCal[:,1]


    def fitFunc(x,A):
        return A*x

    params, cov = curve_fit(fitFunc,xData,yData)

    xDataFit = np.linspace(min(xData),max(xData),200)
    yDataFit = np.array([fitFunc(x,*params) for x in xDataFit])

    if plotCalData_vs_Theoretical:
        fig, ax = plt.subplots()
        fig.suptitle('Theoretical Energy vs Chamber Calibration Estimate')
        ax.scatter(xData,yData,label='Raw',color='blue')
        ax.plot(xDataFit, yDataFit, label=f'Fit Function: y = A*x\n A={params[0]}',color='red')
        ax.set_ylabel('Calibration Estimate [eV]')
        ax.set_xlabel('Energy (Theoretical) [eV]')
        ax.legend()
        plt.show()

    new_EEPAA_Energies = list(Energy_Theory*params[0])[2:]
    new_EEPAA_Energies.append(Energy_Theory[0]*params[0])
    new_EEPAA_Energies = [round(engy,2) for engy in new_EEPAA_Energies]
    print(new_EEPAA_Energies[::-1])








# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
EEPAA_Energy_Cal(inputFile)