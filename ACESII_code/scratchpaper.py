# # --- scratchpaper.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering


# imports
from ACESII_code.myImports import *


inputFile = 'C:\Data\ACESII\L2\high\ACESII_36359_l2_eepaa_fullCal.cdf'

data_dict = loadDictFromFile(inputFile,{})

inputFile1 = 'C:\Data\ACESII\L2\ACESII_36364_l2_eepaa_fullCal_kentons.cdf'

data_dict_kenton = loadDictFromFile(inputFile1,{})
for i in range(10):
    print(data_dict['Epoch_esa'][0][i],data_dict_kenton['Epoch_esa'][0][i])

