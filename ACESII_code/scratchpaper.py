# # --- scratchpaper.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering


# from ACESII_code.class_var_func import loadDictFromFile
# import numpy as np
#
# import os
#
# # os.environ['CMAKE_MAKE_PROGRAM'] =r'C:\msys64\mingw64.exe'
#
# inputFile = r'C:\Data\ACESII\attitude\high\ACESII_36359_Attitude_Solution.cdf'
# data_dict_attitude = loadDictFromFile(inputFile,{})
#
# Lat = data_dict_attitude['Latgd'][0]
# Long = data_dict_attitude['Long'][0]
# Alt = np.array(data_dict_attitude['Alt'][0])/1000
# Epoch = data_dict_attitude['Epoch'][0]

from ACESII_code.supportCode.Support_Libraries import pyglow
print(dir(pyglow))