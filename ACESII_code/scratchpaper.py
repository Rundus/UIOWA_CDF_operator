from myspaceToolsLib.CDF_load import loadDictFromFile
fileLocation = 'https://space.physics.uiowa.edu/rockets/data/ACESII_36359/L2/Flight/ACESII_36359_l2_eepaa_fullCal.cdf'

data_dict = loadDictFromFile(fileLocation)
for key in data_dict.keys():
    print(key)