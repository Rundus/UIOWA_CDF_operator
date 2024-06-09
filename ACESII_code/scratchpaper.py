from ACESII_code.myImports import *

ACESI_EEPAA_path = r'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\PlotExtra\BetaFit\ACESI_EEPAA.cdf'
data_dict_ACESI = loadDictFromFile(ACESI_EEPAA_path)
print(data_dict_ACESI['energy_cal'][0])