# # --- scratchPaper.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering
from ACESII_code.class_var_func import setupPYCDF
setupPYCDF()


input = r'C:\Users\cfeltman\Desktop\sweptCals_Monday\LP_1p000k.csv'
import spacepy.pycdf as cdf

# Load the CSV file
data = cdf.from_csv('data.csv')

print(data)