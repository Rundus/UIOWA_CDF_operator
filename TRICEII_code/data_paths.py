# --- data_paths.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to store the pathing information of the project


# --- imports ---
from glob import glob
from os import path



# --- --- --- --- --- --- ---
# --- USER SPECIFIC DATA ---
# --- --- --- --- --- --- ---
if path.exists(r'D:\Data\\'):
    PATH_TO_DATA_FOLDER = r'D:\Data\\'
elif path.exists(r'C:\Users\cfeltman\PycharmProjects\UIOWA_CDF_operator\Data\\'):
    PATH_TO_DATA_FOLDER = r'C:\Users\cfeltman\PycharmProjects\UIOWA_CDF_operator\Data\\'

CDF_LIB = r"C:/Users/cfeltman/PycharmProjects/UIOWA_CDF_operator/CDF/lib"  # location to directory that contains the dllcdf.dll needed for pycdf the library



HOMEDRIVE = 'C:'
HOMEPATH = 'C:\\'
fliers = ['high','low']

# --- --- --- --- --- ---
# --- TRICE II PATHS ---
# --- --- --- --- --- ---
TRICE_data_folder = fr'{PATH_TO_DATA_FOLDER}TRICEII\\'
TRICE_science_folder = fr'{TRICE_data_folder}\science\\'
TRICE_attitude_folder= fr'{TRICE_data_folder}\attitude\\'
TRICE_ACS_files = [glob(f'{TRICE_data_folder}attitude\{fliers[0]}\*.cdf'), glob(f'{TRICE_data_folder}attitude\{fliers[1]}\*.cdf')]
TRICE_tad_files = [glob(f'{TRICE_data_folder}tad\{fliers[0]}\*.tad'),glob(f'{TRICE_data_folder}tad\{fliers[1]}\*.tad')]
TRICE_tmCDF_files = [glob(f'{TRICE_data_folder}tmCDF\{fliers[0]}\*.cdf'),glob(f'{TRICE_data_folder}tmCDF\{fliers[1]}\*.cdf')]
TRICE_L0_files = [glob(f'{TRICE_data_folder}L0\{fliers[0]}\*.cdf'),glob(f'{TRICE_data_folder}L0\{fliers[1]}\*.cdf')]
TRICE_L1_files = [glob(f'{TRICE_data_folder}L1\{fliers[0]}\*.cdf'),glob(f'{TRICE_data_folder}L1\{fliers[1]}\*.cdf')]
TRICE_L2_files = [glob(f'{TRICE_data_folder}L2\{fliers[0]}\*.cdf'),glob(f'{TRICE_data_folder}L2\{fliers[1]}\*.cdf')]
TRICE_mag_files = [glob(f'{TRICE_data_folder}mag\{fliers[0]}\*.cdf'),glob(f'{TRICE_data_folder}mag\{fliers[1]}\*.cdf')]
TRICE_attitude_csv_files = [glob(f'{TRICE_data_folder}attitude\{fliers[0]}\*.csv'),glob(f'{TRICE_data_folder}attitude\{fliers[1]}\*.csv')]
TRICE_attitude_cdf_files = [glob(f'{TRICE_data_folder}attitude\{fliers[0]}\*.cdf'),glob(f'{TRICE_data_folder}attitude\{fliers[1]}\*.cdf')]

# --- --- --- --- --- ---
# --- Science Files ---
# --- --- --- --- --- ---
COUNTS_p1_files = [glob(f'{TRICE_science_folder}\high\*_p1.cdf'),glob(f'{TRICE_science_folder}\low\*_p1.cdf')]
COUNTS_p2_files = [glob(f'{TRICE_science_folder}\high\*_p2.cdf'),glob(f'{TRICE_science_folder}\low\*_p2.cdf')]
CALC_pitch_angle_p3_files = [glob(f'{TRICE_science_folder}\high\*_p3.cdf'),glob(f'{TRICE_science_folder}\low\*_p3.cdf')]
CHISQUARE_padPairs_p4_files = [glob(f'{TRICE_science_folder}\high\*_p4.cdf'),glob(f'{TRICE_science_folder}\low\*_p4.cdf')]