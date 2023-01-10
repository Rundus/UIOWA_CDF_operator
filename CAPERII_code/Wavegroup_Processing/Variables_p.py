from cdflib import cdfread
import numpy as np

# from epoch_adjusting import epoch_adjusted_1,epoch_adjusted_2,epoch_adjusted_2_sec,epoch_adjusted_2_mins,epoch_adjusted_1_sec,epoch_adjusted_1_mins
from files_p import EEPAA_file_1,EEPAA_file_2,correlator_file,mag_file


# Reduce the data to spencer's 1min 27sec Window 490 - 577seconds after 9:27:00
t_rocket_start = 3772 # Corresponds to 9:27:00.035 in EEPAA DATA
t_start = 10943  # Corresponds to 9:32:59.998 in EEPAA DATA

# t_end = 15283 + 1 # Corresponds to 9:36:36.999, also the +1 is for python indexing
t_end = 17543 + 1 # Corresponds to 9:38:30.00 ( the +1 is for python indexing) in EEPAA DATA

t_rocket_start_correlator = 191825 # Corresponds to 9:27:00.000
t_start_correlator = 550579 # Corresponds to 9:33:00.000
t_end_correlator = 880581 + 1  # Corresponds to 9:38:30.000 ( the +1 is for python indexing)

# adjust = [epoch_adjusted_1,epoch_adjusted_2]
# adjust_min = [epoch_adjusted_1_mins,epoch_adjusted_2_mins]
# adjust_sec = [epoch_adjusted_1_sec,epoch_adjusted_2_sec]




EEPAA_file_1_info = EEPAA_file_1.cdf_info()
EEPAA_file_2_info = EEPAA_file_2.cdf_info()

zvars_1 = EEPAA_file_1_info['zVariables']
zvars_2 = EEPAA_file_1_info['zVariables']

epoch_1 = EEPAA_file_1.varget('Epoch')
epoch_2 = EEPAA_file_2.varget('Epoch')

counts_1 = EEPAA_file_1.varget('eepaa')
counts_2 = EEPAA_file_2.varget('eepaa')

pitches_1 = EEPAA_file_1.varget('Pitch_Angle')
pitches_2 = EEPAA_file_2.varget('Pitch_Angle')

energies_1 = EEPAA_file_1.varget('Energy')
energies_2 = EEPAA_file_2.varget('Energy')

geometric_1 = EEPAA_file_1.varget('geometric_factor')
geometric_2 = EEPAA_file_2.varget('geometric_factor')

count_interval_1 = EEPAA_file_1.varget('Count_Interval')  # Count is in units of micro-seconds and looks like [918 917 917 ...]
count_interval_2 = EEPAA_file_2.varget('Count_Interval')

count_interval_1_ns = [count_interval_1[x]*1000 for x in range(len(count_interval_1))]
count_interval_2_ns = [count_interval_2[x]*1000 for x in range(len(count_interval_2))]

ranges_1 = [range(len(counts_1)),range(len(counts_1[0])),range(len(counts_1[0][0]))]
ranges_2 = [range(len(counts_2)),range(len(counts_2[0])),range(len(counts_2[0][0]))]

# Get the working global attributes
globalAttrs_1 = EEPAA_file_1.globalattsget(expand=True)
globalAttrs_2 = EEPAA_file_2.globalattsget(expand=True)

#Get the FillValues for the eepaa's
fillval = EEPAA_file_1.varattsget(zvars_1[0])
fillval_1 = fillval['FILLVAL']
fillval = EEPAA_file_2.varattsget(zvars_2[0])
fillval_2 = fillval['FILLVAL']



#Corelator Information
correlator_info = correlator_file.cdf_info()
zvars_correlator = correlator_info['zVariables']
epoch_correlator = correlator_file.varget('Epoch')
correlators = correlator_file.varget('Correlators')
phase_bins = correlator_file.varget('Phase_Bin')
energy_correlator = correlator_file.varget('Correlator_Energy')


#Mag File Info
mag_info = mag_file.cdf_info()
mag_data = mag_file.varget('Magnetic_Field')
epoch_mag = mag_file.varget('Epoch')
zvars_mag = mag_info['zVariables']