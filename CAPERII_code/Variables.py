# -----------------------------------------------------------------------------------------------
# A Variables File that contains everything relevant in the parent eepaa1 and eepaa2 counts files
# -----------------------------------------------------------------------------------------------

from cdflib import cdfread
import numpy as np
from files import EEPAA_file_1,EEPAA_file_2,mag_file,correlator_file
from files import diffFlux_file_1,diffFlux_file_2
from files import pitch_actual_file_1,pitch_actual_file_2
from files import flux_aligned_file_1,flux_aligned_file_2
# from files import dist_file_1,dist_file_2
from files import EEPAA_file_high

# ------------------------
# INITIAL VARIABLES CAPERII
# ------------------------

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

# ---------------------------
# PROCESSED VARIABLES CAPERII
# ---------------------------

diffFlux_1_info = diffFlux_file_1.cdf_info()
diffFlux_2_info = diffFlux_file_2.cdf_info()

zvars_dflux_1 = diffFlux_1_info['zVariables']
zvars_dflux_2 = diffFlux_2_info['zVariables']

diffEflux_1 = diffFlux_file_1.varget('Differential_Energy_Flux')
diffNflux_1 = diffFlux_file_1.varget('Differential_Number_Flux')

diffEflux_2 = diffFlux_file_2.varget('Differential_Energy_Flux')
diffNflux_2 = diffFlux_file_2.varget('Differential_Number_Flux')

pitchNAM = 'Magnetometer_Pitch_Angle'

pitch_angle_actual_1 = pitch_actual_file_1.varget(pitchNAM)
pitch_angle_actual_2 = pitch_actual_file_2.varget(pitchNAM)

pitch_actual_1_info = pitch_actual_file_1.cdf_info()
pitch_actual_2_info = pitch_actual_file_2.cdf_info()

zvars_pitch_1 = pitch_actual_1_info['zVariables']
zvars_pitch_2 = pitch_actual_2_info['zVariables']

# Nflux_aligned_1 = flux_aligned_file_1.varget('Differential_Number_Flux_Aligned')
# Eflux_aligned_1 = flux_aligned_file_1.varget('Differential_Energy_Flux_Aligned')
# Nflux_aligned_2 = flux_aligned_file_2.varget('Differential_Number_Flux_Aligned')
# Eflux_aligned_2 = flux_aligned_file_2.varget('Differential_Energy_Flux_Aligned')
#
# flux_aligned_1_info = flux_aligned_file_1.cdf_info()
# flux_aligned_2_info = flux_aligned_file_2.cdf_info()
#
# zvars_flux_aligned_1 = flux_aligned_1_info['zVariables']
# zvars_flux_aligned_2 = flux_aligned_2_info['zVariables']
#
# dist_1 = dist_file_1.varget('Distribution_Function')
# dist_2 = dist_file_2.varget('Distribution_Function')
#
# dist_1_info = dist_file_1.cdf_info()
# dist_2_info = dist_file_2.cdf_info()
#
# zvars_dist_1 = dist_1_info['zVariables']
# zvars_dist_2 = dist_2_info['zVariables']

# ------------------------------
# General Variables MAGNETOMETER
# ------------------------------

mag_file_info = mag_file.cdf_info()
epoch_mag = mag_file.varget('Epoch')
zvars_mag = mag_file_info['zVariables']
mag_dat = mag_file.varget('Magnetic_Field')

# -------------------------
# General Variables TRICEII
# -------------------------
EEPAA_file_high_info = EEPAA_file_high.cdf_info()
zvars_high =EEPAA_file_high_info['zVariables']
count_interval_high = EEPAA_file_high.varget('count_interval')
TRICE_deadtime = EEPAA_file_high.varget('dead_time')
# TRICE_pitch_actual = EEPAA_file_high.varget('Pitch_Angle_Actual')


