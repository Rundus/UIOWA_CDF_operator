# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ---- TRICE Magnetometer Handling ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import cdflib,numpy as np,math,time
from cdflib import cdfwrite


start_time = time.time()
print('importing variables: ', end='')

# ------------
# Import Files
# ------------
from files import Mag_file_high,Mag_file_low,EEPAA_file_high,EEPAA_file_low
from functions import take_closest,UTCtoTT2000, getMagdat
from Variables import EEPAA_file_low_info,EEPAA_file_high_info,zvars_high,zvars_low,EPOCH_High,EPOCH_low,pitch,count_interval_low,count_interval_high,fillval_high,fillval_low


cdf_Pitch_Actual_high = cdfwrite.CDF('C:/Users/Rundus/Desktop/TRICE/pitch_actual_high',cdf_spec=EEPAA_file_high_info,delete=True)
cdf_Pitch_Actual_low = cdfwrite.CDF('C:/Users/Rundus/Desktop/TRICE/pitch_actual_low',cdf_spec=EEPAA_file_low_info,delete=True)

# Assign the magnetometer data in (Bx,By,Bz) form
mag_high = getMagdat(Mag_file_high)
mag_low = getMagdat(Mag_file_low)

MagXdata_high = Mag_file_high.varget('MagX')
MagXdata_low = Mag_file_low.varget('MagX')

MAG_high_time = np.transpose( [MagXdata_high[i][0] for i in range(len(MagXdata_high))])
MAG_low_time = np.transpose([MagXdata_low[i][0] for i in range(len(MagXdata_low))])

print('Done')

# -----------------------------------------------------------------
# -----CONVERT MAGNETOMETER AND TRICE DATA TO READABLE FORMATS-----
# -----------------------------------------------------------------
print('Converting magnetometer times: ', end='')

# TT2000 - year, month, day, hour, minute, second, millisecond, microsecond, nanosecond
breakdown_TRICE_high = cdflib.cdfepoch.breakdown_tt2000(EPOCH_High, to_np=False)
breakdown_TRICE_low = cdflib.cdfepoch.breakdown_tt2000(EPOCH_low, to_np=False)

# UTC to TT2000
MAG_high_time_converted =[UTCtoTT2000(MAG_high_time[i]) for i in range(len(MAG_high_time)) ]
MAG_low_time_converted =[UTCtoTT2000(MAG_low_time[i]) for i in range(len(MAG_low_time)) ]

# ------------------------------------------------------------------------------
# ----------Match the TRICE epoch stamp to the best magnetometer value----------
# ------------------------------------------------------------------------------
# Convert mag time to TT2000 computed values
MAG_high_time_computed = cdflib.cdfepoch.compute_tt2000(MAG_high_time_converted)
MAG_low_time_computed = cdflib.cdfepoch.compute_tt2000(MAG_low_time_converted)



# ----------------------------------------------
# ----- MAKE THE MAGNETOMETER INDEX MATRIX -----
# ----------------------------------------------

# Correlate TRICE and magnetometer timestamps

MAG_matrix_index_high = np.zeros([len(EPOCH_High),49])
MAG_matrix_index_low = np.zeros([len(EPOCH_low),49])

# ASSIGN THE MAGNETIC FIELD VALUES
print('Done')
print('Calculating Magnetometer Indexes for High: ', end='')

for i in range(len(EPOCH_High)):
    MAG_matrix_index_high[i] = [ take_closest(MAG_high_time_computed, EPOCH_High[i] + k*count_interval_high) for k in range(49)]

print('Done')

print('Calculating Magnetometer Indexes for Low: ', end='')
for i in range(len(EPOCH_low)):
    MAG_matrix_index_low[i] = [take_closest(MAG_low_time_computed, EPOCH_low[i] + k*count_interval_low ) for k in range(49)]


print('Done')


# ------------
# PITCH ANGLES
# ------------


print('Begin Calculating True Pitch Angles:', end='')


# ---------------------------------------------------
# ########UNIT VECTOR  --- UNIT IN MAG FRAME#########
# ---------------------------------------------------
#Calculate all the pitch angles
Pitch_output_high = np.ndarray(shape=(len(EPOCH_High),21,49), dtype='float64')
Pitch_output_low = np.ndarray(shape=(len(EPOCH_low),21,49), dtype='float64')

# X,Y,Z unit vectors for the EEPAA frame of reference
unit_vect = np.zeros([21,3])

# Create the unit vector
for i in range(len(pitch)):
    unit_vect[i][0] = math.cos((pitch[i]-90) * math.pi/180)
    unit_vect[i][1] = math.sin((pitch[i]-90)*math.pi/180)
    unit_vect[i][2] = 0


#High
for time_index in range(len(EPOCH_High)):
    for pitch_angle in range(21):
        Pitch_output_high[time_index][pitch_angle] = [(180/math.pi)*math.acos((np.dot(unit_vect[pitch_angle],mag_high[    int(MAG_matrix_index_high[time_index][engy])     ])  )/(np.linalg.norm(mag_high[int(MAG_matrix_index_high[time_index][engy])]))    ) for engy in range(49)]
#Low
for time_index in range(len(EPOCH_low)):
    for pitch_angle in range(21):
        Pitch_output_low[time_index][pitch_angle] = [(180/math.pi)*math.acos(    ( np.dot(unit_vect[pitch_angle],mag_low[    int(MAG_matrix_index_low[time_index][engy])     ])  )/(np.linalg.norm(mag_low[int(MAG_matrix_index_low[time_index][engy])]))    ) for engy in range(49)]







# --------------------------------------------------
# OUTPUT THE 'TRUE' PITCH ANGLES TO THE WORKING FILE
# --------------------------------------------------
pitchNAM = 'Pitch_Angle_Actual'
# ----------------------------
# WRITE OUTPUT PITCH DATA:HIGH
# ----------------------------

tag = 6
vardata= Pitch_output_high

#Var Info
varinfo = EEPAA_file_high.varinq(zvars_high[tag])
varinfo['Variable'] = pitchNAM

#Var Attributes
varattrs_high = EEPAA_file_high.varattsget(zvars_high[tag], expand=True)
varattrs_high['CATDESC'] = pitchNAM
varattrs_high['FIELDNAM'] = 'Pitch_Angle'
varattrs_high['VALIDMIN'] = [fillval_high,'CDF_REAL4'] # THERE APPEARS TO BE SOME ISSUE WITH THIS WHEN TRYING TO READ ITS ATTRIBUTES
varattrs_high['VALIDMAX'] = [180.0,'CDF_REAL4']
varattrs_high['SCALETYP'] = 'linear'
varattrs_high['UNITS'] = 'deg'
varattrs_high['LABLAXIS'] = pitchNAM
varattrs_high['VAR_TYPE'] = ['data','CDF_CHAR']
varattrs_high['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_high['LABL_PTR_2'] = 'LABL_PTR_2'
cdf_Pitch_Actual_high.write_var(varinfo,var_attrs=varattrs_high,var_data=vardata)




# ---------------------------
# WRITE OUTPUT PITCH DATA:LOW
# ---------------------------
pitchNAM = 'Pitch_Angle_Actual'
tag = 6
vardata= Pitch_output_low

#Var Info
varinfo = EEPAA_file_low.varinq(zvars_low[tag])
varinfo['Variable'] = pitchNAM

#Var Attributes
varattrs_low = EEPAA_file_low.varattsget(zvars_low[tag], expand=True)
varattrs_low['CATDESC'] = pitchNAM
varattrs_low['FIELDNAM'] = 'Pitch_Angle'
varattrs_low['VALIDMIN'] = [fillval_low,'CDF_REAL4']
varattrs_low['VALIDMAX'] = [180.0,'CDF_REAL4']
varattrs_low['SCALETYP'] = 'linear'
varattrs_low['UNITS'] = 'deg'
varattrs_low['LABLAXIS'] = pitchNAM
varattrs_low['VAR_TYPE'] = ['data','CDF_CHAR']
varattrs_low['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_low['LABL_PTR_2'] = 'LABL_PTR_2'
cdf_Pitch_Actual_low.write_var(varinfo,var_attrs=varattrs_low,var_data=vardata)





# Close the files
Mag_file_high.close()
Mag_file_low.close()
EEPAA_file_low.close()
EEPAA_file_high.close()
cdf_Pitch_Actual_high.close()
cdf_Pitch_Actual_low.close()



print('Done\n')




print("--- %s seconds for Magnetometer_rotations---" % (time.time() - start_time),'\n')