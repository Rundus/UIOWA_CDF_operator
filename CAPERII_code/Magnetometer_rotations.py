# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ---- TRICE Magnetometer Handling ----
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np,math,time
from cdflib import cdfwrite

start_time = time.time()
print('importing variables: ', end='')





# ------------------------------------
# Import Files & Functions & Variables
# ------------------------------------
from functions import take_closest
from files import EEPAA_file_1,EEPAA_file_2,root,user_path
from Variables import EEPAA_file_1_info,EEPAA_file_2_info,zvars_1,zvars_2,epoch_1,epoch_2,mag_dat,epoch_mag,pitches_1,ranges_1,ranges_2,count_interval_1_ns,count_interval_2_ns,fillval_1,fillval_2,pitchNAM

# ------------
# Output Files
# ------------
pitch_1_tag = 'pitch_actual_1'
pitch_2_tag = 'pitch_actual_2'

cdf_Pitch_Actual_1 = cdfwrite.CDF(user_path + root + "/" + pitch_1_tag + '.cdf',cdf_spec=EEPAA_file_1_info,delete=True)
cdf_Pitch_Actual_2 = cdfwrite.CDF(user_path + root + "/" + pitch_2_tag + '.cdf',cdf_spec=EEPAA_file_2_info,delete=True)
print('Done')


# --------------------------------------------------------------------------------------------------------------------
# Find the Magnetometer values which best match the times of the EEPAA data and create reduced mag_data for each eepaa
# --------------------------------------------------------------------------------------------------------------------

MAG_matrix_index_1 = np.zeros(shape= (len(epoch_1),49))
MAG_matrix_index_2 = np.zeros(shape =(len(epoch_2),49))


# ASSIGN THE MAGNETIC FIELD VALUES

print('Calculating closest match Magnetometer Indexes for eepaa1: ', end='')

for i in ranges_1[0]:
    MAG_matrix_index_1[i] = [ take_closest(epoch_mag, epoch_1[i] + (48 - k)*count_interval_1_ns[i] ) for k in range(49)]

print('Done')

print('Calculating closest match Magnetometer Indexes for eepaa2: ', end='')
for i in ranges_2[0]:
    MAG_matrix_index_2[i] = [ take_closest(epoch_mag, epoch_2[i] + (48 - k)*count_interval_2_ns[i] ) for k in range(49)]

print('Done')



# ----------------------------------------------------------------------------------------------------------
# -----Convert the Magnetometer Coordinates to the detector coordinates AND calculate True Pitch Angle -----
# ----------------------------------------------------------------------------------------------------------

# ---------------------------------------------------
# ########UNIT VECTOR  --- UNIT IN MAG FRAME#########
# ---------------------------------------------------
print('Rotating Matrices: ',end='')

# ------
# EEPAA1
# ------
unit_1 = np.zeros(shape=(21,3))

for i in ranges_1[1]:
    unit_1[i][0] = 0
    unit_1[i][1] = np.sin(  math.radians( pitches_1[i]))
    unit_1[i][2] = - np.cos(math.radians(pitches_1[i])  )



rot = -45
Theta = math.radians(rot)
R_z = np.array([ [math.cos(Theta),-math.sin(Theta),0],[math.sin(Theta),math.cos(Theta),0],[0,0,1]])
mag_rot_1 = [np.matmul(R_z,unit_1[i]) for i in ranges_1[1]]



# ------
# EEPAA2
# ------

unit_2 = np.zeros(shape=(21,3))

for i in ranges_2[1]:
    unit_2[i][0] = 0
    unit_2[i][1] = np.sin(  math.radians( pitches_1[i]))
    unit_2[i][2] = - np.cos(math.radians(pitches_1[i])  )

rot = 135
Theta = math.radians(rot)
R_z = np.array([ [math.cos(Theta),-math.sin(Theta),0],[math.sin(Theta),math.cos(Theta),0],[0,0,1]])
mag_rot_2 = [np.matmul(R_z,unit_2[i]) for i in ranges_2[1]]

print('Done')



# -------------------------------
# CALCULATE THE TRUE PITCH ANGLES
# -------------------------------
print('Begin Calculating True Pitch Angles for eepaa1: ', end='')


#Calculate all the pitch angles
Pitch_output_1 = np.ndarray(shape=(len(epoch_1),21,49), dtype='float64')
Pitch_output_2 = np.ndarray(shape=(len(epoch_1),21,49), dtype='float64')


# #eepaa1
for time_index in ranges_1[0]:
    for pitch_angle in ranges_1[1]:
        Pitch_output_1[time_index][pitch_angle] = [math.degrees(math.acos((np.dot(unit_1[pitch_angle],mag_dat[    int(MAG_matrix_index_1[time_index][engy])     ])  )/(np.linalg.norm(mag_dat[int(MAG_matrix_index_1[time_index][engy])]))    )) for engy in ranges_1[2]]

print('Done')

print('Begin Calculating True Pitch Angles for eepaa2: ', end='')
#eepaa2
for time_index in ranges_2[0]:
    for pitch_angle in ranges_2[1]:
        Pitch_output_2[time_index][pitch_angle] = [math.degrees(math.acos(    ( np.dot(unit_2[pitch_angle],mag_dat[    int(MAG_matrix_index_2[time_index][engy])     ])  )/(np.linalg.norm(mag_dat[int(MAG_matrix_index_2[time_index][engy])]))    )) for engy in ranges_2[2]]




# --------------------------------------------------
# OUTPUT THE 'TRUE' PITCH ANGLES TO THE WORKING FILE
# --------------------------------------------------



# ------------------------------
# WRITE OUTPUT PITCH DATA:eepaa1
# ------------------------------

tag = 18
vardata= Pitch_output_1

#Var Info
varinfo = EEPAA_file_1.varinq(zvars_1[tag])
varinfo['Variable'] = pitchNAM

#Var Attributes
varattrs_high = EEPAA_file_1.varattsget(zvars_1[tag], expand=True)
varattrs_high['CATDESC'] = [pitchNAM,'CDF_CHAR']
varattrs_high['FIELDNAM'] = [pitchNAM,'CDF_CHAR']
varattrs_high['VALIDMIN'] = [fillval_1[0],'CDF_FLOAT'] # THERE APPEARS TO BE SOME ISSUE WITH THIS WHEN TRYING TO READ ITS ATTRIBUTES
varattrs_high['VALIDMAX'] = [180.0,'CDF_REAL4']
varattrs_high['SCALETYP'] = ['linear','CDF_CHAR']
varattrs_high['UNITS'] = ['deg','CDF_CHAR']
varattrs_high['LABLAXIS'] = [pitchNAM,'CDF_CHAR']
varattrs_high['VAR_TYPE'] = ['data','CDF_CHAR']
varattrs_high['LABL_PTR_1'] = ['LABL_PTR_1','CDF_CHAR']
varattrs_high['LABL_PTR_2'] = ['LABL_PTR_2','CDF_CHAR']
cdf_Pitch_Actual_1.write_var(varinfo,var_attrs=varattrs_high,var_data=vardata)


# ------------------------------
# WRITE OUTPUT PITCH DATA:eepaa2
# ------------------------------
vardata= Pitch_output_2

#Var Info
varinfo = EEPAA_file_2.varinq(zvars_2[tag])
varinfo['Variable'] = pitchNAM

#Var Attributes
varattrs_low = EEPAA_file_2.varattsget(zvars_2[tag], expand=True)
varattrs_low['CATDESC'] = [pitchNAM,'CDF_CHAR']
varattrs_low['FIELDNAM'] = [pitchNAM,'CDF_CHAR']
varattrs_low['VALIDMIN'] = [fillval_2[0],'CDF_REAL4']
varattrs_low['VALIDMAX'] = [180.0,'CDF_REAL4']
varattrs_low['SCALETYP'] = ['linear','CDF_CHAR']
varattrs_low['UNITS'] = ['deg','CDF_CHAR']
varattrs_low['LABLAXIS'] = [pitchNAM,'CDF_CHAR']
varattrs_low['VAR_TYPE'] = ['data','CDF_CHAR']
varattrs_low['LABL_PTR_1'] = ['LABL_PTR_1','CDF_CHAR']
varattrs_low['LABL_PTR_2'] = ['LABL_PTR_2','CDF_CHAR']

cdf_Pitch_Actual_2.write_var(varinfo,var_attrs=varattrs_low,var_data=vardata)









# Close the files
cdf_Pitch_Actual_1.close()
cdf_Pitch_Actual_2.close()







print('Done\n')
print("--- %s seconds for Magnetometer_rotations---" % (time.time() - start_time),'\n')