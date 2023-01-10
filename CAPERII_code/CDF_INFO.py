# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GENERAL INFORMATION
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#   This will print a directory of all the main features of the CDF File. The Follow are the keys:

# 1) CDF -  the name and windows path of the CDF
# 2) Version  - the version of the CDF
# 3) Encoding - the endianness of the CDF
# 4) Majority - the row/column majority
# 5) zVariables - the dictionary for zVariable numbers and their corresponding names
# 6) rVariables - the dictionary for rVariable numbers and their corresponding names
# 7) Attributes - the dictionary for attribute numbers and their corresponding names and scopes
# 8) Checksum - the checksum indicator
# 9) Num_rdim - the number of dimensions, applicable only to rVariables
# 10) rDim_sizes - the dimensional sizes, applicable only to rVariables
# 11) Compressed - CDF is compressed at the file-level
# 12) LeapSecondUpdated - The last updated for the leap second table, if applicable


#------2. INQUIRE ABOUT ONE VARIABLE INFORMATION--------

#Note: Only works for zVariables
# print(file.varinq("Energy"))


#------3. INQUIRE ABOUT AN ATTRIBUTE--------
# print(file.attinq("Instrument_type"))

#------4. GET A VARIABLES ATTRIBUTE--------
#print(file.varattsget("Energy",expand=True))
# print(file.varattsget("Energy",expand=False))


#------5. GET ALL GLOBAL ATTRIBUTES--------
# print(file.globalattsget())


#------6. GET VARIABLE DATA--------

# print( file.varget("Differential_Energy_Flux",expand=True)  )
# print( file.varget("Pitch_Angle",expand=False)  )

#------7. EPOCHRANGE--------

# print(info["zVariables"])
# print(file.epochrange("Differential_Energy_Flux"))
# print(file.epochrange("Epoch"))

#------8. GET THE VERSION--------

# print(file.getVersion())




# ------- PRINT SUB-ATTRIBUTE INFORMATION -------

# Attriinquire = file.attinq(attribute=attnum)
# print(file.attinq(attribute=attnum))
#
# for i in range(0, Attriinquire['max_gr_entry']+1):
#     print("\nSub-attribute number " + str(i+1) + " is:")
#     print(file.attget(attribute=attnum, entry=i))

########################################################################################################################################################################################################



# -------------------------------------------------------------
# General File to Look at the Variables and Data of the program
# -------------------------------------------------------------

import cdflib
from cdflib import cdfread, cdfwrite
import itertools, numpy as np,time
#
# from files import mag_file,Mag_file_high,EEPAA_file_high
# from functions import printvarinfo
# from Variables import mag_file_info,zvars_mag,mag_dat,epoch_1,epoch_2,epoch_mag,zvars_high,count_interval_1,count_interval_2,count_interval_high,pitches_1,pitches_2
#


# ----------------------------------------------------------
# Print the variable information corresponding to indicies v
# ----------------------------------------------------------

v = [1,6,13]
# printvarinfo(file1,zvars1,v)
# printvarinfo(file2,zvars2,v)

# -----------------------------------------------------------
# Look at the Contents of the CAPERII EEPAA File and its data
# -----------------------------------------------------------
# print(count_interval_1,count_interval_2,count_interval_high)
# count_int = [count_interval_1[x]*1000 for x in range(len(count_interval_1))]
# print(count_interval_1[0:10])
# print(count_int[0:10])

# print(epoch_1[0])
# print(cdflib.cdfepoch.breakdown_tt2000(epoch_1[0],to_np = False))
# print(cdflib.cdfepoch.breakdown_tt2000(epoch_2[0],to_np = False))
# print(cdflib.cdfepoch.breakdown_tt2000(epoch_mag[0],to_np = False))




# ------------------------------------------------------------------
# Look at the Contents of the CAPERII Magnetometer File and its data
# ------------------------------------------------------------------

# print(mag_file_info)
# print(zvars_mag)
# print(mag_dat)
# print(np.transpose(mag_dat))
# print(len(epoch_1),len(epoch_2))
# print(len(epoch_mag))
# b = [3]
# printvarinfo(mag_file,zvars_mag,b)




# ----------------------------------------------------------------
# Look at the Contents of the TRICE Magnetometer File and its data
# ----------------------------------------------------------------

# mag_info = Mag_file_high.cdf_info()
# print(mag_info)


# ---------------------------------------------------------
# Look at the Contents of the TRICE EEPAA File and its data
# ---------------------------------------------------------
# placeholder = EEPAA_file_high.varattsget(zvars_high[1])
# placeholder1 = placeholder['FILLVAL']
# fillval_high = placeholder1[0]


# --------------------------------------------------------
# Look at the Contents of the Correlator File and its data
# --------------------------------------------------------
# from Variables import correlator_info,correlators,energy_correlator,phase_bins
#
# print(len(correlators),len(correlators[0]),len(correlators[0][0]))









