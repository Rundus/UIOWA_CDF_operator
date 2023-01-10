# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE Final output file code--------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import time
from cdflib import cdfwrite
from cdflib import cdfread

start_time = time.time()
print('importing ALL variables: ', end='')


# -----------------------------
# GET THE VARIOUS FILE ELEMENTS
# -----------------------------


from files import EEPAA_file_high,EEPAA_file_low,Flux_aligned_high,Flux_aligned_low,pitch_file_actual_high,pitch_file_actual_low,J_par_file_high,J_par_file_low,dist_file_high,dist_file_low,cal_file_counts_high,cal_file_counts_low
from files import output_path,root,cal_diff_file_high,cal_diff_file_low
# from Variables import pitch_actual_high_info,pitch_actual_low_info,Flux_aligned_low_info,Flux_aligned_high_info
from Variables import zvars_low,zvars_high,zvars_pitch_high,zvars_pitch_low,zvars_flux_high,zvars_flux_low,zvars_J_par_high,J_par_high_info
from Variables import J_par_low_info,zvars_J_par_low,zvars_cal_diff_high,zvars_cal_diff_low
from Variables import dist_high,dist_low,dist_high_info,dist_low_info,zvars_high_dist,zvars_low_dist
from Variables import EEPAA_file_low_info,EEPAA_file_high_info,cal_counts_high_info,cal_counts_low_info,zvars_cal_counts_high,zvars_cal_counts_low



# Get the working global attributes
globalAttrs_high = EEPAA_file_high.globalattsget(expand=True)
globalAttrs_low = EEPAA_file_low.globalattsget(expand=True)



# ------------------------------
# WRITING THE FINAL OUTPUT FILES
# ------------------------------

# Output Files
cdf_output_high = cdfwrite.CDF(root + output_path + 'TRICE_52003_l2_eepaa_20181208T082239_processed&_v1.1.2.cdf',cdf_spec=EEPAA_file_high_info, delete =True)
cdf_output_low = cdfwrite.CDF(root + output_path + 'TRICE_52004_l2_eepaa_20181208T082243_processed&J_v1.1.2.cdf',cdf_spec=EEPAA_file_low_info, delete=True)

# Write global attributes
cdf_output_high.write_globalattrs(globalAttrs_high)
cdf_output_low.write_globalattrs(globalAttrs_low)

print('Done')




# -------------------------------------------------------------------------
# Copy all relevant variable data from original input files to output files
# -------------------------------------------------------------------------


# ----------High-----------
print('Copying original input file for HIGH: ', end='')

# Assigning the High data
for num_var in range(0,len(zvars_high)):
    if zvars_high[num_var] == 'eepaa_LABL_1' :
        Label_1 = ['-10 deg', '00 deg', '10 deg', '20 deg', '30 deg', '40 deg', '50 deg',
                  '60 deg', '70 deg', '80 deg', '90 deg', '100 deg', '110 deg', '120 deg',
                  '130 deg', '140 deg', '150 deg', '160 deg', '170 deg', '180 deg', '190 deg']
        vardata = Label_1

    elif zvars_high[num_var] =='eepaa_LABL_2':
        Label_2 =['11486 eV            ', '9445 eV             ', '8252 eV             ',
  '7405 eV             ', '6748 eV             ', '6211 eV             ',
  '5758 eV             ', '5365 eV             ', '5018 eV             ',
  '4708 eV             ', '4427 eV             ', '4171 eV             ',
  '3935 eV             ', '3718 eV             ', '3515 eV             ',
  '3325 eV             ', '3146 eV             ', '2978 eV             ',
  '2818 eV             ', '2668 eV             ', '2524 eV             ',
  '2387 eV             ', '2256 eV             ', '2131 eV             ',
  '2010 eV             ', '1895 eV             ', '1784 eV             ',
  '1677 eV             ', '1574 eV             ', '1474 eV             ',
  '1377 eV             ', '1284 eV             ', '1193 eV             ',
  '1106 eV             ', '1020 eV             ', '938 eV              ',
  '875 eV              ', '800 eV              ', '700 eV              ',
  '625 eV              ', '550 eV              ', '500 eV              ',
  '485 eV              ', '450 eV              ', '400 eV              ',
  '350 eV              ', '300 eV              ', '175 eV              ',
  '60 eV               ']
        vardata= Label_2

    else:
        vardata = EEPAA_file_high.varget(zvars_high[num_var])

    varinfo = EEPAA_file_high.varinq(zvars_high[num_var])                         # Get variable specification
    varattrs_high = EEPAA_file_high.varattsget(zvars_high[num_var], expand=True)  # Get the variable's attributes
    cdf_output_high.write_var(varinfo,var_attrs=varattrs_high,var_data=vardata)

print('Done')


# ---------- Low ----------
print('Copying original input file for LOW: ', end='')

# Assigning the High data
for num_var in range(0,len(zvars_low)):
    if zvars_low[num_var] == 'eepaa_LABL_1' :
        Label_1 = ['-10 deg', '00 deg', '10 deg', '20 deg', '30 deg', '40 deg', '50 deg',
                  '60 deg', '70 deg', '80 deg', '90 deg', '100 deg', '110 deg', '120 deg',
                  '130 deg', '140 deg', '150 deg', '160 deg', '170 deg', '180 deg', '190 deg']
        vardata = Label_1

    elif zvars_low[num_var] =='eepaa_LABL_2':
        Label_2 =['11486 eV            ', '9445 eV             ', '8252 eV             ',
  '7405 eV             ', '6748 eV             ', '6211 eV             ',
  '5758 eV             ', '5365 eV             ', '5018 eV             ',
  '4708 eV             ', '4427 eV             ', '4171 eV             ',
  '3935 eV             ', '3718 eV             ', '3515 eV             ',
  '3325 eV             ', '3146 eV             ', '2978 eV             ',
  '2818 eV             ', '2668 eV             ', '2524 eV             ',
  '2387 eV             ', '2256 eV             ', '2131 eV             ',
  '2010 eV             ', '1895 eV             ', '1784 eV             ',
  '1677 eV             ', '1574 eV             ', '1474 eV             ',
  '1377 eV             ', '1284 eV             ', '1193 eV             ',
  '1106 eV             ', '1020 eV             ', '938 eV              ',
  '875 eV              ', '800 eV              ', '700 eV              ',
  '625 eV              ', '550 eV              ', '500 eV              ',
  '485 eV              ', '450 eV              ', '400 eV              ',
  '350 eV              ', '300 eV              ', '175 eV              ',
  '60 eV               ']
        vardata= Label_2
    else:
        vardata = EEPAA_file_low.varget(zvars_low[num_var])

    varinfo = EEPAA_file_low.varinq(zvars_low[num_var])                        # Get variable specification
    varattrs_low = EEPAA_file_low.varattsget(zvars_low[num_var], expand=True)  # Get the variable's attributes
    cdf_output_low.write_var(varinfo,var_attrs=varattrs_low,var_data=vardata)

print('Done')


# # ------------------------------------------
# # Write out all the newly computed variables
# # ------------------------------------------
#
#
# # ----------------------------
# # WRITE OUTPUT PITCH DATA:HIGH
# # ----------------------------
# print('Writing out pitch actual data: ', end='')
#
# vardata = pitch_file_actual_high.varget('Pitch_Angle_Actual')
# varinfo = pitch_file_actual_high.varinq(zvars_pitch_high[0])
# varattrs = pitch_file_actual_high.varattsget(zvars_pitch_low[0], expand=True)
# cdf_output_high.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
#
#
# # ---------------------------
# # WRITE OUTPUT PITCH DATA:LOW
# # ---------------------------
#
# vardata = pitch_file_actual_low.varget('Pitch_Angle_Actual')
# varinfo = pitch_file_actual_low.varinq(zvars_pitch_low[0])
# varattrs = pitch_file_actual_low.varattsget(zvars_pitch_low[0], expand=True)
# cdf_output_low.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
#
# print('Done')
#
# # -----------------------------------
# # WRITE OUTPUT ALIGNED FLUX DATA:HIGH
# # -----------------------------------
#
# print('Writing out flux data: ', end='')
#
# vardata = Flux_aligned_high.varget('Differential_Energy_Flux_Aligned')
# varinfo = Flux_aligned_high.varinq(zvars_flux_high[0])
# varattrs = Flux_aligned_high.varattsget(zvars_flux_high[0], expand=True)
# cdf_output_high.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
# vardata = Flux_aligned_high.varget('Differential_Number_Flux_Aligned')
# varinfo = Flux_aligned_high.varinq(zvars_flux_high[1])
# varattrs = Flux_aligned_high.varattsget(zvars_flux_high[1], expand=True)
# cdf_output_high.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
# # ----------------------------------
# # WRITE OUTPUT ALIGNED FLUX DATA:LOW
# # ----------------------------------
#
# vardata = Flux_aligned_low.varget('Differential_Energy_Flux_Aligned')
# varinfo = Flux_aligned_low.varinq(zvars_flux_low[0])
# varattrs = Flux_aligned_low.varattsget(zvars_flux_low[0], expand=True)
# cdf_output_low.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
# vardata = Flux_aligned_low.varget('Differential_Number_Flux_Aligned')
# varinfo = Flux_aligned_low.varinq(zvars_flux_low[1])
# varattrs = Flux_aligned_low.varattsget(zvars_flux_low[1], expand=True)
# cdf_output_low.write_var(varinfo,var_attrs=varattrs,var_data=vardata)




# print('Done')
#
#
# # ------------------------------
# # WRITE OUTPUT DISTRIBUTION:HIGH
# # ------------------------------
# print('Writing out Distribution Function data: ', end='')
#
# vardata = dist_file_high.varget('Distribution_Function')
# varinfo = dist_file_high.varinq(zvars_high_dist[0])
# varattrs = dist_file_high.varattsget(zvars_high_dist[0], expand=True)
# cdf_output_high.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
#
# # -----------------------------
# # WRITE OUTPUT DISTRIBUTION:LOW
# # -----------------------------
#
# vardata = dist_file_low.varget('Distribution_Function')
# varinfo = dist_file_low.varinq(zvars_low_dist[0])
# varattrs = dist_file_low.varattsget(zvars_low_dist[0], expand=True)
# cdf_output_low.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
#
# print('Done')
# print('Writing out J_parallel data: ', end='')
# # ----------------------------
# # WRITE OUTPUT J_PARALLEL:HIGH
# # ----------------------------
#
# vardata = J_par_file_high.varget('Parallel_Current_Density')
# varinfo = J_par_file_high.varinq(zvars_J_par_high[0])
# varattrs = J_par_file_high.varattsget(zvars_J_par_high[0], expand=True)
# cdf_output_high.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
#
# # ---------------------------
# # WRITE OUTPUT J_PARALLEL:LOW
# # ---------------------------
#
# vardata = J_par_file_low.varget('Parallel_Current_Density')
# varinfo = J_par_file_low.varinq(zvars_J_par_low[0])
# varattrs = J_par_file_low.varattsget(zvars_J_par_low[0], expand=True)
# cdf_output_low.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
#
# # -----------------------------
# # WRITE OUTPUT CAL COUNTS: HIGH
# # -----------------------------
#
# vardata = cal_file_counts_high.varget('Calibrated_Counts')
# varinfo = cal_file_counts_high.varinq(zvars_cal_counts_high[0])
# varattrs = cal_file_counts_high.varattsget(zvars_cal_counts_high[0], expand=True)
# cdf_output_high.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
# # ----------------------------
# # WRITE OUTPUT CAL COUNTS: LOW
# # ----------------------------
#
# vardata = cal_file_counts_low.varget('Calibrated_Counts')
# varinfo = cal_file_counts_low.varinq(zvars_cal_counts_low[0])
# varattrs = cal_file_counts_low.varattsget(zvars_cal_counts_low[0], expand=True)
# cdf_output_low.write_var(varinfo,var_attrs= varattrs,var_data=vardata)



# -------------------------------------------
# WRITE OUTPUT PROCESSED DIFFN and DIFFE: LOW
# -------------------------------------------
vardata = cal_diff_file_low.varget('Differential_Energy_Flux')
varinfo = cal_diff_file_low.varinq(zvars_cal_diff_low[0])
varinfo['Variable'] = 'Differential_Energy_Flux_noise_removed'
varattrs = cal_diff_file_low.varattsget(zvars_cal_diff_low[0], expand=True)
varattrs['CATDESC'] = 'Differential_Energy_Flux_noise_removed'
varattrs['FIELDNAM'] = 'Differential_Energy_Flux_noise_removed'
cdf_output_low.write_var(varinfo,var_attrs=varattrs,var_data=vardata)

vardata = cal_diff_file_low.varget('Differential_Number_Flux')
varinfo = cal_diff_file_low.varinq(zvars_cal_diff_low[0])
varinfo['Variable'] = 'Differential_Number_Flux_noise_removed'
varattrs = cal_diff_file_low.varattsget(zvars_cal_diff_low[0], expand=True)
varattrs['CATDESC'] = 'Differential_Energy_Flux_noise_removed'
varattrs['FIELDNAM'] = 'Differential_Energy_Flux_noise_removed'
cdf_output_low.write_var(varinfo,var_attrs=varattrs,var_data=vardata)


# --------------------------------------------
# WRITE OUTPUT PROCESSED DIFFN and DIFFE: HIGH
# --------------------------------------------
vardata = cal_diff_file_high.varget('Differential_Energy_Flux')
varinfo = cal_diff_file_high.varinq(zvars_cal_diff_high[0])
varinfo['Variable'] = 'Differential_Energy_Flux_noise_removed'
varattrs = cal_diff_file_high.varattsget(zvars_cal_diff_high[0], expand=True)
varattrs['CATDESC'] = 'Differential_Energy_Flux_noise_removed'
varattrs['FIELDNAM'] = 'Differential_Energy_Flux_noise_removed'
cdf_output_high.write_var(varinfo,var_attrs=varattrs,var_data=vardata)

vardata = cal_diff_file_high.varget('Differential_Number_Flux')
varinfo = cal_diff_file_high.varinq(zvars_cal_diff_high[0])
varinfo['Variable'] = 'Differential_Number_Flux_noise_removed'
varattrs = cal_diff_file_high.varattsget(zvars_cal_diff_high[0], expand=True)
varattrs['CATDESC'] = 'Differential_Energy_Flux_noise_removed'
varattrs['FIELDNAM'] = 'Differential_Energy_Flux_noise_removed'
cdf_output_high.write_var(varinfo,var_attrs=varattrs,var_data=vardata)



print('Done')







# Close all files
EEPAA_file_high.close()
EEPAA_file_low.close()
pitch_file_actual_high.close()
pitch_file_actual_low.close()
cdf_output_high.close()
cdf_output_low.close()
J_par_file_high.close()
J_par_file_low.close()
cal_file_counts_high.close()
cal_file_counts_low.close()


print("--- %s seconds for program---" % (time.time() - start_time) ,'\n')

print('END PROGRAM')