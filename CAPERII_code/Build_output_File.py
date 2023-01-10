# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -------- CAPERII Final output file code-------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import time, numpy as np
from cdflib import cdfwrite

start_time = time.time()
print('importing ALL variables: ', end='')


# -----------------------------
# GET THE VARIOUS FILE ELEMENTS
# -----------------------------

from files import EEPAA_file_1,EEPAA_file_2, output_root,user_path
# from files import diffFlux_file_2,pitch_actual_file_1,pitch_actual_file_2,flux_aligned_file_1,flux_aligned_file_2,dist_file_1,dist_file_2
from Variables import EEPAA_file_1_info as EEPAA_working_1_info, EEPAA_file_2_info as EEPAA_working_2_info,zvars_1,zvars_2,globalAttrs_1,globalAttrs_2
from Variables import fillval_1,fillval_2
# from Variables import zvars_pitch_2,pitchNAM,diffFlux_1_info,diffFlux_2_info,flux_aligned_1_info,flux_aligned_2_info
# from Variables import zvars_flux_aligned_1,zvars_flux_aligned_2,zvars_dflux_1,zvars_dflux_2,dist_1_info,dist_2_info,zvars_dist_1,zvars_dist_2



# ------------------------------
# WRITING THE FINAL OUTPUT FILES
# ------------------------------

# Output Files
cdf_output_1 = cdfwrite.CDF(user_path+ output_root +  '/' + '/Output/CAPERII_eepaa1_pitch&dist.cdf',cdf_spec=EEPAA_working_1_info, delete =True)
cdf_output_2 = cdfwrite.CDF(user_path + output_root +  '/' + '/Output/CAPERII_eepaa2_pitch&dist.cdf',cdf_spec=EEPAA_working_2_info, delete =True)

# Write global attributes
cdf_output_1.write_globalattrs(globalAttrs_1)
cdf_output_2.write_globalattrs(globalAttrs_2)

print('Done')


# ----------------------------------------------------------------------------------------------------------------------------------------------------
#Compute the adjusted epoch values and load them into this file to replace the normal epoch files for the wave-group comparison to Spencer's Wave data
# ----------------------------------------------------------------------------------------------------------------------------------------------------
# from Wavegroup_Processing.epoch_adjusting import epoch_adjusted_1, epoch_adjusted_2


# -------------------------------------------------------------------------
# Copy all relevant variable data from original input files to output files
# -------------------------------------------------------------------------


# ----------eepaa1-----------
print('Copying original input file for eepaa1: ', end='')

# Assigning the High data
for num_var in range(0,len(zvars_1)):
    if zvars_1[num_var] == 'eepaa_LABL_1' :
        Label_1 = ['-10 deg', '00 deg', '10 deg', '20 deg', '30 deg', '40 deg', '50 deg',
                  '60 deg', '70 deg', '80 deg', '90 deg', '100 deg', '110 deg', '120 deg',
                  '130 deg', '140 deg', '150 deg', '160 deg', '170 deg', '180 deg', '190 deg']
        vardata = Label_1

    elif zvars_1[num_var] =='eepaa_LABL_2':
        Label_2 = ['12600 eV            ', '10350 eV            ', '9050 eV             ',
                   '8100 eV             ', '7400 eV             ', '6800 eV             ',
                   '6300 eV             ', '5900 eV             ', '5500 eV             ',
                   '5150 eV             ', '4850 eV             ', '4550 eV             ',
                   '4300 eV             ', '4100 eV             ', '3850 eV             ',
                   '3650 eV             ', '3450 eV             ', '3250 eV             ',
                   '3100 eV             ', '2900 eV             ', '2800 eV             ',
                   '2600 eV             ', '2500 eV             ', '2350 eV             ',
                   '2200 eV             ', '2100 eV             ', '1950 eV             ',
                   '1850 eV             ', '1750 eV             ', '1600 eV             ',
                   '1500 eV             ', '1400 eV             ', '1300 eV             ',
                   '1200 eV             ', '1100 eV             ', '1050 eV             ',
                   '950 eV              ', '850 eV              ', '750 eV              ',
                   '700 eV              ', '650 eV              ', '600 eV              ',
                   '550 eV              ', '500 eV              ', '450 eV              ',
                   '400 eV              ', '350 eV              ', '300 eV              ',
                   '250 eV              ']
        vardata= Label_2
    elif zvars_1[num_var] in ['BIAS_DAC','minor_frame_counter']:
        #Because no records of these variables exist, a fillvalue is replacing them
        vardata= np.zeros(1)
        vardata[0] = fillval_1

    #Use the adjusted epoch values
    # elif zvars_1[num_var] == 'Epoch':
    #     vardata = epoch_adjusted_1

    else:
        vardata = EEPAA_file_1.varget(zvars_1[num_var])

    varinfo = EEPAA_file_1.varinq(zvars_1[num_var])                         # Get variable specification
    varattrs_1 = EEPAA_file_1.varattsget(zvars_1[num_var], expand=True)  # Get the variable's attributes
    cdf_output_1.write_var(varinfo,var_attrs=varattrs_1,var_data=vardata)

print('Done')


# ---------- eepaa2 ----------
print('Copying original input file for eepaa2: ', end='')

# Assigning the High data
for num_var in range(0,len(zvars_2)):
    if zvars_2[num_var] == 'eepaa_LABL_1' :
        Label_1 = ['-10 deg', '00 deg', '10 deg', '20 deg', '30 deg', '40 deg', '50 deg',
                  '60 deg', '70 deg', '80 deg', '90 deg', '100 deg', '110 deg', '120 deg',
                  '130 deg', '140 deg', '150 deg', '160 deg', '170 deg', '180 deg', '190 deg']
        vardata = Label_1

    elif zvars_2[num_var] =='eepaa_LABL_2':
        Label_2 = ['12600 eV            ', '10350 eV            ', '9050 eV             ',
                   '8100 eV             ', '7400 eV             ', '6800 eV             ',
                   '6300 eV             ', '5900 eV             ', '5500 eV             ',
                   '5150 eV             ', '4850 eV             ', '4550 eV             ',
                   '4300 eV             ', '4100 eV             ', '3850 eV             ',
                   '3650 eV             ', '3450 eV             ', '3250 eV             ',
                   '3100 eV             ', '2900 eV             ', '2800 eV             ',
                   '2600 eV             ', '2500 eV             ', '2350 eV             ',
                   '2200 eV             ', '2100 eV             ', '1950 eV             ',
                   '1850 eV             ', '1750 eV             ', '1600 eV             ',
                   '1500 eV             ', '1400 eV             ', '1300 eV             ',
                   '1200 eV             ', '1100 eV             ', '1050 eV             ',
                   '950 eV              ', '850 eV              ', '750 eV              ',
                   '700 eV              ', '650 eV              ', '600 eV              ',
                   '550 eV              ', '500 eV              ', '450 eV              ',
                   '400 eV              ', '350 eV              ', '300 eV              ',
                   '250 eV              ']
        vardata= Label_2
    elif zvars_2[num_var] in ['BIAS_DAC','minor_frame_counter']:
        #Because no records of these variables exist, a fillvalue is replacing them
        vardata= np.zeros(1)
        vardata[0] = fillval_2
    # Use the adjusted epoch values
    # elif zvars_1[num_var] == 'Epoch':
    #     vardata = epoch_adjusted_2

    else:
        vardata = EEPAA_file_2.varget(zvars_2[num_var])

    varinfo = EEPAA_file_2.varinq(zvars_2[num_var])                        # Get variable specification
    varattrs_2 = EEPAA_file_2.varattsget(zvars_2[num_var], expand=True)  # Get the variable's attributes
    cdf_output_2.write_var(varinfo,var_attrs=varattrs_2,var_data=vardata)

print('Done')




#
# # ------------------------------------------
# # Write out all the newly computed variables
# # ------------------------------------------
#
#
# # -----------------------------
# # WRITE OUTPUT FLUX DATA:eepaa1
# # -----------------------------
#
# print('Writing out flux data for eepaa1: ', end='')
#
# # DiffEFlux
# vardata = diffFlux_file_1.varget('Differential_Energy_Flux')
# varinfo = diffFlux_file_1.varinq(zvars_dflux_1[1])
# varattrs = diffFlux_file_1.varattsget(zvars_dflux_1[1], expand=True)
# cdf_output_1.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
#
# # DiffNFlux
# vardata = diffFlux_file_1.varget('Differential_Number_Flux')
# varinfo = diffFlux_file_1.varinq(zvars_dflux_1[0])
# varattrs = diffFlux_file_1.varattsget(zvars_dflux_1[0], expand=True)
# cdf_output_1.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
#
# print('Done')
#
# # -----------------------------
# # WRITE OUTPUT FLUX DATA:eepaa2
# # -----------------------------
#
# print('Writing out flux data for eepaa2: ', end='')
#
# # DiffEFlux
# vardata = diffFlux_file_2.varget('Differential_Energy_Flux')
# varinfo = diffFlux_file_2.varinq(zvars_dflux_2[1])
# varattrs = diffFlux_file_2.varattsget(zvars_dflux_2[1], expand=True)
# cdf_output_2.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
#
# # DiffNFlux
# vardata = diffFlux_file_2.varget('Differential_Number_Flux')
# varinfo = diffFlux_file_2.varinq(zvars_dflux_2[0])
# varattrs = diffFlux_file_2.varattsget(zvars_dflux_2[0], expand=True)
# cdf_output_2.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
# print('Done')
#
#
#
#
#
#
# # ------------------------------
# # WRITE OUTPUT PITCH DATA:eepaa1
# # ------------------------------
# print('Writing out magnetometer pitch data: ', end='')
#
#
# vardata = pitch_actual_file_1.varget(pitchNAM)
# varinfo = pitch_actual_file_1.varinq(zvars_pitch_1[0])
# varattrs = pitch_actual_file_1.varattsget(zvars_pitch_1[0], expand=True)
# cdf_output_1.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
#
# # ------------------------------
# # WRITE OUTPUT PITCH DATA:eepaa2
# # ------------------------------
#
# vardata = pitch_actual_file_2.varget(pitchNAM)
# varinfo = pitch_actual_file_2.varinq(zvars_pitch_2[0])
# varattrs = pitch_actual_file_2.varattsget(zvars_pitch_2[0], expand=True)
# cdf_output_2.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
# print('Done')
#
#
#
# # -----------------------------
# # WRITE OUTPUT FLUX DATA:eepaa1
# # -----------------------------
#
# print('Writing out flux data: ', end='')
#
# vardata = flux_aligned_file_1.varget('Differential_Energy_Flux_Aligned')
# varinfo = flux_aligned_file_1.varinq(zvars_flux_aligned_1[0])
# varattrs = flux_aligned_file_1.varattsget(zvars_flux_aligned_1[0], expand=True)
# cdf_output_1.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
# vardata = flux_aligned_file_1.varget('Differential_Number_Flux_Aligned')
# varinfo = flux_aligned_file_1.varinq(zvars_flux_aligned_1[1])
# varattrs = flux_aligned_file_1.varattsget(zvars_flux_aligned_1[1], expand=True)
# cdf_output_1.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
# # -----------------------------
# # WRITE OUTPUT FLUX DATA:eepaa2
# # -----------------------------
#
# vardata = flux_aligned_file_2.varget('Differential_Energy_Flux_Aligned')
# varinfo = flux_aligned_file_2.varinq(zvars_flux_aligned_2[0])
# varattrs = flux_aligned_file_2.varattsget(zvars_flux_aligned_2[0], expand=True)
# cdf_output_2.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
# vardata = flux_aligned_file_2.varget('Differential_Number_Flux_Aligned')
# varinfo = flux_aligned_file_2.varinq(zvars_flux_aligned_2[1])
# varattrs = flux_aligned_file_2.varattsget(zvars_flux_aligned_2[1], expand=True)
# cdf_output_2.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
# print('Done')
#
#
# # ------------------------------
# # WRITE OUTPUT DISTRIBUTION:HIGH
# # ------------------------------
# print('Writing out Distribution Function data: ', end='')
#
# vardata = dist_file_1.varget('Distribution_Function')
# varinfo = dist_file_1.varinq(zvars_dist_1[0])
# varattrs = dist_file_1.varattsget(zvars_dist_1[0], expand=True)
# cdf_output_1.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
#
# # -----------------------------
# # WRITE OUTPUT DISTRIBUTION:LOW
# # -----------------------------
#
# vardata = dist_file_2.varget('Distribution_Function')
# varinfo = dist_file_2.varinq(zvars_dist_2[0])
# varattrs = dist_file_2.varattsget(zvars_dist_2[0], expand=True)
# cdf_output_2.write_var(varinfo,var_attrs=varattrs,var_data=vardata)
#
#
#
#
# print('Done')
#
#
# #--------------------------------------------------------------------------------------------------------------------------------------------
#
# # Close all files
# EEPAA_file_1.close()
# EEPAA_file_2.close()
# diffFlux_file_1.close()
# diffFlux_file_2.close()
# pitch_actual_file_1.close()
# pitch_actual_file_2.close()
# dist_file_1.close()
# dist_file_2.close()
cdf_output_1.close()
cdf_output_2.close()


print("--- %s seconds for program---" % (time.time() - start_time) ,'\n')

print('END PROGRAM')