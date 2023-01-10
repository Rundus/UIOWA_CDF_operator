# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE REMOVE SUN SPIKES --------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This code gathers the data from the raw counts files and bins them into variables for file storage.
# The data is taken from the regions where the noise in the counts appears to be nominal to the
# regions where the main signal/phenomena is.

import itertools
import time
import matplotlib.pyplot as plt
print('importing variables: ', end='')
start_time = time.time()

# -----------
# IMPORT INFO
# -----------
import numpy as np
from functions import write_var_to_file
from files import counts_file_high,counts_file_low,root,attitude_control_file_high,attitude_control_file_low
from cdflib import cdfwrite
from Variables import counts_high_info,counts_low_info,Energies, Maskvals,rocket_names,zvars_counts_high,zvars_counts_low
from Variables import EPOCH_high_counts,EPOCH_low_counts,counts_high_raw,counts_low_raw,ranges_counts_high,ranges_counts_low,noise_threshholds
from Variables import t_noise_start1_high,t_noise_start2_high,t_noise_end1_high,t_noise_end2_high,t_noise_start1_low,t_noise_start2_low,t_noise_end1_low,t_noise_end2_low

print('Done')

# ------------
# OUTPUT FILES
# ------------
sun_spike_noise_file_high = cdfwrite.CDF(root + 'sun_spike_noise_dat_high',cdf_spec=counts_high_info,delete=True)
sun_spike_noise_file_low = cdfwrite.CDF(root + 'sun_spike_noise_dat_low',cdf_spec=counts_low_info,delete=True)

#Some Variables
ranges = [ranges_counts_low,ranges_counts_high]
raw_data = [counts_low_raw,counts_high_raw]
selects = [0,1]

# -------------------------------------------
# CODE TO CHARACTERIZE THE SPIN of the rocket
# -------------------------------------------
roll_output = [[dat for dat in EPOCH_low_counts], [dat for dat in EPOCH_high_counts]]
roll_output_epoch = [[dat for dat in EPOCH_low_counts], [dat for dat in EPOCH_high_counts]]
noise_data_high = []
noise_roll_high = []
noise_energy_high = []
noise_data_low = []
noise_roll_low = []
noise_energy_low = []
noise_total_data_low = []
noise_total_data_high = []

attitude_roll = [attitude_control_file_low.varget('Roll (Euler)'), attitude_control_file_high.varget('Roll (Euler)')]
attitude_epoch = [attitude_control_file_low.varget('Epoch'), attitude_control_file_high.varget('Epoch')]
epochs = [EPOCH_low_counts, EPOCH_high_counts]



for selector in selects:
    w_atti_roll = attitude_roll[selector]
    w_atti_roll_epoch = attitude_epoch[selector]
    w_epoch = epochs[selector]
    w_roll_output = roll_output[selector]
    w_roll_epoch_output = roll_output_epoch[selector]
    wattitude_epoch = np.array(attitude_epoch[selector])


    print('\n')
    for i in range(len(w_epoch)):
        if i%200 == 0:
            print('Matching the epoch/roll indices for ' + rocket_names[selector] + ': ' +  str(round(100 * (i/len(epochs[selector])),1)) + '%',end='\r')
        elif i == (len(epochs[selector])-1):
            print('Matching the epoch/roll indices for ' + rocket_names[selector] + ':' + ' 100%')

        index = np.abs(wattitude_epoch- w_epoch[i]).argmin()
        w_roll_output[i] = w_atti_roll[index]
        w_roll_epoch_output[i] = w_atti_roll_epoch[index]

    roll_output[selector] = w_roll_output
    roll_output_epoch[selector] = w_roll_epoch_output

#------------------
# OUTPUT NOISE DATA
#------------------
for select in selects:
    print('Preparing data for filtering ' + rocket_names[select],end='')
    which_range = ranges[select]
    raw_dat = raw_data[select]
    output_raw = np.zeros(shape=(len(raw_dat), len(raw_dat[0]), len(raw_dat[0][0])))
    output = np.zeros(shape=(len(raw_dat), len(raw_dat[0]), len(raw_dat[0][0])))
    wnoise_threshhold = noise_threshholds[selector]

    # Apply a Mask on the raw counts data
    for tme0,ptch0,engy0 in itertools.product(*which_range):
        if tme0 % 500 == 0:
            print('Preparing data for filtering ' + rocket_names[select] + ': ' + str(round(100 * (tme0 / len(epochs[select])), 1)) + '%', end='\r')
        elif tme0 == (len(epochs[select]) - 1) and ptch0 == 20 and engy0 == 48:
            print('Preparing data for filtering ' + rocket_names[select] + ':' + ' 100%')

        dat_checker = (raw_dat[tme0][ptch0][engy0] - Maskvals[select])
        if dat_checker <= 0:
            output_raw[tme0][ptch0][engy0] = 0
        else:
            output_raw[tme0][ptch0][engy0] = dat_checker

    #------
    # low #
    #------
    if select == 0:
        # APPLY THE LOW COUNTS MASK
        rollchecker = roll_output[select]
        for tme1, ptch1, engy1 in itertools.product(*which_range):

            if tme1%200 == 0:
                print('Filtering bad data for ' + rocket_names[select] + ': ' +  str(round(100 * (tme1/len(epochs[select])),1)) + '%',end='\r')
            elif tme1 == (len(epochs[select])-1) and ptch1 ==20 and engy1 == 48:
                print('Filtering bad data for ' + rocket_names[select] + ':' + ' 100%')

            if (tme1>= t_noise_start1_low and tme1<= t_noise_end1_low) or (tme1 >= t_noise_start2_low and tme1<= t_noise_end2_low):
                dat = raw_dat[tme1][ptch1][engy1]
                if dat >= wnoise_threshhold:
                    rchecker = float(rollchecker[tme1])
                    noise_total_data_low.append([dat,rchecker,Energies[engy1],ptch1])

        # -------------------------------------
        # Write out Noise Data
        # -------------------------------------
        vardata = np.array(output_raw, dtype='float64')
        varattributes = counts_file_low.varattsget(zvars_counts_low[18],expand=True)
        varattributes['CATDESC'] = 'counts_masked'
        varattributes['FIELDNAM'] = 'counts_masked'
        varattributes['VALIDMIN'] = [vardata.min(), 'CDF_FLOAT']
        varattributes['VALIDMAX'] = [vardata.max(), 'CDF_FLOAT']

        varinfo = counts_file_low.varinq(zvars_counts_low[18])
        varinfo['Variable'] = 'counts_masked'
        sun_spike_noise_file_low.write_var(varinfo, var_attrs=varattributes, var_data=vardata)


        vardata = np.array(roll_output[select],dtype='float64')
        attrs = ['Roll_Reduced', [-1.e+31], [vardata.min()], [vardata.max()], 'linear', 'deg', 'series']
        infos = [4, 44, len(vardata), attrs[0], [-9223372036854775807]]
        varattributes = {'CATDESC': attrs[0], 'DEPEND_0':'Roll_Epoch', 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0],'FILLVAL': np.array(attrs[1], dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0],'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2],'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP': attrs[4]}
        varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1],'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [],'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [],'Pad': np.array(infos[4], dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
        sun_spike_noise_file_low.write_var(varinfo, var_attrs=varattributes, var_data=vardata)


        vardata = np.array(roll_output_epoch[select],dtype='float64')
        attributes = ['Roll_Epoch', 'ns', 'linear', vardata.min(), vardata.max(),counts_file_low.varattsget(zvars_counts_low[0], expand=True)]
        attrs = attributes[5]
        attrs['VAR_TYPE'] = 'data'
        attributes[5] = attrs
        varinfo = counts_file_low.varinq(zvars_counts_low[0])
        write_var_to_file(sun_spike_noise_file_low, varinfo, vardata, attributes)

        vardata = np.array(noise_total_data_low,dtype='float64')
        attrs = ['noise_vals', [-1.e+31], [vardata.min()], [vardata.max()], 'linear', 'counts', 'series']
        infos = [4, 44, len(vardata), attrs[0], [-9223372036854775807]]
        varattributes = {'CATDESC': attrs[0], 'DEPEND_0': 'Roll_Epoch', 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0],'FILLVAL': np.array(attrs[1], dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0],'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2],'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP': attrs[4]}
        varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1],'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 2, 'Dim_Sizes': [len(noise_total_data_low[0])],'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [],'Pad': np.array(infos[4], dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
        sun_spike_noise_file_low.write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #-------
    # High #
    #-------
    elif select == 1:
        # APPLY THE HIGH COUNTS MASK
        rollchecker = roll_output[select]
        for tme1, ptch1, engy1 in itertools.product(*which_range):

            if tme1%200 == 0:
                print('Filtering bad data for ' + rocket_names[select] + ': ' +  str(round(100 * (tme1/len(epochs[select])),1)) + '%',end='\r')
            elif tme1 == (len(epochs[select]) - 1) and ptch1 == 20 and engy1 == 48:
                print('Filtering bad data for ' + rocket_names[select] + ':' + ' 100%')

            if (tme1>= t_noise_start1_high and tme1<= t_noise_end1_high) or (tme1 >= t_noise_start2_high and tme1<= t_noise_end2_high):
                dat = raw_dat[tme1][ptch1][engy1]
                if dat >= wnoise_threshhold:
                    rchecker = float(rollchecker[tme1])
                    noise_total_data_high.append([dat,rchecker,Energies[engy1],ptch1])

        # -------------------------------------
        # Write out Noise Data
        # -------------------------------------
        vardata = np.array(output_raw, dtype='float64')
        varattributes = counts_file_high.varattsget(zvars_counts_low[18], expand=True)
        varattributes['CATDESC'] = 'counts_masked'
        varattributes['FIELDNAM'] = 'counts_masked'
        varattributes['VALIDMIN'] = [vardata.min(), 'CDF_FLOAT']
        varattributes['VALIDMAX'] = [vardata.max(), 'CDF_FLOAT']
        varinfo = counts_file_high.varinq(zvars_counts_high[18])
        varinfo['Variable'] = 'counts_masked'
        sun_spike_noise_file_high.write_var(varinfo, var_attrs=varattributes, var_data=vardata)


        vardata = np.array(roll_output[select],dtype='float64')
        attrs = ['Roll_Reduced', [-1.e+31], [vardata.min()], [vardata.max()], 'linear', 'deg', 'series']
        infos = [4, 44, len(vardata), attrs[0], [-9223372036854775807]]
        varattributes = {'CATDESC': attrs[0], 'DEPEND_0': 'Roll_Epoch','DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0],'FILLVAL': np.array(attrs[1], dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0],'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2],'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP': attrs[4]}
        varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1],'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [],'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [],'Pad': np.array(infos[4], dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
        sun_spike_noise_file_high.write_var(varinfo, var_attrs=varattributes, var_data=vardata)

        vardata = np.array(roll_output_epoch[select], dtype='float64')
        attributes = ['Roll_Epoch', 'ns', 'linear', vardata.min(), vardata.max(),counts_file_high.varattsget(zvars_counts_high[0], expand=True)]
        attrs = attributes[5]
        attrs['VAR_TYPE'] = 'data'
        attributes[5] = attrs
        varinfo = counts_file_high.varinq(zvars_counts_high[0])
        write_var_to_file(sun_spike_noise_file_high, varinfo, vardata, attributes)

        vardata = np.array(noise_total_data_high,dtype='float64')
        attrs = ['noise_vals', [-1.e+31], [vardata.min()], [vardata.max()], 'linear', 'counts', 'series']
        infos = [4, 44, len(vardata), attrs[0], [-9223372036854775807]]
        varattributes = {'CATDESC': attrs[0], 'DEPEND_0': 'Roll_Epoch', 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0],'FILLVAL': np.array(attrs[1], dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0],'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2],'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP': attrs[4]}
        varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1],'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 2, 'Dim_Sizes': [len(noise_total_data_high[0])],'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [],'Pad': np.array(infos[4], dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
        sun_spike_noise_file_high.write_var(varinfo, var_attrs=varattributes, var_data=vardata)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print("--- %s seconds for Remove_Sun_spikes---" % (time.time() - start_time) ,'\n')