# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE SUN SPIKES STATISTICS- -------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np
import time
import random

print('importing variables: ', end='')
start_time = time.time()

from files import root
from cdflib import cdfwrite
from functions import write_var_to_file
from Variables import counts_high_info,counts_low_info,Energies,noise_vals_low,noise_vals_high,rocket_names
from Variables import counts_file_high, zvars_counts_high,bin_spacing,pitch
print('Done')


# ------------
# OUTPUT FILES
# ------------
spike_rem_file_statistics_high = cdfwrite.CDF(root + 'spike_rem_statistics_high',cdf_spec=counts_high_info,delete=True)
spike_rem_file_statistics_low = cdfwrite.CDF(root + 'spike_rem_statistics_low',cdf_spec=counts_low_info,delete=True)

#Some variables
selects = [0,1]

#--------------------------------------------
# Noise data used to eliminate noise in data:
#--------------------------------------------
print('Constructing variables for Noise data:', end='')
noise = [noise_vals_low,noise_vals_high]
counts = [[],[]]
Roll = [[],[]]
EnergiesBin = [[],[]]
EnergiesIndexes = [[],[]]
PitchA = [[],[]]

for j in range(2):
    ndat = noise[j]
    for i in range(len(ndat)):
        counts[j].append(ndat[i][0])
        Roll[j].append(ndat[i][1])
        EnergiesBin[j].append(ndat[i][2])
        EnergiesIndexes[j].append(np.where(Energies == ndat[i][2])[0][0])
        PitchA[j].append(ndat[i][3])


countsa = np.array(counts,dtype='object')
Rolla = np.array(Roll,dtype='object')
Energiesa = np.array(EnergiesBin,dtype='object')
EnergiesIndexesa = np.array(EnergiesIndexes,dtype='object')
Pitcha = np.array(PitchA,dtype='object')
print('Done')



# -------------------------------------
# CREATE NOISE MAP FOR DATA TO SUBTRACT
# -------------------------------------

# Create bins for data to be sorted to. Xaxis is roll angle, yaxis is energy:
no_of_points = int(360 / (bin_spacing ))
roll_bins = np.linspace(-180,180,no_of_points+1)


for select in selects:
    Engyindicies = EnergiesIndexesa[select]
    roll_vals = Rolla[select]
    count_vals = countsa[select]
    pitch_vals = Pitcha[select]
    noise_means_list = [[[] for i in range(len(Energies))] for j in range(len(roll_bins))]
    noise_stddev_list = [[[] for i in range(len(Energies))] for j in range(len(roll_bins))]
    noise_list = [[[[] for i in range(len(Energies))] for j in range(len(roll_bins))] for k in range(len(pitch))]
    noise_means_total = [[[[] for i in range(len(Energies))] for j in range(len(roll_bins))] for k in range(len(pitch))]
    noise_stds_total = [[[[] for i in range(len(Energies))] for j in range(len(roll_bins))] for k in range(len(pitch))]

    #Place all the count values in a sorted bin
    for i in range(len(count_vals)):
        if i % 200 == 0:
            print('Sorting count values for ' + rocket_names[select] + ': ' + str(round(100 * (i / len(count_vals)), 1)) + '%', end='\r')
        elif i == (len(count_vals) - 1):
            print('Sorting count values for ' + rocket_names[select] + ':' + ' 100%')

        #Determine which roll bin (x val) the data belongs to (PART THAT TAKE SO LONG)
        x_index_pos = np.abs(roll_bins - roll_vals[i]).argmin()
        y_index_pos = Engyindicies[i]
        z_index_pos = int(pitch_vals[i])
        noise_list[z_index_pos][x_index_pos][y_index_pos].append(count_vals[i])

    #Convert all the lists into averages

    for ptch in range(len(pitch)):
        for i in range(len(roll_bins)):

            if i != (len(roll_bins) - 1):
                print('Creating noise map for ' + rocket_names[select] + ' at '+ str(pitch[ptch]) +' deg: ' + str(  round(   (100 * i / len(roll_bins)), 1)) + '%', end='\r' )
            elif i == (len(roll_bins) - 1):
                print('Creating noise map for ' + rocket_names[select] + ' at '+ str(pitch[ptch]) +' deg: '+ ' 100%')

            for j in range(len(Energies)):
                val = noise_list[ptch][i][j]

                if len(val) == 0:
                    noise_means_total[ptch][i][j] = 0
                    noise_stds_total[ptch][i][j] = 0
                elif len(val) == 1:
                    noise_means_total[ptch][i][j] = val[0]
                    noise_stds_total[ptch][i][j] = 0
                elif len(val) > 1:
                    noise_means_total[ptch][i][j] = round(np.mean(val),0)
                    noise_stds_total[ptch][i][j] = np.std(val)

    print(len(pitch))
    # WRITE DATA OUT TO FILE:

    if select ==0:
        vardata = np.array(noise_means_total, dtype='float64')
        attrs = ['noise_means', [-1.e+31], [vardata.min()], [vardata.max()], 'linear', 'counts', 'series']
        infos = [4, 44, len(vardata), attrs[0], [-9223372036854775807]]
        varattributes = {'CATDESC': attrs[0], 'DEPEND_0':'Pitch_Angle','DEPEND_1':'roll_bins', 'DEPEND_2': 'Energy', 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0],
                         'FILLVAL': np.array(attrs[1], dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0],
                         'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2],
                         'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP': attrs[4]}
        varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1],
                   'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 2, 'Dim_Sizes': [len(roll_bins),49],
                   'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [-1,-1],
                   'Pad': np.array(infos[4], dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
        spike_rem_file_statistics_low.write_var(varinfo, var_attrs=varattributes, var_data=vardata)



        vardata = np.array(noise_stds_total, dtype='float64')
        attrs = ['noise_stds', [-1.e+31], [vardata.min()], [vardata.max()], 'linear', 'counts', 'series']
        infos = [4, 44, len(vardata), attrs[0], [-9223372036854775807]]
        varattributes = {'CATDESC': attrs[0],'DEPEND_0':'Pitch_Angle' , 'DEPEND_1':'roll_bins', 'DEPEND_2': 'Energy', 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0],
                         'FILLVAL': np.array(attrs[1], dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0],
                         'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2],
                         'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP': attrs[4]}
        varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1],
                   'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 2, 'Dim_Sizes': [len(roll_bins),49],
                   'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [-1,-1],
                   'Pad': np.array(infos[4], dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
        spike_rem_file_statistics_low.write_var(varinfo, var_attrs=varattributes, var_data=vardata)

        vardata = Energies
        attributes = ['Energy', 'eV', 'log', Energies.min(), Energies.max(),counts_file_high.varattsget(zvars_counts_high[16], expand=True)]
        varinfo = counts_file_high.varinq(zvars_counts_high[16])
        write_var_to_file(spike_rem_file_statistics_low, varinfo, vardata, attributes)

        vardata = np.array(pitch)
        attributes = ['Pitch_Angle', 'deg', 'linear', vardata.min(), vardata.max(),counts_file_high.varattsget(zvars_counts_high[9], expand=True)]
        varinfo = counts_file_high.varinq(zvars_counts_high[9])
        write_var_to_file(spike_rem_file_statistics_low, varinfo, vardata, attributes)

        vardata = roll_bins
        attrs = ['roll_bins', [-1.e+31], [vardata.min()], [vardata.max()], 'linear', 'deg', 'series']
        infos = [4, 44, len(vardata), attrs[0], [-9223372036854775807]]
        varattributes = {'CATDESC': attrs[0], 'DEPEND_0': 'roll_bins', 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0],
                         'FILLVAL': np.array(attrs[1], dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0],
                         'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2],
                         'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP': attrs[4]}
        varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1],
                   'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 1, 'Dim_Sizes': [],
                   'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [],
                   'Pad': np.array(infos[4], dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
        spike_rem_file_statistics_low.write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    elif select ==1:
        vardata = np.array(noise_means_total, dtype='float64')
        attrs = ['noise_means', [-1.e+31], [vardata.min()], [vardata.max()], 'linear', 'counts', 'series']
        infos = [4, 21, len(vardata), attrs[0], [-9223372036854775807]]
        varattributes = {'CATDESC': attrs[0], 'DEPEND_0':'Pitch_Angle','DEPEND_1': 'roll_bins', 'DEPEND_2': 'Energy', 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0],
                         'FILLVAL': np.array(attrs[1], dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0],
                         'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2],
                         'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP': attrs[4]}
        varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1],
                   'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 2, 'Dim_Sizes': [len(roll_bins),49],
                   'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [-1,-1],
                   'Pad': np.array(infos[4], dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
        spike_rem_file_statistics_high.write_var(varinfo, var_attrs=varattributes, var_data=vardata)

        vardata = np.array(noise_stds_total, dtype='float64')
        attrs = ['noise_stds', [-1.e+31], [vardata.min()], [vardata.max()], 'linear', 'counts', 'series']
        infos = [4, 21, len(vardata), attrs[0], [-9223372036854775807]]
        varattributes = {'CATDESC': attrs[0],'DEPEND_0':'Pitch_Angle', 'DEPEND_1': 'roll_bins', 'DEPEND_2': 'Energy', 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0],
                         'FILLVAL': np.array(attrs[1], dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0],
                         'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2],
                         'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP': attrs[4]}
        varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1],
                   'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 2, 'Dim_Sizes': [len(roll_bins),49],
                   'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [-1,-1],
                   'Pad': np.array(infos[4], dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
        spike_rem_file_statistics_high.write_var(varinfo, var_attrs=varattributes, var_data=vardata)

        vardata = Energies
        attributes = ['Energy', 'eV', 'log', Energies.min(), Energies.max(),counts_file_high.varattsget(zvars_counts_high[16], expand=True)]
        varinfo = counts_file_high.varinq(zvars_counts_high[16])
        write_var_to_file(spike_rem_file_statistics_high, varinfo, vardata, attributes)

        vardata = np.array(pitch)
        attributes = ['Pitch_Angle', 'deg', 'linear', vardata.min(), vardata.max(),counts_file_high.varattsget(zvars_counts_high[9], expand=True)]
        varinfo = counts_file_high.varinq(zvars_counts_high[9])
        write_var_to_file(spike_rem_file_statistics_high, varinfo, vardata, attributes)

        vardata = roll_bins
        attrs = ['roll_bins', [-1.e+31], [vardata.min()], [vardata.max()], 'linear', 'deg', 'series']
        infos = [4, 44, len(vardata), attrs[0], [-9223372036854775807]]
        varattributes = {'CATDESC': attrs[0], 'DEPEND_0': 'roll_bins', 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0],
                         'FILLVAL': np.array(attrs[1], dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0],
                         'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2],
                         'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP': attrs[4]}
        varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1],
                   'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 1, 'Dim_Sizes': [],
                   'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [],
                   'Pad': np.array(infos[4], dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
        spike_rem_file_statistics_high.write_var(varinfo, var_attrs=varattributes, var_data=vardata)






print("--- %s seconds for Sun_spike_noise_statistics ---" % (time.time() - start_time) ,'\n')