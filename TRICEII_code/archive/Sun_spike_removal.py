# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE REMOVE SUN SPIKES --------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import itertools, time, numpy as np
print('importing variables: ', end='')
start_time = time.time()

from functions import write_var_to_file
from files import counts_file_high,counts_file_low,root
from cdflib import cdfwrite
from Variables import pitch,zvars_counts_high,zvars_counts_low,counts_high_info,counts_low_info,bin_spacing
from Variables import EPOCH_high_counts,EPOCH_low_counts,ranges_counts_high,ranges_counts_low
from Variables import Energies,t_end_low,t_end_high,t_start_high,t_start_low
from Variables import roll_reduced_low,roll_reduced_high,masked_counts_low,masked_counts_high
from Variables import rocket_names,spike_statistics_means_low
from Variables import spike_statistics_means_high,spike_statistics_stddevs_low,spike_statistics_stddevs_high
# from Variables import LOW_ranges,HIGH_ranges,
print('Done')
print('Preparing data for filtering: ',end='')

# ------------
# OUTPUT FILES
# ------------
selects = [0,1]

#angle bounds for roll data characterization
spike_rem_file_counts_high = cdfwrite.CDF(root + 'spike_rem_counts_high',cdf_spec=counts_high_info,delete=True)
spike_rem_file_counts_low = cdfwrite.CDF(root + 'spike_rem_counts_low',cdf_spec=counts_low_info,delete=True)

#--------------------------------
# APPLY NOISE REDUCING FILTER MAP
#--------------------------------
margin = 100 #Controls the range the start/end points for applying the threshold. Set to =0 if you want t_start_time/t_end_time only
std_dev_threshold = 10
std_dev_multiplier = 1
mean_multiplier = 1
adjusted_filtering = True

# Create bins for data to be sorted to. Xaxis is roll angle, yaxis is energy:

no_of_points = int(360 / (bin_spacing ))
roll_bins = np.linspace(-180,180,no_of_points+1)

roll_reduced = [roll_reduced_low,roll_reduced_high]
masked_counts = [masked_counts_low,masked_counts_high]
ranges = [ranges_counts_low,ranges_counts_high]
percents = [[],[]]
starts_ends = [[t_start_low,t_end_low],[t_start_high,t_end_high] ]
# Rollranges = [LOW_ranges,HIGH_ranges]
count_epochs = [EPOCH_low_counts,EPOCH_high_counts]
noise_means_data = [spike_statistics_means_low,spike_statistics_means_high]
noise_std_data = [spike_statistics_stddevs_low,spike_statistics_stddevs_high]


print('Done')

for select in selects:
    which_range = ranges[select]
    wstarts_ends =starts_ends[select]
    wmask_counts = masked_counts[select]
    output = np.zeros(shape=(len(wmask_counts), len(wmask_counts[0]), len(wmask_counts[0][0])))
    rollchecker = roll_reduced[select]
    # wrollranges = Rollranges[select]
    wcount_epoch = count_epochs[select]
    wnoise_means_list = noise_means_data[select]
    wnoise_std_list = noise_std_data[select]

    # APPLY THE LOW COUNTS THRESHOLD
    for tme1, ptch1, engy1 in itertools.product(*which_range):
        if tme1%200 == 0:
            print('Applying noise map for ' + rocket_names[select] + ': ' +  str(round(100 * (tme1/len(wcount_epoch)),1)) + '%',end='\r')
        elif tme1 == (len(wcount_epoch)-1) and ptch1 ==20 and engy1 == 48:
            print('Applying noise map for ' + rocket_names[select] + ':' + ' 100%')

        dat = wmask_counts[tme1][ptch1][engy1]
        rchecker = float(rollchecker[tme1])

        # ------------------------------
        # Code to apply noise filter map
        # ------------------------------

        if tme1>= (wstarts_ends[0] - margin) and tme1<= (wstarts_ends[1]+margin):
            # if (rchecker >= wrollranges[1][0] and rchecker <= wrollranges[1][1]) or (rchecker >= wrollranges[0][0] and rchecker <= wrollranges[0][1]):

            roll_index = np.abs(roll_bins - rchecker).argmin()
            mapped_noise_mean = wnoise_means_list[ptch1][roll_index][engy1]
            mapped_noise_std = wnoise_std_list[ptch1][roll_index][engy1]

            # CONDITIONAL STATEMENT ABOUT THE STANDARD DEVIATION

            #Ajust some correction factors for the pads
            if adjusted_filtering:
                if select == 0:
                    if ptch1 == 9:
                        mean_multiplier = 1.5
                        std_dev_multiplier = 2
                    elif ptch1 == 11 and engy1 < 25 :
                        mean_multiplier = 4.3
                        std_dev_multiplier = 2
                    elif ptch1 == 11 and engy1 >= 25 :
                        mean_multiplier = 3
                        std_dev_multiplier = 3
                    elif ptch1 == 10 and engy1 < 25 :
                        mean_multiplier = 4
                        std_dev_multiplier = 2
                    elif ptch1 == 10 and engy1 >= 25 :
                        mean_multiplier = 2.5
                        std_dev_multiplier = 2.5
                    elif ptch1==12 and engy1 <25:
                        mean_multiplier = 3.3
                        std_dev_multiplier = 3
                    elif ptch1==12 and engy1 >=25:
                        mean_multiplier = 3
                        std_dev_multiplier = 3
                    else:
                        mean_multiplier = 1
                        std_dev_multiplier = 1

                elif select == 1:
                    if (ptch1 == 10 or ptch1 in range(14,21)) and engy1 <25:
                        mean_multiplier = 2
                        std_dev_multiplier = 2
                    elif (ptch1 == 10 or ptch1 in range(14,21)) and engy1 >=25:
                        mean_multiplier = 1.6
                        std_dev_multiplier = 1.5
                    elif ptch1 in range(0,10):
                        mean_multiplier = 1.7
                        std_dev_multiplier = 1.7
                    elif (ptch1 == 11 ) and engy1 <25:
                        mean_multiplier = 2
                        std_dev_multiplier = 3
                    elif (ptch1==12 or ptch1 ==13) and engy1 < 25:
                        mean_multiplier = 2.8
                        std_dev_multiplier = 3
                    elif (ptch1==11 or ptch1==12 or ptch1 ==13) and engy1 >= 25:
                        mean_multiplier = 1.9
                        std_dev_multiplier = 1.5

                    else:
                        mean_multiplier = 1
                        std_dev_multiplier = 1



            if mapped_noise_mean !=0:
                percents[select].append((100 * std_dev_threshold / mapped_noise_mean))
                if mapped_noise_std >= (100 *std_dev_threshold / mapped_noise_mean):
                    if (dat - (mean_multiplier*mapped_noise_mean + std_dev_multiplier*mapped_noise_std)) < 0:
                        dat = 0
                    else:
                        dat = (dat - (mean_multiplier*mapped_noise_mean + std_dev_multiplier*mapped_noise_std))
                else:
                    if (dat - (mean_multiplier*mapped_noise_mean)) < 0:
                        dat = 0
                    else:
                        dat = (dat - (mean_multiplier*mapped_noise_mean))
            elif mapped_noise_mean ==0:
                percents[select].append('N/A')

        output[tme1][ptch1][engy1] = dat

    # print(percents[select])

    # low
    if select == 0:

        # -------------------------------------
        # Write out Counts with cutoffs applied
        # -------------------------------------
        vardata = output
        attributes = ['removed_spike_counts_threshed_data', 'Counts', 'log', output.min(), output.max(),counts_file_low.varattsget(zvars_counts_low[18], expand=True)]
        varinfo = counts_file_low.varinq(zvars_counts_low[18])
        write_var_to_file(spike_rem_file_counts_low, varinfo, vardata, attributes)

        # -------------------------------------
        # Write out EPOCH,pitch,energy
        # -------------------------------------
        vardata = EPOCH_low_counts
        attributes = ['Epoch', 'ns', 'linear', EPOCH_low_counts.min(), EPOCH_low_counts.max(),counts_file_low.varattsget(zvars_counts_low[0], expand=True)]
        varinfo = counts_file_low.varinq(zvars_counts_low[0])
        write_var_to_file(spike_rem_file_counts_low, varinfo, vardata, attributes)

        vardata = pitch
        attributes = ['Pitch_Angle', 'deg', 'linear', pitch.min(), pitch.max(),counts_file_low.varattsget(zvars_counts_low[9], expand=True)]
        varinfo = counts_file_low.varinq(zvars_counts_low[9])
        write_var_to_file(spike_rem_file_counts_low, varinfo, vardata, attributes)

        vardata = Energies
        attributes = ['Energy', 'eV', 'log', Energies.min(), Energies.max(),counts_file_low.varattsget(zvars_counts_low[16], expand=True)]
        varinfo = counts_file_low.varinq(zvars_counts_low[16])
        write_var_to_file(spike_rem_file_counts_low, varinfo, vardata, attributes)

    # High
    elif select == 1:
        # -------------------------------------
        # Write out Counts with cutoffs applied
        # -------------------------------------
        vardata = output
        attributes = ['removed_spike_counts_threshed_data', 'Counts', 'log', output.min(), output.max(),counts_file_high.varattsget(zvars_counts_high[18], expand=True)]
        varinfo = counts_file_high.varinq(zvars_counts_high[18])
        write_var_to_file(spike_rem_file_counts_high, varinfo, vardata, attributes)

        # -------------------------------------
        # Write out EPOCH,pitch,energy
        # -------------------------------------
        vardata = EPOCH_high_counts
        attributes = ['Epoch', 'ns', 'linear', EPOCH_high_counts.min(), EPOCH_high_counts.max(),counts_file_high.varattsget(zvars_counts_high[0], expand=True)]
        varinfo = counts_file_high.varinq(zvars_counts_high[0])
        write_var_to_file(spike_rem_file_counts_high, varinfo, vardata, attributes)

        vardata = pitch
        attributes = ['Pitch_Angle', 'deg', 'linear', pitch.min(), pitch.max(),counts_file_high.varattsget(zvars_counts_high[9], expand=True)]
        varinfo = counts_file_high.varinq(zvars_counts_high[9])
        write_var_to_file(spike_rem_file_counts_high, varinfo, vardata, attributes)

        vardata = Energies
        attributes = ['Energy', 'eV', 'log', Energies.min(), Energies.max(),counts_file_high.varattsget(zvars_counts_high[16], expand=True)]
        varinfo = counts_file_high.varinq(zvars_counts_high[16])
        write_var_to_file(spike_rem_file_counts_high, varinfo, vardata, attributes)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print("--- %s seconds for Sun_spike_removal---" % (time.time() - start_time) ,'\n')