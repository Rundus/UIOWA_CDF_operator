# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE CDF INFOMATION FILE --------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import time
import itertools
import cdflib,math,numpy as np
from files import counts_file_low,counts_file_high,root
from Variables import EPOCH_low_counts,t_start_low,t_end_low,pitch,t_start_high,t_end_high,Energies,ranges_counts_low,counts_low_info,zvars_counts_low,zvars_counts_high,ranges_counts_high,EPOCH_high_counts
from cdflib import cdfread,cdfwrite
from scipy.fft import rfft,irfft,rfftfreq
from scipy.signal import detrend
import matplotlib.pyplot as plt
from functions import write_var_to_file

print('Importing variables: ',end='')



# file1 = root + 'spike_rem_counts_low.cdf'
# file2 = root +'spike_rem_counts_high.cdf'
# data_file1 = cdfread.CDF(file1)
# data_file2 = cdfread.CDF(file2)
#
# output_low_filter = data_file1.varget('Spike_counts_thresh_data')
# output_high_filter = data_file2.varget('Spike_counts_thresh_data')
# output_low_diff = data_file1.varget('Spike_counts_thresh_data')
# output_high_diff = data_file2.varget('Spike_counts_thresh_data')


file_low = root + 'TRICE_52004_l1_eepaa_20181208T082243_v1.1.2_COUNTS.cdf'
file_high = root +'TRICE_52003_l1_eepaa_20181208T082239_v1.1.2_COUNTS.cdf'
data_file1 = cdfread.CDF(file_low)
data_file2 = cdfread.CDF(file_high)

output_low_filter = data_file1.varget('eepaa')
output_high_filter = data_file2.varget('eepaa')
output_low_diff = data_file1.varget('eepaa')
output_high_diff = data_file2.varget('eepaa')






spike_rem_file_counts_low = cdfwrite.CDF(root + 'spike_rem_counts_low_CDFINFO',cdf_spec=counts_low_info,delete=True)
spike_rem_file_counts_high = cdfwrite.CDF(root + 'spike_rem_counts_high_CDFINFO',cdf_spec=counts_low_info,delete=True)

# -----------------
# TOGGLES
# -----------------
DelMethod = False #Implement the deletion method or not
select = 0 #Which file to write out to
files = [spike_rem_file_counts_low,spike_rem_file_counts_high]
wrange = [range(t_start_low,t_end_low),range(t_start_high,t_end_high)]
wtrange = [ranges_counts_low,ranges_counts_high]
wattrs = [counts_file_low.varattsget(zvars_counts_low[18],expand=True),counts_file_high.varattsget(zvars_counts_high[18],expand=True)]
wofilter = [output_low_filter,output_high_filter]
output_filter = wofilter[select]
wepoch = [EPOCH_low_counts,EPOCH_high_counts]
EPOCH = wepoch[select]
wvarinfo = [counts_file_low.varinq(zvars_counts_low[18]),counts_file_high.varinq(zvars_counts_high[18])]
wdiff = [output_low_diff,output_high_diff]
output_diff = wdiff[select]
titles = ['LOW','HIGH']

#GET NOISEY FREQ DATA
pslicer = [i for i in range(len(pitch))]
eslicer = [j for j in range(10)]
# sr = 100 #Sample RAtes for hte FFTs
sr = 20 # 1 / 50ms per sweep means => 20 sweeps a second
#-------------

if select == 0:
    vardata = EPOCH
    attributes = ['Epoch', 'ns', 'linear', EPOCH.min(), EPOCH.max(),counts_file_low.varattsget(zvars_counts_low[0], expand=True)]
    varinfo = counts_file_low.varinq(zvars_counts_low[0])
    write_var_to_file(files[select], varinfo, vardata, attributes)

    vardata = pitch
    attributes = ['Pitch_Angle', 'deg', 'linear', pitch.min(), pitch.max(),counts_file_low.varattsget(zvars_counts_low[9], expand=True)]
    varinfo = counts_file_low.varinq(zvars_counts_low[9])
    write_var_to_file(files[select], varinfo, vardata, attributes)

    vardata = Energies
    attributes = ['Energy', 'eV', 'log', Energies.min(), Energies.max(),counts_file_low.varattsget(zvars_counts_low[16], expand=True)]
    varinfo = counts_file_low.varinq(zvars_counts_low[16])
    write_var_to_file(files[select], varinfo, vardata, attributes)

elif select ==1:
    vardata = EPOCH
    attributes = ['Epoch', 'ns', 'linear', EPOCH.min(), EPOCH.max(),counts_file_high.varattsget(zvars_counts_high[0], expand=True)]
    varinfo = counts_file_high.varinq(zvars_counts_high[0])
    write_var_to_file(files[select], varinfo, vardata, attributes)

    vardata = pitch
    attributes = ['Pitch_Angle', 'deg', 'linear', pitch.min(), pitch.max(),counts_file_high.varattsget(zvars_counts_high[9], expand=True)]
    varinfo = counts_file_high.varinq(zvars_counts_high[9])
    write_var_to_file(files[select], varinfo, vardata, attributes)

    vardata = Energies
    attributes = ['Energy', 'eV', 'log', Energies.min(), Energies.max(),counts_file_high.varattsget(zvars_counts_high[16], expand=True)]
    varinfo = counts_file_high.varinq(zvars_counts_high[16])
    write_var_to_file(files[select], varinfo, vardata, attributes)


output_fft_data = np.zeros(shape=(len(output_filter),21,49))

for tmee,pitche,engye in itertools.product(*wtrange[select]):
    output_fft_data[tmee][pitche][engye] =output_filter[tmee][pitche][engye]

print('Done')
print('Performing FFT and sorting data to get the NOISE: ',end='' )


# -----------------------------------------
# USE ENERGIES 0-10 TO GET THE NOISE SIGNAL
# -----------------------------------------



power = []
freq = []
power1 = []
freq1 = []


for pslice in pslicer:
    for eslice in eslicer:
        hdat = detrend(np.array([output_filter[tme2][pslice][eslice] for tme2 in wrange[select]]),type='linear')
        duration = sum([EPOCH_low_counts[i+1] - EPOCH_low_counts[i] for i in wrange[select]])/(10**9) # seconds
        N = int(sr * duration)
        hdat_fft = rfft(hdat)
        hdat_freq = rfftfreq(N, 1/sr)

        #Gather all data
        for i in range(len(hdat_fft)):
            y_val = np.abs(hdat_fft[i])
            x_val = hdat_freq[i]
            power.append(y_val)
            freq.append(x_val)


# To characterize the sun-spike noise we noticed it is the spiky lines on the FFT higher frequencies. We see to find frequencies with larger max values and delete that noise

noisey_freq = []
noisey_indicies = []

if DelMethod:
    # Sift through frequencies and find the max value of each freq and compare to max next to it
    for i in range(len(freq)):
        if freq[i] > 0.4 and freq[i] < 0.8 :
            if power[i] >= 10000:
                noisey_freq.append(freq[i])

    noisey_freqs = np.unique(noisey_freq)

    for freqval in noisey_freqs:
        noisey_indicies.append(np.where(freq == freqval))

    for array in noisey_indicies:
        for list in array:
            for val in list:
                power[val] = 0



# -----------------------
# GETTING THE NOISE - END
# -----------------------
print('Done')
print('Performing FFT and writing DATA: ',end='' )

# ----------------------------------------
# USE THE NOISE SIGNATURE TO CLEAN UP DATA
# ----------------------------------------
pslicer1 = [i for i in range(len(pitch))]
eslicer1 = [j for j in range(len(Energies))]

for pslice in pslicer1:
    for eslice in eslicer1:
        hdat = detrend(np.array([output_filter[tme2][pslice][eslice] for tme2 in wrange[select]]),type='linear')
        duration = sum([EPOCH[i+1] - EPOCH[i] for i in wrange[select]])/(10**9) # seconds
        N = int(sr * duration)
        hdat_fft = rfft(hdat)
        hdat_freq = rfftfreq(N, 1/sr)

        #Delete data in boxed range
        for i in range(len(hdat_fft)):
            y_val = np.abs(hdat_fft[i])
            x_val = hdat_freq[i]

            if DelMethod:
                if x_val in noisey_freqs :
                    hdat_fft[i] = 0

            power1.append(np.abs(hdat_fft[i]))
            freq1.append(x_val)


        #Reverse Transform
        hdat_ifft = irfft(hdat_fft)


        for tme4 in range(len(hdat_ifft)):
            dat = hdat_ifft[tme4]

            if dat < 0:
                output_fft_data[tme4 + t_start_low][pslice][eslice] = 0
            else:
                output_fft_data[tme4 + t_start_low][pslice][eslice] = round(dat)


vardata = np.array(output_fft_data)
attributes = ['Spike_counts_fft_data', 'Counts', 'log', vardata.min(), vardata.max() ,wattrs[select]]
varinfo = wvarinfo[select]
write_var_to_file(files[select], varinfo, vardata, attributes)


for tmee,pitche,engye in itertools.product(*wtrange[select]):
    dat = output_filter[tmee][pitche][engye] - output_fft_data[tmee][pitche][engye]
    if dat <=0:
        output_diff[tmee][pitche][engye] = 0
    else:
        output_diff[tmee][pitche][engye] = dat

vardata = np.array(output_diff)
attributes = ['Spike_counts_filter_minus_fft', 'Counts', 'log',vardata.min(), vardata.max() ,wattrs[select]]
varinfo = wvarinfo[select]
write_var_to_file(files[select], varinfo, vardata, attributes)


print('Done')



plt.subplot(2,1,1)
plt.plot(np.array(freq),power,'.')
# plt.title(titles[select] + '\nPitch: ' + str(pitch[pslice]) + 'deg')
plt.title(titles[select] + '\n All Pitches\n 11486eV to 4427eV ')
# plt.ylim(0,20000)
plt.ylabel('Power')
plt.xlabel('Freq [Hz]')

plt.subplot(2,1,2)
plt.plot(np.array(freq1),power1,',')
plt.title(titles[select] + '\n All Pitches, All Enegies')
# plt.ylim(0,20000)
plt.ylabel('Power')
plt.xlabel('Freq [Hz]')
plt.tight_layout(-0.3)
plt.show()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%















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



