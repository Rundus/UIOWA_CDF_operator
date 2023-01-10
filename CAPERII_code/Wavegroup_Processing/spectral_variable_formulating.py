# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -------- Reading DC Spectral data and correlation-------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import cdflib
from matplotlib import pyplot as plt
import numpy as np
from files_p import file_freqs,file_LH_mag,file_RH_mag,file_LR_ratio,file_times,file_integrate
from Variables_p import epoch_1,epoch_2,t_end,t_start,t_rocket_start


# -----------------
# Get the Time Data
# -----------------

file = open(file_times,'r')
file_str_dat =file.readlines()
file_dat = []
for i in range(len(file_str_dat)):
    temp_list = file_str_dat[i].strip('\n').split(',')
    file_dat.append(float(temp_list[0]))
dat_times = np.array(file_dat)

# -----------------------------------------
# Get the Integrated DC data and manipulate
# -----------------------------------------


# Read in the Integrated Data and assign it to column array dat_integrate array
file = open(file_integrate,'r')
file_str_dat =file.readlines()
file_dat = []

for i in range(len(file_str_dat)):
    if i != 0:
        temp_list = file_str_dat[i].strip('\n').split(',')
        file_dat.append([float(string) for string in temp_list])
    else:
        headers = file_str_dat[i].strip('\n').strip('#').split(',')
        file_dat.append(headers)

dat = file_dat
dat_integrate = [[],[],[],[]]
dat_integrate_headers = []
for i in range(len(dat)):
    if i == 0:
        for j in range(4):
            dat_integrate_headers.append(dat[i][j])
    else:
        dat_integrate[0].append(dat[i][0])
        dat_integrate[1].append(dat[i][1])
        dat_integrate[2].append(dat[i][2])
        dat_integrate[3].append(dat[i][3])




# -----------------------------
# normalize the integrated data
# -----------------------------
dat_integrate_norm = [[],[],[],[]]

for row in range(4):
    max_val = max(dat_integrate[row])
    for i in range(len(dat_integrate[row])):
        dat_integrate_norm[row].append( dat_integrate[row][i]/max_val)






# Close all Files
file.close()
