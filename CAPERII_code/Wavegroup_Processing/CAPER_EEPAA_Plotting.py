# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -------- CAPERII graphing one pitch, one energy -------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



import numpy as np
from matplotlib import pyplot as plt
from cdflib import cdfread
from epoch_adjusting import epoch_adjusted_1_mins,epoch_adjusted_2_mins,epoch_adjusted_1_sec,epoch_adjusted_2_sec,epoch_adjusted_1,epoch_adjusted_2
from files_p import files
from Variables_p import adjust,adjust_min,adjust_sec

# -----------------------
# ------- TOGGLES -------
# -----------------------

ans_file = input('\n\neepaa 1 or eepaa 2?\n (1 or 2)\n Ans: ')
select = int(ans_file) - 1
ans = input('\n\nSlice in?\n 0 - time\n 1 - pitch\n 2- energy\n Ans: ')

# -------------------- TOGGLES --------------------0
File =  cdfread.CDF(files[select])
slice = ['time','pitch','energy']
names = ['eepaa1','eepaa2']
plot_min = -100
plot_max = 1000
FilterQ = False
Filtert = 1000
Filterb = 250
# ------------------------------------------------



# Get the appropriate variables
counts = File.varget('eepaa')
epoch = File.varget('Epoch')
epoch_adjusted = adjust[select].tolist()
pitches = File.varget('Pitch_Angle').tolist()
energies = File.varget('Energy').tolist()
# energies.reverse()







# ----------------------------------------------------------------
# ------- PLOT THE CAPER COUNTS DATA - FILTERED/UNFILTERED -------
# ----------------------------------------------------------------

if ans == '0':
    slice_time = input('\nWhat slice in time?\nAns:')
    # Plotting - slices in time
    dat = counts[int(slice_time)]
    temp = np.transpose(dat)
    Z = temp.tolist()


    x = pitches
    y = energies
    plt.pcolormesh(x, y, Z, shading='auto', vmin=plot_min, vmax=plot_max, cmap='rainbow')
    plt.title(names[select] +' Time Slice t=' + slice_time)
    plt.xlabel('Pitch Angle (deg)')
    plt.ylabel('Energy (eV)')
    plt.colorbar(label='Counts')



elif ans == '1':
    slice_pitch = input('\nWhat slice in pitch? (0-20)\nAns:')
    data_slice = []

    for i in range(len(adjust[select])):
        dat = counts[i][int(slice_pitch)]
        dat.tolist()
        data_slice.append(dat)
    Z = np.transpose(np.array(data_slice)).tolist()

    # Apply a Filter for 0 and 1
    if FilterQ:
        plot_min = 0
        plot_max = 1
        for i in range(len(Z)):
            for j in range(len(Z[0])):
                if Z[i][j] >= Filtert or Z[i][j] <= Filterb:
                    Z[i][j] = 0
                else:
                    Z[i][j] = 1
        # Plotting - slices in time
        x = adjust[select]
        y = energies
        plt.pcolormesh(x, y, Z, shading='gouraud', vmin=plot_min, vmax=plot_max, cmap='rainbow')
        plt.title(names[select] + ' Pitch Slice =' + str(
            pitches[int(slice_pitch)]) + ' deg\n' + 'Thresholds - min:' + str(Filterb) + ' max:' + str(Filtert))
        plt.ylabel('Energy (eV)')
        plt.xlabel('Epoch (min)')
        plt.colorbar(label='Counts')
        plt.yscale('log')

    else:
        # Plotting - slices in time
        x = adjust[select]
        y = energies
        plt.pcolormesh(x, y, Z, shading='gouraud', vmin=plot_min, vmax=plot_max, cmap='rainbow')
        plt.title(names[select] + ' Pitch Slice =' + str(pitches[int(slice_pitch)]) + ' deg')
        plt.ylabel('Energy (eV)')
        plt.xlabel('Epoch (min)')
        plt.colorbar(label='Counts')
        plt.yscale('log')




elif ans =='2':
    slice_energy = input('\nWhat slice in energy?(0 to 48)\nAns:')
    data_slice = []

    for i in range(len(adjust[select])):
        temp_list = []
        for j in range(21):
            temp_list.append(counts[i][j][int(slice_energy)])
        data_slice.append(temp_list)

    Z = np.transpose(np.array(data_slice)).tolist()
    # Plotting - slices in time
    x = adjust[select]
    y = pitches
    plt.pcolormesh(x,y, Z, shading='gouraud', vmin=plot_min, vmax=plot_max, cmap='rainbow')
    plt.title(names[select] + ' Pitch Slice =' + str(energies[int(slice_energy)]) + ' eV')
    plt.ylabel('Pitches (deg)')
    plt.xlabel('Epoch (min)')
    plt.colorbar(label='Counts')


plt.show()



