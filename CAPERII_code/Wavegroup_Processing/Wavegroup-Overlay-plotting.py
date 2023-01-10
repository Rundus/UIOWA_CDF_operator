# ------------------------------------------------------------------------
# File designed to overlay the CAPERII EEPAA data on top of Spencer's Data
# ------------------------------------------------------------------------

import cdflib
import matplotlib.dates as mdates
import numpy as np,matplotlib.pyplot as plt
import matplotlib.colors as colors
from Variables_p import counts_1,counts_2,epoch_1,epoch_2,energies_1,energies_2,pitches_1,t_end,t_start
from Variables_p import  correlators,energy_correlator,phase_bins,t_start_correlator,t_rocket_start_correlator,t_end_correlator,epoch_correlator,t_rocket_start
from spectral_variable_formulating import dat_integrate,dat_integrate_headers,dat_times


# %%%%%%%%%%%%%%%%%%%%%%%
# PLOTTING THE EEPAA DATA
# %%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%
# TOGGLES
#%%%%%%%%

#EEPAA GRAPH

changerange = True

if changerange:
    t_start = 13000 # control when EEPAA spectrogram starts
    t_end = 15000 # Control when EEPAA spectrogram ends

#Spectral Wave
which_freq = [0] # (pick values between 0 and 3)
num_of_ticks = 50 # Number of x-ticks on plot
myFmt = mdates.DateFormatter('%H:%M:%S:%MS')  # Format the labels on x-tiks
ylim_low,ylim_high = 0,200 # set limits on Spectral data plot
line_colors = ['black','red','brown','purple'] # Sets the colors that the overlayed spectral data can be

#General
tight_layout = False
xinches,yinches = 20,10 # X,Y dimensions of plot

#Colorbar - Log/Min/Max
colorbarlog = False
logy = True

#%%%%%%%%%%%%
# END TOGGLES
#%%%%%%%%%%%%



# -------------------------------------------------------
# ALIGINGING THE SPECTRAL INTEGRATED DATA WITH EEPAA DATA
# -------------------------------------------------------
# Assuming t=0 for spectral occurs at 9:27:00.034.590.000 EEPAA, which is t_rocket_start = 3772

dat_times_EEPAA = cdflib.epochs.CDFepoch.to_datetime(    np.array( [    (dat_times[i]*(10**(9))) +epoch_1[t_rocket_start] for i in range(len(dat_times))], dtype='int64'))



# Get the slice in pitch for EEPAA
eepaa_1_dat = []
eepaa_2_dat = []
correlator_dat = []

select = int(input('\n\neepaa 1 or eepaa 2 or correlator ?\n (1 or 2 or 3)\n Ans: '))

if select ==1 or select == 2:
    w_pitch = int(input('\n\nWhich slice in Pitch?\n (0 through 20) \n Ans: '))
    # Get the slice in pitch for plotting
    maximum1 = 0
    maximum2 = 0
    for i in range(t_start, t_end):
        eepaa_1_dat.append([counts_1[i][w_pitch][j] for j in range(49)])


    for i in range(t_start, t_end):
        eepaa_2_dat.append([counts_2[i][w_pitch][j] for j in range(49)])

    # -------

elif select == 3:
    w_pitch = int(input('\n\nWhich Phase Bin?\n (0 through 15) \n Ans: '))

    maximum3 = 0
    #get the slice in phase bin
    for i in range(t_start_correlator,t_end_correlator):
        correlator_dat.append([correlators[i][w_pitch][j] for j in range(8)])







# Get the data ready for plotting
if select == 1 or select == 2:
    eepaa_ar = [np.array(eepaa_1_dat),np.array(eepaa_2_dat)]
    counts = [counts_1,counts_2]
    epochs = [epoch_1[t_start:t_end],epoch_2[t_start:t_end]]
    energies = [energies_1,energies_2]
    maxs = [maximum1,maximum2]
    plot_energies = energies[select-1]
    epoch =epochs[select-1]
    plot_eepaa = eepaa_ar[select - 1]
    plot_max = maxs[select-1]
    plot_titles = ['EEPAA_1 Counts','EEPAA_2 Counts']
    plot_pitch = str(pitches_1[w_pitch])
    plot_epoch = cdflib.epochs.CDFepoch.to_datetime(epoch)
    # Make the EEPAA plot (X,Y,Z) coordiantes
    X = plot_epoch
    Y = plot_energies
    Z = plot_eepaa
    Xm,Ym = np.meshgrid(X,Y)

elif select == 3:
    plot_titles = 'Correlator data'
    plot_bin = str(phase_bins[w_pitch])
    X = cdflib.epochs.CDFepoch.to_datetime(epoch_correlator[t_start_correlator:t_end_correlator])
    Y = energy_correlator
    Z = np.array(correlator_dat)
    Xm,Ym = np.meshgrid(X,Y)
    plot_epoch = X




#%%%%%%%%%%%%%%%%%%
# PLOTTING THE DATA
#%%%%%%%%%%%%%%%%%%

Zmax = Z.max()
Zmin = Z.min()

fig,ax0 = plt.subplots()
fig.autofmt_xdate() #Bend the epoch times sideways
fig.set_size_inches(xinches,yinches)

# ----------------
# Spectrogram Plot
# ----------------


if Z.min() ==0:
    Zmin = 1
else:
    Zmin = Z.min()

#EEPAA
if select == 1 or select == 2:
    ax0.set_title(plot_titles[select - 1] + '\n' + 'Pitch ' + plot_pitch + ' deg')

    if colorbarlog:
        im = ax0.pcolormesh(Xm,Ym,np.transpose( Z),norm=colors.LogNorm(vmin=Zmin,vmax=Zmax), shading='auto',cmap='jet')
    else:
        im = ax0.pcolormesh(Xm, Ym, np.transpose(Z), shading='auto', vmin=Zmin, vmax=Zmax,cmap='jet')

    cbar = fig.colorbar(im,ax=ax0)
    cbar.set_label('Counts')
    ratio = round((t_end - t_start) / num_of_ticks)
    ticks = [plot_epoch[i] for i in range(0, (t_end - t_start)) if i % ratio == 0]

    if logy:
        plt.yscale('log')
    ax0.set_ylabel('Energy (eV)')
    ax0.set_title('EEPAA '+ str(select) + ' Counts Data & Integrated PSDs \n' + 'Pitch: ' + plot_pitch + 'deg\n'  )

    if tight_layout:
        plt.tight_layout(-0.3)

#Correlator
elif select == 3:

    if colorbarlog:
        im = ax0.pcolormesh(Xm,Ym,np.transpose( Z),norm=colors.LogNorm(vmin=Zmin,vmax=Zmax), shading='auto',cmap='jet')
    else:
        im = ax0.pcolormesh(Xm, Ym, np.transpose(Z), shading='auto', vmin=Zmin, vmax=Zmax,cmap='jet')

    cbar = fig.colorbar(im, ax=ax0)
    cbar.set_label('Correlators (#)')
    ratio = round((t_end_correlator - t_start_correlator) / num_of_ticks)
    ticks = [plot_epoch[i] for i in range(0, (t_end - t_start)) if i % ratio == 0]

    if logy:
        plt.yscale('log')
    ax0.set_ylabel('Energy (eV)')
    ax0.set_title( 'Correlator Data & Integrated PSDs \n' + 'Phase Bin: ' + plot_bin)

    if tight_layout:
        plt.tight_layout(-0.3)







# -----------------------
# SPECTRAL WAVE PLOT DATA
# -----------------------

ax1 = ax0.twinx()
for freq in which_freq:
    ax1.plot(dat_times_EEPAA,dat_integrate[freq],color=line_colors[freq])



# plt.xticks([tick for tick in dat_times_EEPAA if dat_times_EEPAA.index(tick)%num_of_ticks == 0])# Set the tick resolution on X-scale
# ax1.xaxis.set_major_formatter(myFmt) # Set the format of the xticks to myfmt

# Create legend for plot
which_head = [dat_integrate_headers[i] for i in which_freq]
which_title = ''
for header in which_head:
    which_title += header + ', '
ax1.legend(which_head)

plt.ylim(ylim_low,ylim_high)
ax1.set_xlabel('epoch')
ax1.set_ylabel('Spectral Power [(mV/m)$^{2}$]')

ax0.xaxis.set_major_formatter(myFmt)
ax0.set_xticks(ticks)  # Set the tick resolution on X-scale
plt.show()







