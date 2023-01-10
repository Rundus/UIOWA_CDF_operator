# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE CHARACTERIZE SUN SPIKES --------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

from files import attitude_control_file_high,attitude_control_file_low
import matplotlib.pyplot as plt,numpy as np
from Variables import Energies
import matplotlib.colors as colors
import statistics as st
from Variables import roll_reduced_low,roll_reduced_high,roll_epoch_low,roll_epoch_high,noise_vals_low,noise_vals_high,hist_val_engy
from Variables import LOW_ranges,HIGH_ranges,LOW_ranges_inner,HIGH_ranges_inner

noise = [noise_vals_low,noise_vals_high]
counts = [[],[]]
Raw_Roll_dat = [ attitude_control_file_low.varget('Roll (Euler)'), attitude_control_file_high.varget('Roll (Euler)')]
Roll = [[],[]]
EnergiesBin = [[],[]]
PitchA = [[],[]]

for j in range(2):
    ndat = noise[j]
    for i in range(len(ndat)):
        counts[j].append(ndat[i][0])
        Roll[j].append(ndat[i][1])
        EnergiesBin[j].append(ndat[i][2])
        PitchA[j].append(ndat[i][3])


countsa = np.array(counts,dtype='object')
Rolla = np.array(Roll,dtype='object')
Energiesa = np.array(EnergiesBin,dtype='object')
PitchAa = np.array(PitchA,dtype='object')


#--------------------------
# -------PLOTTING----------
#--------------------------

plot_roll_spectrogram = True
calc_hist_noise_reduction = False
plot_rollreduced_spectrogram = False
clean_up_data = True

#----------------
# --- Toggles ---
#----------------
hist_val_threshold = 3
select = 1 #For finding optimal roll ranges
select2 =  0 #for Roll_Reduced

Xdat = np.array(Rolla[select],dtype='float64')
Ydat = np.array(Energies)
Zdat = np.array(countsa[select],dtype='float64')

x = np.array(Rolla[select],dtype='object')
y = np.array(countsa[select],dtype='object')

filenam = ['LOW','HIGH']
xrange = [90,130]
# yrange = [y.min(),100]
xrange = [x.min(),x.max()]
# yrange = [y.min(),y.max()]
yrange = [0,200]
ylabel = 'counts'
zlabel = 'Freq'
xlabel = 'Roll Angle (deg)'
filenam = ['LOW','HIGH']
cbarmax = 100
xtickspacing = 5
no_of_bins = 800
Roll_lower_range = [[-80,-57],[-84,-62]]
Roll_higher_range = [[100,116],[96,118]]

LOW_ranges = [[-82,-57],[97,125]]
HIGH_ranges = [[-93,-59],[92,125]]
LOW_ranges_inner = [[-79,-59],[101,114]]
HIGH_ranges_inner = [[-84,-62],[95,118]]
#--------------------------------------------------------------------


#Sort the roll data
Raw_Roll = Raw_Roll_dat[select]
Raw_Roll_cleaned =np.sort( np.unique(Raw_Roll[1:-1]))


#------------------------
####### PLOTTING ########
#------------------------

##### HISTOGRAM ####

if plot_roll_spectrogram:
    fig,ax = plt.subplots()
    hist, xedges, yedges = np.histogram2d(x, y, bins=no_of_bins, range=[xrange, yrange])
    plt.title(filenam[select]+"\nThreshold: " + str(hist_val_threshold) + '\n Energies: 11486eV-' + str(Energies[hist_val_engy]) + 'eV')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Construct arrays for the anchor positions of the bars.
    boxingy = 0
    boxingx = 0
    xpos, ypos = np.meshgrid(xedges[:-1]+boxingx, yedges[:-1]+boxingy, indexing="ij")
    xposr = xpos.ravel()
    yposr = ypos.ravel()

    if clean_up_data:
        xpos = np.transpose(xpos)
        ypos = np.transpose(ypos)
        hist = np.transpose(hist)

        xpos_temp = []
        ypos_temp = []
        hist_temp = []

        for i in range(len(hist)):
            if list(hist[i]) != list(np.zeros(len(hist))):
                hist_temp.append(hist[i])
                xpos_temp.append(xpos[i])
                ypos_temp.append(ypos[i])

        xpos = np.transpose(np.array(xpos_temp))
        ypos = np.transpose(np.array(ypos_temp))
        hist = np.transpose(np.array(hist_temp))

    # im = ax.pcolormesh(xpos,ypos,hist,cmap='jet',shading='auto',norm=colors.LogNorm(vmin=1,vmax=hist.max()))
    im = ax.pcolormesh(xpos,ypos,hist,cmap='jet',shading='auto',vmin=hist.min(),vmax=cbarmax)

    cbar = fig.colorbar(im , ax = ax)
    cbar.set_label('Freq.')
    plt.xticks(range(round(xposr.min()),round(xposr.max()),xtickspacing),rotation=45)
    plt.show()
#----------------------------------------------------------------------------------------

#------------------------------------------
#Roll data used to eliminate noise in data:
#------------------------------------------

roll_averages = [[],[]]
Xbins = []
selects = [0,1]

if calc_hist_noise_reduction:
    for select in selects:
        x = np.array(Rolla[select], dtype='object')
        y = np.array(countsa[select], dtype='object')
        xrange = [x.min(), x.max()]
        yrange = [y.min(),y.max()]
        hist, xedges, yedges = np.histogram2d(x, y, bins=no_of_bins, range=[xrange, yrange])

        # Construct arrays for the anchor positions of the bars.
        boxingy = 0
        boxingx = 0
        xpos, ypos = np.meshgrid(xedges[:-1] + boxingx, yedges[:-1] + boxingy, indexing="ij")

        if clean_up_data:
            xpos = np.transpose(xpos)
            ypos = np.transpose(ypos)
            hist = np.transpose(hist)

            xpos_temp = []
            ypos_temp = []
            hist_temp = []

            for i in range(len(hist)):
                if list(hist[i]) != list(np.zeros(len(hist))):
                    hist_temp.append(hist[i])
                    xpos_temp.append(xpos[i])
                    ypos_temp.append(ypos[i])

            xpos = np.transpose(np.array(xpos_temp))
            ypos = np.transpose(np.array(ypos_temp))
            hist = np.transpose(np.array(hist_temp))

        Xbins.append( [xpos[j][0] for j in range(len(xpos[0]))])

        for row in hist:
            temp_data = [dat for dat in row if dat != 0]
            if len(temp_data) != 0:
                roll_averages[select].append(sum(temp_data) / len(temp_data))
            else:
                roll_averages[select].append(0)

print(Xbins)
print(roll_averages)


# -------------------------------------------------------
##### USE ROLL  INFORMATION TO PLOT COUNTS VS ENERGY ####
# -------------------------------------------------------

#Create the Roll Reduced variables

fliers = [0,1]
RR_counts = [[],[]]
RR_Roll = [[],[]]
RR_Energies = [[],[]]
RR_Pitch = [[],[]]


for wrockets in fliers:
    Rolls = Rolla[wrockets]
    Engy = Energiesa[wrockets]
    upper = Roll_higher_range[wrockets]
    lower = Roll_lower_range[wrockets]
    for i in range(len(countsa[wrockets])):
        if (Rolls[i] >= upper[0] and Rolls[i] <= upper[1]) or (Rolls[i] >= lower[0] and Rolls[i] <= lower[1]):
            if Engy[i] >= float(Energies[hist_val_engy]):
                RR_counts[wrockets].append(countsa[wrockets][i])
                RR_Roll[wrockets].append(Rolla[wrockets][i])
                RR_Energies[wrockets].append(Energiesa[wrockets][i])
                RR_Pitch[wrockets].append(PitchAa[wrockets][i])

#--------------------
# ---- TOGGLES 2 ----
#--------------------
xRR = np.array(RR_Energies[select2],dtype='object')
yRR= np.array(RR_counts[select2],dtype='object')
# xrange= [xRR.min(),xRR.max()]
xrange= [1500,12000]
# yrange= [yRR.min(),yRR.max()]
yrange= [0,70]
cbarmax = 500
xlabel = 'Energy (eV)'
ylabel = 'Counts'
no_of_bins = 200
#--------------------



if plot_rollreduced_spectrogram:
    fig, ax = plt.subplots()
    hist, xedges, yedges = np.histogram2d(xRR, yRR, bins=no_of_bins, range=[xrange, yrange])

    filenam = ['LOW', 'HIGH']
    plt.title(filenam[select2] + "\nThreshold: " + str(hist_val_threshold) + '\n Roll Angles' + str(Roll_lower_range[select2]) +','+str(Roll_higher_range[select2]) + 'deg')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Construct arrays for the anchor positions of the bars.
    boxingy = 0
    boxingx = 0
    xpos, ypos = np.meshgrid(xedges[:-1] + boxingx, yedges[:-1] + boxingy, indexing="ij")

    if clean_up_data:
        xpos = np.transpose(xpos)
        ypos = np.transpose(ypos)
        hist = np.transpose(hist)

        xpos_temp = []
        ypos_temp = []
        hist_temp = []

        for i in range(len(hist)):
            if list(hist[i]) != list(np.zeros(len(hist))):
                hist_temp.append(hist[i])
                xpos_temp.append(xpos[i])
                ypos_temp.append(ypos[i])

        xpos = np.transpose(np.array(xpos_temp))
        ypos = np.transpose(np.array(ypos_temp))
        hist = np.transpose(np.array(hist_temp))

    im = ax.pcolormesh(xpos, ypos, hist, cmap='jet', shading='auto', vmin=hist.min(), vmax=cbarmax)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Freq.')
    plt.xticks(Energies[0:hist_val_engy +1],rotation=45)
    plt.show()



#--------------------------------
#### COMPUTE DATA STATISTICS ####
#--------------------------------
selects = [0,1]
ranges = [LOW_ranges,HIGH_ranges]
ranges_inner = [LOW_ranges_inner,HIGH_ranges_inner]
Inners_lower = [[],[]]
Inners_higher = [[],[]]
Outters_lower = [[],[]]
Outters_higher = [[],[]]


Energies_inner = [ [[] for i in range(hist_val_engy+1)] ,  [[] for i in range(hist_val_engy+1)]  ]
Energies_outter = [ [[] for i in range(hist_val_engy+1)] ,  [[] for i in range(hist_val_engy+1)]  ]
avgs_inner = [[],[] ]
avgs_outter = [[],[] ]

names = ['LOW','HIGH']
Energylist = list(Energies)

for select in selects:
    wcounts = countsa[select]
    wroll = Rolla[select]
    wengy = Energiesa[select]
    w_ranges = ranges[select]
    w_ranges_inner = ranges_inner[select]

    for i in range(len(wcounts)):
        if i%200 == 0:
            print('Sorting the noise data for ' + names[select] + ': ' +  str(round(100 * (i/len(wcounts)),1)) + '%',end='\r')
        elif i == (len(wcounts)-1):
            print('Sorting the noise data for ' + names[select] + ':' + ' 100%')

        engyindex = Energylist.index(wengy[i])

        if (wroll[i] >= w_ranges_inner[0][0] and wroll[i] <= w_ranges_inner[0][1]):
            Energies_inner[select][engyindex].append(wcounts[i])
        elif (wroll[i] >= w_ranges_inner[1][0] and wroll[i] <= w_ranges_inner[1][1]):
            Energies_inner[select][engyindex].append(wcounts[i])
        elif (wroll[i] >= w_ranges[0][0] and wroll[i] <= w_ranges[0][1]):
            Energies_outter[select][engyindex].append(wcounts[i])
        elif (wroll[i] >= w_ranges[1][0] and wroll[i] <= w_ranges[1][1]):
            Energies_outter[select][engyindex].append(wcounts[i])

    for j in range(len(Energies_inner[select])):
        avgs_inner[select].append([sum(Energies_inner[select][j])/len(Energies_inner[select][j]),st.median(Energies_inner[select][j])] )

    for j in range(len(Energies_outter[select])):
        avgs_outter[select].append([sum(Energies_outter[select][j])/len(Energies_outter[select][j]),st.median(Energies_outter[select][j])])

#
# print(avgs_inner[0])
# print(avgs_inner[1])
# print(avgs_outter[0])
# print(avgs_outter[1])






print('Done')
