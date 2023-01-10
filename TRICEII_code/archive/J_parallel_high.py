# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE Parallel Current Density J --------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np,time,math
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from cdflib import cdfwrite


print('importing variables: ', end='')

start_time = time.time()

# ------------
# Import Files
# ------------

from files import dist_file_high,EEPAA_file_high,root
from Variables import Energies,EPOCH_High,EEPAA_file_high_info,dist_jpar_high
from Variables import pitch_output_high,zvars_high,J_high_tag,pitch


print('Done')



# ######################################
# Calculate the parallel current density
# ######################################



# Toggles-------------------------------------
tme = 16500  # Which time slice
N = 10      # how many points on the mesh grid
fillval = 0
# Nvals = [1500] # ideal value for semi-accuacy vs  speed of computation is Nvals = [1500]...ACTUALLY it could be N=500
Nvals = [500]
# Toggles-------------------------------------

# ----------------------------
# Get parallel/Perp velocities
# ----------------------------

#Constants
m = 9.11 * 10**(-31)
q = 1.60217662 * 10**(-19)
e = -1*q


# ------------
# Output Files
# ------------

J_par_high_output = cdfwrite.CDF(root  + J_high_tag + '.cdf',cdf_spec=EEPAA_file_high_info,delete=True)



print('Calculating J_parallel for HIGH: ', end='')
high_time = time.time()
Js_high = []

for i in range(len(EPOCH_High)):
    print(i)
    v_par_high = []
    v_perp_high = []
    dist_high_1Dlist = []
    # -----------------
    # Organize the data
    # -----------------

    for ptch in range(21):
        if ptch == 0:
            # need to make the pitch value -10
            for engy in range(49):
                v_par_high.append(math.cos(     math.radians( pitch_output_high[i][ptch][engy] - 20 )) * ( ((2 * Energies[engy] * q) / m) ** (1 / 2)))
                v_perp_high.append(math.sin(math.radians(pitch_output_high[i][ptch][engy] - 20)) * (((2 * Energies[engy] * q) / m) ** (1 / 2)))
                dist_high_1Dlist.append(dist_jpar_high[i][ptch][engy])

        elif ptch == 20:
            # need to make the pitch value -10
            for engy in range(49):
                v_par_high.append(math.cos(     math.radians( pitch_output_high[i][ptch][engy] + 20 )) * ( ((2 * Energies[engy] * q) / m) ** (1 / 2)))
                v_perp_high.append(math.sin(math.radians(pitch_output_high[i][ptch][engy] + 20)) * (((2 * Energies[engy] * q) / m) ** (1 / 2)))
                dist_high_1Dlist.append(dist_jpar_high[i][ptch][engy])
        else:
            for engy in range(49):
                v_par_high.append( math.cos(math.radians( pitch_output_high[i][ptch][engy] )) * ( ((2 * Energies[engy] * q) / m) ** (1 / 2)))
                v_perp_high.append(math.sin(math.radians(pitch_output_high[i][ptch][engy])) * (((2 * Energies[engy] * q) / m) ** (1 / 2)) )
                dist_high_1Dlist.append(dist_jpar_high[i][ptch][engy])



    # ------------------
    # Interpolation high
    # ------------------

    # times = []
    x= np.array(v_perp_high)
    y= np.array(v_par_high)
    z= np.array(dist_high_1Dlist)
    #
    for num in Nvals:
        start = time.time()
        X = np.linspace(min(v_perp_high),max(v_perp_high),num)
        Y = np.linspace(min(v_par_high),max(v_par_high),num)
        X,Y = np.meshgrid(X,Y, indexing = 'xy')
        interp = LinearNDInterpolator(list(zip(x,y)),z,fill_value = fillval)
        Z = interp(X,Y)
    #     plt.pcolormesh(X, Y, Z, shading='auto')
    #     plt.plot(x,y)
    #     plt.legend()
    #     plt.colorbar()
    #     plt.xlabel('$v_{\perp}$')
    #     plt.ylabel('$v_{\parallel}$')
    #     plt.axis("equal")
    #     plt.show()

        # -----------------------------
        # Integrate to get First moment
        # -----------------------------
        #
        # METHOD: Trapezoidal Integration
        G = []
        J = []

        # The variable data looks like:
        # x,y describe the coordinates on the graph of plot(X,Y,Z)
        # linspace'd X,Y data looks like X[0][x], Y[y][0] because there are no variation across the first index in X or the second index in Y (they're all the same because of meshgrind)
        # The distribution value should be worked as Z[y][x]
        #
        # High
        for k in range(num):
            suming = [(Y[j+1][k] - Y[j][k]) * (Z[j+1][k] + Z[j][k])/2 for j in range(num-1)]
            G.append(sum(suming))
        sum1ing = [(X[0][l+1] - X[0][l]) * (G[l + 1] - G[l])/2 for l in range(num-1)]
        Js_high.append( -1*sum(sum1ing)*2*math.pi*q)


print('Done')

print("--- %s seconds---" % (time.time() - high_time),'\n\n')


#Convet data to format that can write out to file
Js_high_dat = np.ndarray(shape=(len(Js_high)),dtype='float64')


for i in range(len(Js_high)):
    Js_high_dat[i] = Js_high[i]



# ---------------
# OUTPUT THE DATA
# ---------------

#HIGH
vardata= Js_high_dat
varinfo = EEPAA_file_high.varinq(zvars_high[6])
varinfo['Variable'] = 'Parallel_Current_Density'
varinfo['Num_Dims'] = 0
varinfo['Dim_Sizes'] = []


varattrs_1 = EEPAA_file_high.varattsget(zvars_high[6], expand=True)
varattrs_1['CATDESC'] = 'Parallel_Current_Density'
varattrs_1['FIELDNAM'] = 'Parallel_Current_Density'
varattrs_1['UNITS'] = 'A!N m!U-2!'
varattrs_1['SCALETYP'] = 'log'
varattrs_1['VALIDMIN'] = [min(Js_high),'CDF_FLOAT']
varattrs_1['VALIDMAX'] = [max(Js_high),'CDF_FLOAT']
varattrs_1['LABLAXIS'] = 'Parallel_Current_Density'
#Things that need changing
# varattrs_high['LABL_PTR_1'] = Label_1
# varattrs_high['LABL_PTR_2'] = Label_2
varattrs_1['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_1['LABL_PTR_2'] = 'LABL_PTR_2'

J_par_high_output.write_var(varinfo,var_attrs=varattrs_1,var_data=vardata)



j_par_time = time.time()

print("--- %s seconds---" % (time.time() - j_par_time),'\n\n')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dist_file_high.close()
EEPAA_file_high.close()




print("--- %s seconds for J_parallel_high---" % (time.time() - start_time) ,'\n')






# -----------------
# CODE ANALYSIS
# -----------------
#
# print('\n')
# print(Js)
# print(Nvals)
#
# plt.title('J_{$\parallel$} Code Analysis')
# plt.subplot(3,1,1)
# plt.plot(Nvals,Js,'ok')
#
# plt.ylabel('Curent Density (J)')
#
# plt.subplot(3,1,2)
# plt.plot(Nvals,times,'ok')
# plt.ylabel('time (s)')
#
# times_tot = []
# for tme in times:
#     times_tot.append( ((tme * len(EPOCH_High))/60)  )
#
# j_par_time = time.time()
#
# plt.subplot(3,1,3)
# plt.plot(Nvals,times_tot,'ok')
# plt.xlabel('Number of grindpoints, NxN (N) ')
# plt.ylabel('Total time to run code (min)')
# plt.show()

