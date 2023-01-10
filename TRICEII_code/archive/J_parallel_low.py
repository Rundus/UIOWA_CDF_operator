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

from files import dist_file_low,EEPAA_file_low,root
from Variables import Energies,EPOCH_low,EEPAA_file_low_info
from Variables import pitch_output_low,zvars_low,J_low_tag,dist_jpar_low


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

J_par_low_output = cdfwrite.CDF(root  + J_low_tag +'.cdf',cdf_spec=EEPAA_file_low_info,delete=True)





Js_low = []

print('Calculating J_parallel for LOW: ', end='')
low_time = time.time()

bins = [tme,tme+1,tme+100,tme+1000]
#

for i in range(len(EPOCH_low)):
    v_perp_low = []
    v_par_low = []
    dist_low_1Dlist = []
    # -----------------
    # Organize the data
    # -----------------
    for ptch in range(21):
        if ptch == 0:
            # need to make the pitch value -10
            for engy in range(49):
                v_par_low.append(math.cos(math.radians(pitch_output_low[i][ptch][engy] - 20)) * (((2 * Energies[engy] * q) / m) ** (1 / 2)))
                v_perp_low.append(math.sin(math.radians(pitch_output_low[i][ptch][engy] - 20)) * (((2 * Energies[engy] * q) / m) ** (1 / 2)))
                dist_low_1Dlist.append(dist_jpar_low[i][ptch][engy])

        elif ptch == 20:
            # need to make the pitch value -10
            for engy in range(49):
                v_par_low.append(math.cos(math.radians(pitch_output_low[i][ptch][engy] + 20)) * (((2 * Energies[engy] * q) / m) ** (1 / 2)))
                v_perp_low.append(math.sin(math.radians(pitch_output_low[i][ptch][engy] + 20)) * (((2 * Energies[engy] * q) / m) ** (1 / 2)))
                dist_low_1Dlist.append(dist_jpar_low[i][ptch][engy])
        else:
            for engy in range(49):
                v_par_low.append(math.cos(math.radians(pitch_output_low[i][ptch][engy])) * (((2 * Energies[engy] * q) / m) ** (1 / 2)))
                v_perp_low.append(math.sin(math.radians(pitch_output_low[i][ptch][engy])) * (((2 * Energies[engy] * q) / m) ** (1 / 2)))
                dist_low_1Dlist.append(dist_jpar_low[i][ptch][engy])


    # -----------------
    # Interpolation low
    # -----------------
    x = np.array(v_perp_low)
    y = np.array(v_par_low)
    z = np.array(dist_low_1Dlist)
    print(i)
    for num in Nvals:
        start = time.time()
        X = np.linspace(min(v_perp_low), max(v_perp_low), num)
        Y = np.linspace(min(v_par_low), max(v_par_low), num)
        X, Y = np.meshgrid(X, Y, indexing='xy')
        interp = LinearNDInterpolator(list(zip(x, y)), z, fill_value=fillval)
        Z = interp(X, Y)
        # plt.pcolormesh(X, Y, Z, shading='auto')
        # plt.plot(x,y)
        # plt.legend()
        # plt.colorbar()
        # plt.xlabel('$v_{\perp}$')
        # plt.ylabel('$v_{\parallel}$')
        # plt.axis("equal")
        # plt.show()

        # -----------------------------
        # Integrate to get First moment
        # -----------------------------

        # METHOD: Trapezoidal Integration
        G = []
        suming = []
        sum1ing = []
        # High
        for k in range(num):
            suming = [(Y[j + 1][k] - Y[j][k]) * (Z[j + 1][k] + Z[j][k]) / 2 for j in range(num - 1)]
            G.append(sum(suming))
        sum1ing = [(X[0][l + 1] - X[0][l]) * (G[l + 1] - G[l]) / 2 for l in range(num - 1)]
        Js_low.append(-1 * sum(sum1ing) * 2 * math.pi * q)



print('Done')


print("--- %s seconds---" % (time.time() - low_time),'\n\n')

#Convet data to format that can write out to file

Js_low_dat = np.ndarray(shape=(len(Js_low)),dtype='float64')




for i in range(len(Js_low)):
    Js_low_dat[i] = Js_low[i]

print(len(Js_low))
print(Js_low)
print(Js_low_dat)

# ---------------
# OUTPUT THE DATA
# ---------------


#LOW
vardata= Js_low_dat

varinfo = EEPAA_file_low.varinq(zvars_low[6])
varinfo['Variable'] = 'Parallel_Current_Density'
varinfo['Num_Dims'] = 0
varinfo['Dim_Sizes'] = []

varattrs_1 = EEPAA_file_low.varattsget(zvars_low[6], expand=True)
varattrs_1['CATDESC'] = 'Parallel_Current_Density'
varattrs_1['FIELDNAM'] = 'Parallel_Current_Density'
varattrs_1['UNITS'] = 'A!N m!U-2!'
varattrs_1['SCALETYP'] = 'log'
varattrs_1['VALIDMIN'] = [min(Js_low),'CDF_FLOAT']
varattrs_1['VALIDMAX'] = [max(Js_low),'CDF_FLOAT']
varattrs_1['LABLAXIS'] = 'Parallel_Current_Density'
#Things that need changing
# varattrs_high['LABL_PTR_1'] = Label_1
# varattrs_high['LABL_PTR_2'] = Label_2
varattrs_1['LABL_PTR_1'] = 'LABL_PTR_1'
varattrs_1['LABL_PTR_2'] = 'LABL_PTR_2'

J_par_low_output.write_var(varinfo,var_attrs=varattrs_1,var_data=vardata)







j_par_time = time.time()

print("--- %s seconds---" % (time.time() - j_par_time),'\n\n')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dist_file_low.close()
EEPAA_file_low.close()



print("--- %s seconds for J_parallel_low---" % (time.time() - start_time) ,'\n')



