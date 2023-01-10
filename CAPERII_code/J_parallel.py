# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE Parallel Current Density J --------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



import numpy as np,time,itertools,math
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from cdflib import cdfread,cdfwrite
from iteration_utilities import duplicates


print('importing variables: ', end='')

start_time = time.time()

# ------------
# Import Files
# ------------

from files import dist_file_high,dist_file_low,EEPAA_file_low,EEPAA_file_high,pitch_actual_high_new,pitch_actual_low_new

# dist_file_high = cdfread.CDF("C:/Users/Rundus/Desktop/TRICE/Dist_high.cdf")
# dist_file_low = cdfread.CDF("C:/Users/Rundus/Desktop/TRICE/Dist_low.cdf")
# EEPAA_file_high = cdfread.CDF("C:/Users/Rundus/Desktop/TRICE/TRICE_52003_l2_eepaa_20181208T082239_v1.1.2.cdf")
# EEPAA_file_low = cdfread.CDF("C:/Users/Rundus/Desktop/TRICE/TRICE_52004_l2_eepaa_20181208T082243_v1.1.2.cdf")
# pitch_actual_high = cdfread.CDF('C:/Users/Rundus/Desktop/TRICE/pitch_actual_high.cdf')
# pitch_actual_low = cdfread.CDF('C:/Users/Rundus/Desktop/TRICE/pitch_actual_low.cdf')

# Get the input file variable information
EEPAA_file_high_info = EEPAA_file_high.cdf_info()
EEPAA_file_low_info = EEPAA_file_low.cdf_info()
dist_high_info = dist_file_high.cdf_info()
dist_low_info = dist_file_low.cdf_info()


# Get the input file zVariables
zvars_dist_high = dist_high_info['zVariables']
zvars_dist_low = dist_low_info['zVariables']
zvars_high = EEPAA_file_high_info['zVariables']
zvars_low = EEPAA_file_low_info['zVariables']


# ---------
# VARIABLES
# ---------
dist_high = dist_file_high.varget('Distribution_Function')
dist_low = dist_file_low.varget('Distribution_Function')
EPOCH_High = EEPAA_file_high.varget('Epoch')
EPOCH_low = EEPAA_file_low.varget('Epoch')
Energies = EEPAA_file_high.varget('Energy')
pitches = EEPAA_file_high.varget('Pitch_Angle')
pitches_high = pitch_actual_high_new.varget('Pitch_Angle_Actual')
pitches_low = pitch_actual_low_new.varget('Pitch_Angle_Actual')
count_interval_high = EEPAA_file_high.varget('count_interval')* (10**9)
count_interval_low = EEPAA_file_low.varget('count_interval')* (10**9)



ranges_high_avg = [
    range(len(EPOCH_High)),
    range(21),
    range(49)]

ranges_low_avg = [
    range(len(EPOCH_low)),
    range(21),
    range(49)]

print('Done')



# ######################################
# Calculate the parallel current density
# ######################################



# ----------------------------
# Get parallel/Perp velocities
# ----------------------------

m = 9.11 * 10**(-31)
q = 1.60217662 * 10**(-19)
e = -1.60217662 * 10**(-19)


v_par_high=[]
v_par_low=[]
v_perp_high= []
v_perp_low = []
dist_high_list = []
dist_low_list = []


print('Calculating parallel Velocities for HIGH: ', end='')
vel_high_time = time.time()



tme = 16106

#High (time slice)
for ptch in range(21):
    for engy in range(49):
        v_par_high.append( math.cos((math.pi / 180) * pitches_high[tme][ptch][engy]) * ((2 * Energies[engy] * q) / m) ** (1 / 2))
        v_perp_high.append( math.sin((math.pi / 180) * pitches_high[tme][ptch][engy]) * ((2 * Energies[engy] * q) / m) ** (1 / 2))
        v_par_low.append(math.cos((math.pi / 180) * pitches_low[tme][ptch][engy]) * ((2 * Energies[engy] * q) / m) ** (1 / 2))
        v_perp_low.append(math.sin((math.pi / 180) * pitches_low[tme][ptch][engy]) * ((2 * Energies[engy] * q) / m) ** (1 / 2))
        dist_high_list.append(dist_high[tme][ptch][engy])
        dist_low_list.append(dist_low[tme][ptch][engy])



print('Done')


# ------------------
# Interpolation high
# ------------------

N = 150

print('Calculating J_parallel for HIGH: ', end='')
j_par_time = time.time()
X = np.linspace(min(v_perp_high),max(v_perp_high),N)
Y = np.linspace(min(v_par_high),max(v_par_high),N)
z = np.hypot(v_perp_high,v_par_high)
Xm,Ym = np.meshgrid(X,Y)

nancheck = 0

interp = LinearNDInterpolator(list(zip(v_perp_high,v_par_high)),z,fill_value = nancheck)
Z = interp(Xm,Ym)
# plt.pcolormesh(Xm, Ym, Z, shading='auto')
# plt.plot(v_perp_high, v_par_high, "ok", label="input point")
# plt.legend()
# plt.xlabel('$v_{\perp}$')
# plt.ylabel('$v_{\parallel}$')
# plt.colorbar()
# plt.axis("equal")
# plt.show()




# ------------------------------------------
# Code to search for best integration points
# ------------------------------------------

#Calculates the v_perp from the "idea" pitch angles from the pitch angles bins

v_perp_target_high = []
v_perp_target_low = []



for ptch in range(21):
    for engy in range(49):
        v_perp_target_high.append( math.sin((math.pi / 180) * pitches[ptch]) * ((2 * Energies[engy] * q) / m) ** (1 / 2))
        v_perp_target_low.append(math.sin((math.pi / 180) * pitches[ptch]) * ((2 * Energies[engy] * q) / m) ** (1 / 2))



v_perp_target_high = list(set([x for x in v_perp_target_high if v_perp_target_high.count(x) > 1]))
v_perp_target_low = list(set([x for x in v_perp_target_low if v_perp_target_low.count(x) > 1]))

#Finds the indices in the interpolation which correspond to the ideal v_perp

Indexes = []





# -----------------------------
# Integrate to get First moment
# -----------------------------

# METHOD: Trapezoidal Integration
# print(X)
# print(Y)


G = []
J = []

# High
for i in range(N):
    sum = 0
    for j in range(N-1):
        sum += (Y[j+1] - Y[j]) * (Z[i][j+1] + Z[i][j])/2
    G.append(sum)



sum1 = 0
for i in range(N-1):
    sum1 += (X[i+1] - X[i]) * (G[i + 1] - G[i])/2

J = -1*sum1*2*math.pi*q


print('Done')

print("--- %s seconds---" % (time.time() - j_par_time),'\n\n')












#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dist_file_high.close()
dist_file_low.close()
EEPAA_file_high.close()
EEPAA_file_low.close()



print("--- %s seconds for J_parallel---" % (time.time() - start_time) ,'\n')