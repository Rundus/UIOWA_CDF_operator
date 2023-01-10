# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------- TRICE Magnetometer,Pitch Angle and Distribution Funct. output/sorting Code--------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Time the program
import time
start_time = time.time()

# Run the main program
import Sun_spike_noise_collection
import Sun_spike_noise_statistics
import Sun_spike_removal
import Magnetometer_rotations
import Chi_square_pad_pair_investigation
import Calc_Differential_Flux
import Sorting_DiffFlux
import Dist_Funcs
# import J_parallel_high
# import J_parallel_low
import Build_output_File


print("--- %s seconds for full program---" % (time.time() - start_time) )