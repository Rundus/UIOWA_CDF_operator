# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# -------- CAPERII Magnetometer,Pitch Angle and Distribution Funct. output/sorting Code-------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Time the program
import time
start_time = time.time()


# --------------------
# Run the main program
# --------------------

# import Chi_square_corrections
# import Calc_Differential_Flux
import Magnetometer_rotations
import Sorting_DiffFlux
import Dist_Funcs
# import J_parallel
import Build_output_File


print("--- %s seconds for full program---" % (time.time() - start_time) )