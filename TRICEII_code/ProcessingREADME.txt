--- NOTES ON THE DATA PROCESSING FOR THE TRICEII DATA ---

-The order the code is to be executed is:

[N/A] csv_to_cdf.attitude.py
[p1] sun_glint_collection.py
[p2] glint_noise_map.py
[p3] magnetometer_pitch&transf.py
[p4] chiSquare_padPair.py
[p5] chiSquare_analysis.py
[p6] chiSquare_corrections.py
[p7] L1_to_L2

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


- The data files are reported by the order of their processing, similar to the L0,L1,L2 style.
  I use the identifier 'p#', where 'p' means processed data and # is the index of the code which produced it.
  For example, 'p5' means it is processed data produced by code #5 in the above list.