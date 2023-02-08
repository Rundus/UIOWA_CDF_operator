# # --- scratchPaper.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering

import numpy as np


a =np.array( [
    [ #t0
        [0,2,3], #p1
        [4,5,6], #p2
        [7,8,9]  #p3
    ],

    [ #t1
        [1,2,3],
        [4,5,6],
        [7,8,9]
    ]
])

print(a.min())
