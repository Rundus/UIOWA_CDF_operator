# --- dispersionAttributes.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Location to store the information for the dispersion times

import numpy as np
import datetime as dt


###############################################################
# --- PLOT 1: BOTH FLYERS, ALT VS LAT (w/ CONNECTING LINES) ---
###############################################################

class dispersionAttributes:

    # taken from the 10deg bin on eepaa_fullCal. These are the most noticable
    keyDispersionTimes = [
        [dt.datetime(2022, 11, 20, 17, 24, 56, 000000), dt.datetime(2022, 11, 20, 17, 24, 57, 600000)],# s1
        [dt.datetime(2022, 11, 20, 17, 24, 57, 710000), dt.datetime(2022, 11, 20, 17, 24, 59, 186000)],# s2
        [dt.datetime(2022, 11, 20, 17, 24, 58, 558000), dt.datetime(2022, 11, 20, 17, 24, 59, 958000)],# s3
        [dt.datetime(2022, 11, 20, 17, 24, 59, 965000), dt.datetime(2022, 11, 20, 17, 25, 00, 758000)],# s4  right before Inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 00, 501000), dt.datetime(2022, 11, 20, 17, 25, 1, 513000)], # s5 inside inverted V's left edge
        [dt.datetime(2022, 11, 20, 17, 25, 00, 909000), dt.datetime(2022, 11, 20, 17, 25, 2, 16000)],  # s6 under V, it's faint
        [dt.datetime(2022, 11, 20, 17, 25, 3, 800000), dt.datetime(2022, 11, 20, 17, 25, 5, 200000)],  # s7 in choatic region
        [dt.datetime(2022, 11, 20, 17, 25, 5, 457000), dt.datetime(2022, 11, 20, 17, 25, 5, 905000)],   # s8 VERTICAL BAR that doesn't look like it's in an inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 6, 760000), dt.datetime(2022, 11, 20, 17, 25, 7, 606000)]   # s9 High Count Feature that MIGHT be a 1/2 inverted V
    ]

    additionalDispersionTimes = [
        [dt.datetime(2022, 11, 20, 17, 24, 31, 193000), dt.datetime(2022, 11, 20, 17, 24, 32, 656000)], # sa1 near the beginning of the flight
        [dt.datetime(2022, 11, 20, 17, 26, 10, 107000), dt.datetime(2022, 11, 20, 17, 26, 11, 159000)],  # sa2 near the south edge of inverted V, almost undernearth it. NORTH side of Big Inverted V
        [dt.datetime(2022, 11, 20, 17, 26, 20, 603000), dt.datetime(2022, 11, 20, 17, 26, 21, 308000)], # sa3 underneath an inverted  V on the northernmost edge of the big inverted V
        [dt.datetime(2022, 11, 20, 17, 26, 20, 603000), dt.datetime(2022, 11, 20, 17, 26, 21, 308000)]
        # sa3 underneath an inverted  V on the northernmost edge of the big inverted V

    ]