# --- dispersionAttributes.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Location to store the information for the dispersion times

import numpy as np
import datetime as dt


###############################################################
# --- PLOT 1: BOTH FLYERS, ALT VS LAT (w/ CONNECTING LINES) ---
###############################################################
def boxRemove(Data, EnergyMin, EnergyMax, TimeMin, TimeMax, Energy, Time):
    newData = Data
    for tme in range(len(Data)):
        for engy in range(len(Data[tme])):
            if Energy[engy] > EnergyMin and Energy[engy] < EnergyMax and Time[tme] > TimeMin and Time[tme] < TimeMax:
                newData[tme][engy] = 0
    return Data

def diagonalRemove(Data, indexSum ,TimeMin, TimeMax, Energy, Time, upper):
    newData = Data
    for tme in range(len(Data)):
        for engy in range(len(Data[tme])):
            if upper:
                if Time[tme] > TimeMin and Time[tme] < TimeMax and engy - tme <= indexSum:
                    newData[tme][engy] = 0
            else:
                if Time[tme] > TimeMin and Time[tme] < TimeMax and engy - tme >= indexSum:
                    newData[tme][engy] = 0
    return Data

class dispersionAttributes:

    # taken from the 10deg bin on eepaa_fullCal. These are the most noticable
    keyDispersionTimes = [
        [dt.datetime(2022, 11, 20, 17, 24, 56, 000000), dt.datetime(2022, 11, 20, 17, 24, 57, 600000)],# s1
        [dt.datetime(2022, 11, 20, 17, 24, 57, 710000), dt.datetime(2022, 11, 20, 17, 24, 59, 186000)],# s2
        [dt.datetime(2022, 11, 20, 17, 24, 58, 558000), dt.datetime(2022, 11, 20, 17, 24, 59, 958000)],# s3
        [dt.datetime(2022, 11, 20, 17, 24, 59, 965000), dt.datetime(2022, 11, 20, 17, 25, 00, 758000)],# s4  right before Inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 00, 501000), dt.datetime(2022, 11, 20, 17, 25, 1, 513000)], # s5 inside inverted V's left edge
        [dt.datetime(2022, 11, 20, 17, 25, 00, 909000), dt.datetime(2022, 11, 20, 17, 25, 2, 16000)], # s6 under V, it's faint
        [dt.datetime(2022, 11, 20, 17, 25, 3, 500000), dt.datetime(2022, 11, 20, 17, 25, 4, 550000)], # s7 Faint signature after small inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 3, 800000), dt.datetime(2022, 11, 20, 17, 25, 5, 200000)], # s8 in choatic region
        [dt.datetime(2022, 11, 20, 17, 25, 5, 450000), dt.datetime(2022, 11, 20, 17, 25, 5, 920000)], # s9 VERTICAL BAR that doesn't look like it's in an inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 5, 806000), dt.datetime(2022, 11, 20, 17, 25, 6, 506000)], # s10 Last dispersive feature in the chaotic region
        [dt.datetime(2022, 11, 20, 17, 25, 6, 562000), dt.datetime(2022, 11, 20, 17, 25, 6, 945000)], # s11 small feature right before Very Strong feature outside inverted V region
        [dt.datetime(2022, 11, 20, 17, 25, 6, 760000), dt.datetime(2022, 11, 20, 17, 25, 7, 606000)], # s12 Very Strong Wave-like feature on north end of Inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 7, 641000), dt.datetime(2022, 11, 20, 17, 25, 8, 406000)], # s13 Smaller feature, similar to s1
        [dt.datetime(2022, 11, 20, 17, 25, 8, 434000), dt.datetime(2022, 11, 20, 17, 25, 8, 946000)] # s14 Tiny Feature, similar to s10 and s11
    ]

    # the overhead noise on the dispersion features shouldn't be included. Here we have the energy indicies for the ~max energy the feature appears to be at when looking at pitch 0deg
    # NOTE: these indicies assume you're using the full 49 length energy array
    keyDispersionEnergyLimits = np.array([
        24,  # s1
        27,  # s2
        26,  # s3
        20,  # s4
        19,  # s5
        18,  # s6
        25,  # s7
        24,  # s8
        23,  # s9
        22,  # s10
        24,  # s11
        22,  # s12
        27,  # s13
        31  # s14
        ])

    additionalDispersionTimes = [
        [dt.datetime(2022, 11, 20, 17, 24, 31, 193000), dt.datetime(2022, 11, 20, 17, 24, 32, 656000)], # sa1 near the beginning of the flight
        [dt.datetime(2022, 11, 20, 17, 25, 23, 500000), dt.datetime(2022, 11, 20, 17, 24, 24, 461000)], # sa2 at the beginnig edge of the largest inverted V
        [dt.datetime(2022, 11, 20, 17, 26, 10, 107000), dt.datetime(2022, 11, 20, 17, 26, 11, 159000)],  # sa3 near the south edge of inverted V, almost undernearth it. NORTH side of Big Inverted V
        [dt.datetime(2022, 11, 20, 17, 26, 20, 603000), dt.datetime(2022, 11, 20, 17, 26, 21, 308000)], # sa4 underneath an inverted  V on the northernmost edge of the big inverted V
        [dt.datetime(2022, 11, 20, 17, 26, 20, 603000), dt.datetime(2022, 11, 20, 17, 26, 21, 308000)]
        # sa3 underneath an inverted  V on the northernmost edge of the big inverted V
    ]


    ##############################
    # --- CLEAN THE SIGNATURES ---
    ##############################

    # I need to isolate the dispersion signatures so I can T.O.F. fit them. To do this,
    # have to manually isolate them. Lets create a special set of functions for each dispersion to do just this.


    @staticmethod
    def cleanS1(inputData,Energy,Time): # doesn't need anything
        return inputData

    @staticmethod
    def cleanS2(inputData,Energy,Time): # remove data in top right
        newData = inputData

        # --- BOX REMOVE ---
        EnergyMin, EnergyMax,TimeMin,TimeMax = 90, 200, 1, 1.5
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax,Energy,Time)

        return newData

    @staticmethod
    def cleanS3(inputData, Energy, Time):  # doesn't need anything
        newData = inputData

        # --- BOX REMOVE ---
        # remove lower left
        EnergyMin, EnergyMax, TimeMin, TimeMax = 0, 70, -0.1, 0.58
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax,Energy,Time)

        # --- BOX REMOVE ---
        # upper right
        EnergyMin, EnergyMax, TimeMin, TimeMax = 95, 300, 1, 2
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax, Energy, Time)

        return newData

    @staticmethod
    def cleanS4(inputData, Energy, Time):  # doesn't need anything
        newData = inputData

        # --- BOX REMOVE ---
        # remove past 0.7s
        EnergyMin, EnergyMax, TimeMin, TimeMax = 128, 1000, 0.6, 2
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax, Energy, Time)

        # --- BOX REMOVE ---
        # remove bottom left
        EnergyMin, EnergyMax, TimeMin, TimeMax = 0, 90, 0, 0.25
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax, Energy, Time)

        return newData

    @staticmethod
    def cleanS5(inputData, Energy, Time):  # doesn't need anything
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 7, -0.1, 0.3
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove upper right (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = -5, 0.3, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # remove upper right (FINE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = 0, 0.4, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # remove upper right (FINE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = 1, 0.48, 0.7
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS6(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 5, -0.1, 0.6
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = -2, 0.3, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS7(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 9, -0.1, 0.35
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove upper right (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = -3, 0.45, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)


        return newData

    @staticmethod
    def cleanS8(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 3, -0.1, 0.7
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove upper right (COARSE)
        # upper = True
        # diagonalCounter, TimeMin, TimeMax = -5, 0.6, 2
        # newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove upper right (Fine)
        upper = True
        diagonalCounter, TimeMin, TimeMax = -6, 0.4, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS9(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 13, -0.1, 0.2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove upper right (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = -5, 0.15, 1
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData


    @staticmethod
    def cleanS10(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 11, -0.1, 0.4
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (FINE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 6, -0.1, 0.15
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (FINE FINE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 8, 0.1, 0.3
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = -3, 0.3, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS11(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 9, -0.1, 0.16
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (FINE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 4, -0.1, 0.06
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove upper right (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = -2, 0.2, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS12(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 9, -0.1, 0.25
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax =-3, 0.55, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS13(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 5, -0.1, 0.2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = -5, 0.6, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS14(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 4, -0.1, 0.15
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = -4, 0.25, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData


    # Store the functions in a datavariable
    isolationFunctions = {'s1': cleanS1,
                          's2': cleanS2,
                          's3': cleanS3,
                          's4': cleanS4,
                          's5': cleanS5,
                          's6': cleanS6,
                          's7': cleanS7,
                          's8': cleanS8,
                          's9': cleanS9,
                          's10': cleanS10,
                          's11': cleanS11,
                          's12': cleanS12,
                          's13': cleanS13,
                          's14': cleanS14}