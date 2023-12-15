
def timeWindow(wTargetTimes,wRocket):
    import datetime as dt
    # --- use the SSA components to determine which
    if wTargetTimes == 0: # ALFVEN REGION
        scienceRegions = [
            [dt.datetime(2022, 11, 20, 17, 24, 25, 000000), dt.datetime(2022, 11, 20, 17, 25, 50, 000000)],
            [dt.datetime(2022, 11, 20, 17, 24, 25, 000000), dt.datetime(2022, 11, 20, 17, 25, 50, 000000)]]  # ACTUAL WINDOW for WL 501 High Flyer
    else:
        if wTargetTimes == 1: # KENTON's ZOOMED REGION
            scienceRegions = [
                [dt.datetime(2022, 11, 20, 17, 24, 30, 000000), dt.datetime(2022, 11, 20, 17, 25, 00, 000000)],
                [dt.datetime(2022, 11, 20, 17, 24, 30, 000000), dt.datetime(2022, 11, 20, 17, 25, 00, 000000)]]  # ACTUAL WINDOW for WL 501 High Flyer
        else: # THE WHOLE FLIGHT
            scienceRegions = [
                [dt.datetime(2022, 11, 20, 17, 22, 00, 000000), dt.datetime(2022, 11, 20, 17, 29, 00, 000000)],
                [dt.datetime(2022, 11, 20, 17, 24, 00, 000000), dt.datetime(2022, 11, 20, 17, 27, 00, 000000)]
            ]
    return scienceRegions[wRocket-4]



def groupings(wRocket,SSA_window_Size,wUseData):

    # the first component is the original data
    group = [[i for i in range(SSA_window_Size*3)]]


    # --- High Flyer Groupings ---
    if wRocket == 4:
        if SSA_window_Size == 501:

            if wUseData == 0: # THE ALFVEN RGION
                # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
                grouping = [
                    # [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 46, 47, 48, 49, 68, 69, 76,
                    #  79, 82, 83, 85, 86, 88, 89, 94, 95, 98, 99, 100, 101, 102, 103, 104, 106, 108,109, 110, 111, 112, 113, 117,
                    # 119, 120, 121, 122,125, 128, 129, 130]
                    [0,1,2,3,4,5,7,8,9]
                ]

                # investigate other components
                grouping += [[i] for i in range(10,16)]

                # the "noise" data
                limit = 250
                grouping += [[i for i in range(limit, 3 * SSA_window_Size)]]

                # # get the good stuff
                goodStuff = []
                for i in range(limit): # should be 130
                    if i not in grouping[0]:
                        goodStuff.append(i)

                grouping += [goodStuff]

            else:
                # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
                grouping = [
                    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 36, 37, 44, 45, 55, 56]
                ]

                # investigate other components
                grouping += [[i] for i in range(60, 70)]

                # the "noise" data
                # grouping += [[i for i in range(130, 3 * SSA_window_Size)]]

                # get the good stuff
                goodStuff = []
                for i in range(60): # should be 130
                    if i not in grouping[0]:
                        goodStuff.append(i)

                grouping += [goodStuff]
        elif SSA_window_Size == 1201:
            # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
            grouping = [
                [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                 20,21,22,23,24,25,26,27,30,31,32,33,36,37,59,60,
                 87,88,89,101,102,123,124,126,127,128,138,139,149,
                 150,151,157,158,166,167,168,171,172,173,174,178,179,
                 183,184,193,194,195,197,207,208,209,210,220,221,226,227,
                 228,229,230,237,242,243,244,245,246,247,248
                 ]
            ]

            # investigate other components
            grouping += [[i] for i in range(240,250)]

            # the "noise" data
            limit = 1100
            grouping += [[i for i in range(limit, 3 * SSA_window_Size)]]

            # get the good stuff
            goodStuff = []
            for i in range(limit):  # should be 130
                if i not in grouping[0]:
                    goodStuff.append(i)

            grouping += [goodStuff]
        else:
            grouping = [
                [i] for i in range(10)
            ]

    # --- LF Groupings ---
    elif wRocket == 5:

        if SSA_window_Size == 501:

            if wUseData==0: # THE ALFVEN REGION

                # [dt.datetime(2022, 11, 20, 17, 24, 25, 000000), dt.datetime(2022, 11, 20, 17, 25, 18, 000000)]  # ACTUAL WINDOW for WL 501 Low Flyer
                # identify the harmonic components
                grouping = [
                    [0,1,2,3,4,7,8,9]
                ]

                # investigate other components
                grouping += [[i] for i in range(25,35)]

                # show the "noise"
                limit = 150
                grouping += [[i for i in range(limit,3*SSA_window_Size)]]

                # get the good stuff
                goodStuff = []
                for i in range(limit):  # should be 150
                    if i not in grouping[0]:
                        goodStuff.append(i)

                grouping += [goodStuff]
            else:
                # grouping = [
                #     [i] for i in range(20)
                # ]

                grouping = [
                    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
                ]

                # investigate other components
                grouping += [[i] for i in range(15,21)]

                # show the "noise"
                limit = 400
                grouping += [[i for i in range(limit, 3 * SSA_window_Size)]]

                # get the good stuff
                goodStuff = []
                for i in range(limit):  # should be 150
                    if i not in grouping[0]:
                        goodStuff.append(i)

                grouping += [goodStuff]
        elif SSA_window_Size == 1001:
            grouping = [
                [0,1,2,3,4,5, 6, 7, 8, 9, 14, 15, 18, 19,20, 26, 29, 30, 31]
            ]

            # investigate other components
            # grouping += [[i] for i in range(36,45)]

            # show the "noise"
            limit = 250
            grouping += [[i for i in range(limit, 3 * SSA_window_Size)]]

            # get the good stuff
            goodStuff = []
            for i in range(limit):  # should be 150
                if i not in grouping[0]:
                    goodStuff.append(i)

            grouping += [goodStuff]

        else:
            grouping = [
                [i] for i in range(10)
            ]

    return group + grouping