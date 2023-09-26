
def groupings(wRocket,SSA_window_Size):

    # the first component is the original data
    group = [[i for i in range(SSA_window_Size*3)]]


    # --- High Flyer Groupings ---
    if wRocket == 4:

        # --- WINDOWSIZE = 401 groupings ---
        if SSA_window_Size == 401:
            grouping = [
                [0, 1, 2, 3, 4, 5, 6, 7, 14, 15, 16, 17] + [8, 9] + [12, 13] + [10, 11],
                [i for i in range(26, 75)], # my "noise" data (groupings above 100 are truly useless)
                [i for i in range(75, SSA_window_Size * 3)],
                [18, 19, 20, 21, 23, 24, 25]  # THE GOOD ONE + a little noise/signal
            ]

        # --- WINDOWSIZE = 251 groupings ---
        elif SSA_window_Size == 251:
            # [10, 11] may also be noise
            grouping = [
                [0, 1],
                [2, 3],
                [4, 5],
                [6, 7],
                [8, 9],
                [i for i in range(10, SSA_window_Size*3)]
            ]
        elif SSA_window_Size == 351:
            grouping = [
                [0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15],
                [10, 11],
                [8, 9],
                [i for i in range(16, 50)] + [i for i in range(50, SSA_window_Size * 3)],
                [i for i in range(50, SSA_window_Size * 3)]
            ]
        elif SSA_window_Size == 451:
            # is [12, 13] signal or noise?
            # is [8,9] signal or noise?
            grouping = [
                [0, 1, 2, 3, 4 , 5, 6, 7, 10, 11, 14, 15, 16, 17 ] + [12, 13] + [40, 41, 46],
                [i for i in range(18, 40)] + [42, 43, 44, 45, 47, 48] + [i for i in range(49, 54)], # The good signal
                [i for i in range(55, 3*SSA_window_Size)],
                [i for i in range(18, 40)] + [42, 43, 44, 45] + [i for i in range(47, 54)] + [i for i in range(55, 3*SSA_window_Size)]
            ]
        elif SSA_window_Size == 501:

            # identify the noise components
            grouping = [
                [i for i in range(16)] + [17, 18, 19, 43, 44, 45, 47, 48, 49, 62, 63, 64, 79, 80, 81, 85,
                                          86, 88, 89, 92, 93, 96, 97, 98, 99, 100, 101, 102, 105, 106, 110,
                                          111, 113, 114, 117, 118, 119, 120, 123, 124, 128, 129],
            ]

            # investigate noise
            # grouping += [[i] for i in range(131)]
            grouping += [ [i for i in range(131, 3*SSA_window_Size)] ]

            # get the good stuff
            goodStuff = []
            for i in range(130 + 1):
                if i not in grouping[0]:
                    goodStuff.append(i)

            grouping += [goodStuff]



        elif SSA_window_Size == 751:
            # is [10, 11] signal or noise?
            # is [4, 5] signal or noise?
            grouping = [
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21],
                [12, 13] + [i for i in range(22, 3*SSA_window_Size)]
            ]
        else:
            grouping = [
                [i] for i in range(20)
            ]

    # --- LF Groupings ---
    elif wRocket == 5:
        grouping = [
            [0, 1],
            [2, 3],
            [4, 5],
            [i for i in range(6, SSA_window_Size * 3)]
        ]



    return group + grouping