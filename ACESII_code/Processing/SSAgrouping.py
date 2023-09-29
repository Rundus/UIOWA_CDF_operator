
def groupings(wRocket,SSA_window_Size,useAlfvenRegion):

    # the first component is the original data
    group = [[i for i in range(SSA_window_Size*3)]]


    # --- High Flyer Groupings ---
    if wRocket == 4:
        if SSA_window_Size == 501:

            if useAlfvenRegion:
                # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
                grouping = [
                    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11, 12, 13, 15, 16, 18, 19, 46, 47, 48, 49, 68, 69, 76,
                     79, 82, 83, 85, 86, 88, 89, 94, 95, 98, 99, 100, 101, 102, 103,104, 106, 108,109, 110, 111, 112, 113, 117,
                    119, 120, 121, 122,125, 128, 129, 130]
                ]

                # investigate other components
                # grouping += [[i] for i in range(10,18)]

                # the "noise" data
                limit = 130
                grouping += [[i for i in range(limit, 3 * SSA_window_Size)]]

                # get the good stuff
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
        elif SSA_window_Size == 1001:
            # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
            grouping = [
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 18, 19, 46, 47, 48, 49, 68, 69, 76,
                 79, 82, 83, 85, 86, 88, 89, 94, 95, 98, 99, 100, 101, 102, 103, 104, 106, 108, 109, 110, 111, 112, 113,
                 117,
                 119, 120, 121, 122, 125, 128, 129, 130]
            ]

            # investigate other components
            # grouping += [[i] for i in range(10,18)]

            # the "noise" data
            limit = 130
            grouping += [[i for i in range(limit, 3 * SSA_window_Size)]]

            # get the good stuff
            goodStuff = []
            for i in range(limit):  # should be 130
                if i not in grouping[0]:
                    goodStuff.append(i)

            grouping += [goodStuff]
            # # identify the harmonic components (maybe 32,33 harmonic or no)
            grouping = [
                [i for i in range(26)] +
                [30, 31, 32, 33, 36, 46,
                 47, 52, 53, 61, 62, 63,
                 86, 87, 92, 93, 102, 103,
                 108, 109, 118, 119, 121, 122,
                 129, 130, 131, 132, 134, 135, 136,
                 137, 139, 140, 148, 149, 160, 161, 162,
                 163, 166, 167, 169, 170, 171, 172, 173,
                 176, 177, 183, 186, 187, 189, 190,
                 193, 194, 199, 200, 201, 203, 204,
                 209, 211, 212, 214, 215, 218, 219,
                 223, 224, 226, 227, 228, 229,
                 231, 232, 236, 237, 238, 239,
                 243, 244, 245, 246, 247, 248, 249,
                 250, 251, 253, 254, 256, 257
                 ]
            ]

            # investigate other components
            grouping += [[i] for i in range(260, 271)]

            # the "noise" data
            limit = 270
            grouping += [[i for i in range(limit, 3 * SSA_window_Size)]]

            # get the good stuff
            goodStuff = []
            for i in range(limit):  # should be 130
                if i not in grouping[0]:
                    goodStuff.append(i)

            # goodStuff = [34, 35, 37, 38, 39, 40]

            grouping += [goodStuff]

        else:
            grouping = [
                [i] for i in range(20)
            ]

    # --- LF Groupings ---
    elif wRocket == 5:

        if SSA_window_Size == 501:

            if useAlfvenRegion:

                # [dt.datetime(2022, 11, 20, 17, 24, 25, 000000), dt.datetime(2022, 11, 20, 17, 25, 18, 000000)]  # ACTUAL WINDOW for WL 501 Low Flyer
                # identify the harmonic components
                grouping = [
                    [0, 1, 4, 5, 15, 16]
                ]

                # investigate other components
                # grouping += [[i] for i in range(70,81)]

                # show the "noise"
                grouping += [[i for i in range(150,3*SSA_window_Size)]]

                # get the good stuff
                goodStuff = []
                for i in range(150):  # should be 150
                    if i not in grouping[0]:
                        goodStuff.append(i)

                grouping += [goodStuff]
            else:
                grouping = [
                    [i] for i in range(20)
                ]
        elif SSA_window_Size == 1001:
            grouping = [
                [0,1,2,3,4,5, 6, 7, 8, 9, 14, 15, 18, 19,20, 26, 29, 30, 31]
            ]

            # investigate other components
            # grouping += [[i] for i in range(36,45)]

            # show the "noise"
            grouping += [[i for i in range(200, 3 * SSA_window_Size)]]

            # get the good stuff
            goodStuff = []
            for i in range(200):  # should be 150
                if i not in grouping[0]:
                    goodStuff.append(i)

            grouping += [goodStuff]

        else:
            grouping = [
                [i] for i in range(20)
            ]




    return group + grouping