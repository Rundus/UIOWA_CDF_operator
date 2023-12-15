
def groupings(wRocket, SSA_window_Size, subsetNo):

    # the first component is the original data
    group = [[i for i in range(SSA_window_Size*3)]]

    # --- High Flyer Groupings ---
    if wRocket == 4:

        if subsetNo == 0: # THE ALFVEN RGION
            # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
            grouping = [
                [0,1,2,3,4,5,6,7]
            ]

            # investigate other components
            grouping += [[i] for i in range(7,10)]

            # the "noise" data
            limit = 3 * SSA_window_Size-1
            grouping += [[i for i in range(limit, 3 * SSA_window_Size)]]

            # get the good stuff
            goodStuff = []
            for i in range(limit): # should be 130
                if i not in grouping[0]:
                    goodStuff.append(i)

            grouping += [goodStuff]

        elif subsetNo == 1:
            # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
            grouping = [
                []
            ]

            # investigate other components
            grouping += [[i] for i in range(10)]
        elif subsetNo == 2:
            # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
            grouping = [
                []
            ]

            # investigate other components
            grouping += [[i] for i in range(10)]
        elif subsetNo == 3:
            # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
            grouping = [
                []
            ]

            # investigate other components
            grouping += [[i] for i in range(10)]
        elif subsetNo == 4:
            # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
            grouping = [
                []
            ]

            # investigate other components
            grouping += [[i] for i in range(10)]


    # --- LF Groupings ---
    elif wRocket == 5:

        if SSA_window_Size == 501:
            if subsetNo == 0:  # THE ALFVEN RGION
                # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
                grouping = [
                    []
                ]

                # investigate other components
                grouping += [[i] for i in range(10)]

                # the "noise" data
                # limit = 250
                # grouping += [[i for i in range(limit, 3 * SSA_window_Size)]]

                # # get the good stuff
                # goodStuff = []
                # for i in range(limit): # should be 130
                #     if i not in grouping[0]:
                #         goodStuff.append(i)
                #
                # grouping += [goodStuff]

            elif subsetNo == 1:
                # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
                grouping = [
                    []
                ]

                # investigate other components
                grouping += [[i] for i in range(10)]
            elif subsetNo == 2:
                # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
                grouping = [
                    []
                ]

                # investigate other components
                grouping += [[i] for i in range(10)]
            elif subsetNo == 3:
                # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
                grouping = [
                    []
                ]

                # investigate other components
                grouping += [[i] for i in range(10)]
            elif subsetNo == 4:
                # # identify the harmonic components (maybe 10,11,16,17 harmonic or no)
                grouping = [
                    []
                ]

                # investigate other components
                grouping += [[i] for i in range(10)]



    return group + grouping