
def groupings(wRocket,SSA_window_Size,useAlfvenRegion):

    # the first component is the original data
    group = [[i for i in range(SSA_window_Size*3)]]

    # --- High Flyer Groupings ---
    if wRocket == 4:
        print('There aint no HF E-Field, fool!')
    # --- LF Groupings ---
    elif wRocket == 5:

        if SSA_window_Size == 501:

            if useAlfvenRegion:

                # [dt.datetime(2022, 11, 20, 17, 24, 25, 000000), dt.datetime(2022, 11, 20, 17, 25, 18, 000000)]  # ACTUAL WINDOW for WL 501 Low Flyer
                # identify the harmonic components
                grouping = [
                    [0,1, 2]
                ]

                # investigate other components
                grouping += [[i] for i in range(0,4)]

                # show the "noise"
                limit = 200
                grouping += [[i for i in range(limit,3*SSA_window_Size)]]

                # # get the good stuff
                goodStuff = []
                for i in range(limit):  # should be 150
                    if i not in grouping[0]:
                        goodStuff.append(i)
                #
                grouping += [goodStuff]
            else:
                grouping = [
                    [i] for i in range(20)
                ]

        elif SSA_window_Size == 101:
            grouping = [
                [0,2,5]
            ]

            # investigate other components
            grouping += [[i] for i in range(0,10)]

            # show the "noise"
            limit = 50
            grouping += [[i for i in range(limit, 3 * SSA_window_Size)]]

            # get the good stuff
            goodStuff = []
            for i in range(limit):  # should be 150
                if i not in grouping[0]:
                    goodStuff.append(i)

            grouping += [goodStuff]

        else:
            grouping = [
                [i] for i in range(20)
            ]




    return group + grouping