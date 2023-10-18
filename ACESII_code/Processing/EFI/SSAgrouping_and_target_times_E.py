def timeWindow(wTargetTimes,wRocket): # MAKE SURE THESE ARE THE SAME AS THOSE IN SSAgrouping_and_Target_times_B

    from ACESII_code.Processing.Magnetometer.SSAgrouping_and_target_times_B import timeWindow

    return timeWindow(wTargetTimes,wRocket)

def groupings(wRocket,SSA_window_Size,wUseData):

    # the first component is the original data
    group = [[i for i in range(SSA_window_Size*3)]]

    # --- High Flyer Groupings ---
    if wRocket == 4:
        print('There aint no HF E-Field, fool!')
    # --- LF Groupings ---
    elif wRocket == 5:

        if SSA_window_Size == 501:

            if wUseData ==0: # THE ALFVEN REGION

                # [dt.datetime(2022, 11, 20, 17, 24, 25, 000000), dt.datetime(2022, 11, 20, 17, 25, 18, 000000)]  # ACTUAL WINDOW for WL 501 Low Flyer
                # identify the harmonic components
                grouping = [
                    [0,1,2,3,4,5,6,7]
                ]

                # investigate other components
                grouping += [[i] for i in range(30,40)]

                # show the "noise"
                limit = 1000
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
                    [i] for i in range(10)
                ]

        else:
            grouping = [
                [i] for i in range(10)
            ]




    return group + grouping