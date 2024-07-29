import csv
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

path = r'C:\Users\cfelt\OneDrive\Desktop\Edge of Space\2024\HALCEON (Geiger)\geiger_data.csv'
timestamp = []
count1 = []
count2 = []
count3 = []

loopCounter = 0

with open(path, newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    try:
        for row in spamreader:

            if loopCounter >1:
                timestamp.append(row[0])
                count1.append(int(row[1]))
                count2.append(int(row[2]))
                count3.append(int(row[3]))
            loopCounter +=1
    except:
        print('Done')


# --- PROCESS THE DATA ---
# Clean datetimes
timeStamps_dt = [dt.datetime.strptime(date_string, '%Y-%m-%d %H:%M:%S') for date_string in timestamp]

# find the extrema in count3

for i in range(len(count3)-1):
    nextVal = count3[i+1]
    thisVal = count3[i]

    if np.abs(nextVal-thisVal)>600:
        print(thisVal,nextVal,np.abs(nextVal-thisVal),i)




# --- Plot Everything ---
PlotRawData = False
PlotSummedData = True
if PlotRawData:
    fig, ax = plt.subplots(3)
    ax[0].plot(timeStamps_dt, count1)
    ax[0].set_ylim(0,20)
    ax[1].plot(timeStamps_dt, count2)
    ax[1].set_ylim(0,20)
    ax[2].plot(timeStamps_dt, count3)
    plt.show()
if PlotSummedData:
    n = 30
    count1 = [sum(count1[i:i + n]) for i in range(0, len(count1), n)]
    count2 = [sum(count2[i:i + n])for i in range(0, len(count2), n)]
    count3 = [sum(count3[i:i + n]) for i in range(0, len(count3), n)]
    timestampNum = [i for i in range(len(count1))]

    fig, ax = plt.subplots(3)
    ax[0].set_title('Count1')
    ax[0].plot(timestampNum, count1)
    ax[0].set_ylim(-3, 60)
    ax[1].set_title('Count2')
    ax[1].plot(timestampNum, count2)
    ax[1].set_ylim(-4, 60)
    ax[2].set_title('Count3')
    ax[2].plot(timestampNum, count3)
    ax[2].set_xlabel('Sample #')
    plt.tight_layout()
    plt.show()




