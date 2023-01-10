from files_p import mag_file,root,user_path
from functions import write_var_to_file,write_var_to_file_mag
from cdflib import cdfwrite
from Variables_p import mag_data as mag,epoch_mag,mag_info,zvars_mag
import numpy as np

print('Processing mag data: ')

from files_p import mag_file,root,user_path
from functions import write_var_to_file,write_var_to_file_mag
from cdflib import cdfwrite
from Variables_p import mag_data as mag,epoch_mag,mag_info,zvars_mag
import numpy as np, cdflib,matplotlib.pyplot as plt


X_mag = np.zeros(shape=(len(mag)),dtype='float64')
Y_mag =np.zeros(shape=(len(mag)),dtype='float64')
Z_mag = np.zeros(shape=(len(mag)),dtype='float64')
corrected_epoch = np.zeros(shape=(len(epoch_mag)),dtype='float64')

fillval_mag = 4294967295.0

for i in range(len(mag)):
    X_mag[i] =  mag[i][0]
    Y_mag[i] =  mag[i][1]
    Z_mag[i]=  mag[i][2]


#Correct the corrupted Epoch Values
for i in range(len(epoch_mag)-1):
    checker = (epoch_mag[i+1] - epoch_mag[i])
    mag_checker = [X_mag[i],Y_mag[i],Z_mag[i]]

    if i < (len(epoch_mag) /2):
        if mag_checker[0] > 48000 or mag_checker[0] < 12000:
            X_mag[i] = fillval_mag
        elif mag_checker[1] > 50000 or mag_checker[1] < 18000:
            Y_mag[i] = fillval_mag
        elif mag_checker[2] > 27000 or mag_checker[2] < 10000:
            Z_mag[i] = fillval_mag
    elif i >= (len(epoch_mag) - 12000):
        if mag_checker[0] > 35000 or mag_checker[0] < 30000:
            X_mag[i] = fillval_mag
        elif mag_checker[1] > 35000 or mag_checker[1] < 30000:
            Y_mag[i] = fillval_mag
        elif mag_checker[2] > 17500 or mag_checker[2] < 10000:
            Z_mag[i] = fillval_mag

    if checker <0:
        corrected_epoch[i] = 0
        X_mag[i] = fillval_mag
        Y_mag[i] = fillval_mag
        Z_mag[i] = fillval_mag

    elif checker > 1000000:
        corrected_epoch[i+1] = 0
        X_mag[i+1] = fillval_mag
        Y_mag[i+1] = fillval_mag
        Z_mag[i+1] = fillval_mag

    elif (checker > 0 and checker < 1000000) or (checker > 1000000):
        corrected_epoch[i + 1] = 0
        X_mag[i + 1] = fillval_mag
        Y_mag[i + 1] = fillval_mag
        Z_mag[i + 1] = fillval_mag
    else:
        corrected_epoch[i] = epoch_mag[i]



#corect epoch
bad_indicies = np.where(corrected_epoch==0)
corrected_epoch = np.delete(corrected_epoch,bad_indicies)
X_mag = np.delete(X_mag,bad_indicies)
Y_mag = np.delete(Y_mag,bad_indicies)
Z_mag = np.delete(Z_mag,bad_indicies)

#Correct the  corrupted Mag Values
indicies1 = np.where(X_mag == fillval_mag)
indicies2 = np.where(Y_mag == fillval_mag)
indicies3 = np.where(Z_mag == fillval_mag)

index1 = np.append(indicies1,indicies2)
indexes = np.append(index1,indicies3)
bad_indicies = np.unique(indexes)

X_mag = np.delete(X_mag,bad_indicies)
Y_mag = np.delete(Y_mag,bad_indicies)
Z_mag = np.delete(Z_mag,bad_indicies)
corrected_epoch = np.delete(corrected_epoch,bad_indicies)

print(len(corrected_epoch),len(X_mag),len(Y_mag),len(Z_mag))



plt.subplot(3,1,1)
plt.title('X_mag')
plt.plot(corrected_epoch,X_mag)

plt.subplot(3,1,2)
plt.title('Y_mag')
plt.plot(corrected_epoch,Y_mag)

plt.subplot(3,1,3)
plt.title('Z_mag')
plt.plot(corrected_epoch,Z_mag)
plt.show()


Mag_file_output = cdfwrite.CDF(user_path + root + '/output/' + 'CAPERII_magnetometer_data',cdf_spec=mag_info,delete=True)

vardata = X_mag
attributes = ['Mag_X', 'nT', 'linear', X_mag.min(), X_mag.max(),mag_file.varattsget(zvars_mag[3], expand=True)]
varinfo = mag_file.varinq(zvars_mag[3])
write_var_to_file_mag(Mag_file_output, varinfo, vardata, attributes)


vardata = Y_mag
attributes = ['Mag_Y', 'nT', 'linear', Y_mag.min(), Y_mag.max(),mag_file.varattsget(zvars_mag[3], expand=True)]
varinfo = mag_file.varinq(zvars_mag[3])
write_var_to_file_mag(Mag_file_output, varinfo, vardata, attributes)


vardata = Z_mag
attributes = ['Mag_Z', 'nT', 'linear', Z_mag.min(), Z_mag.max(),mag_file.varattsget(zvars_mag[3], expand=True)]
varinfo = mag_file.varinq(zvars_mag[3])
write_var_to_file_mag(Mag_file_output, varinfo, vardata, attributes)

vardata = corrected_epoch
attributes = ['Epoch', 'ns', 'linear', epoch_mag.min(), epoch_mag.max(),mag_file.varattsget(zvars_mag[0], expand=True)]
varinfo = mag_file.varinq(zvars_mag[0])
write_var_to_file(Mag_file_output, varinfo, vardata, attributes)





print('Done',end='')
