import matplotlib.pyplot as plt
from ACESII_code.myImports import *
from numpy.fft import rfft, fftfreq

# inputFiles_elec = glob(f'C:\Data\ACESII\L3\E.cdf')[0]
# inputFiles_mag = glob(f'C:\Data\ACESII\L3\B.cdf')[0]

inputFiles_elec = glob(f'C:\Data\ACESII\L2\low\ACESII_36364_l2_E_Field_ENU_downsampled.cdf')[0]
inputFiles_mag = glob(f'C:\Data\ACESII\L2\low\ACESII_36364_l2_RingCore_ENU.cdf')[0]

data_dict_mag = loadDictFromFile(inputFiles_mag)
data_dict_elec = loadDictFromFile(inputFiles_elec)
Bvec = np.array([data_dict_mag['B_East'][0],data_dict_mag['B_North'][0],data_dict_mag['B_Up'][0]]).T
Evec = np.array([data_dict_elec['E_East'][0],data_dict_elec['E_North'][0],data_dict_elec['E_Up'][0]]).T

# Bvec = np.array([   [data_dict_mag['Mag_East'][0][i], data_dict_mag['Mag_North'][0][i], data_dict_mag['Mag_Up'][0][i]]  for i in range(len(data_dict_mag['Mag_East'][0])) ]).T
# Evec = np.array([   [data_dict_elec['EFI_East'][0][i], data_dict_elec['EFI_North'][0][i], data_dict_elec['EFI_Up'][0][i]]  for i in range(len(data_dict_mag['Mag_East'][0])) ]).T
#
# B =np.array( [ [val[0][0],val[1][0],val[2][0]] for val in Bvec])
# E = np.array([ [val[0][0],val[1][0],val[2][0]] for val in Evec])


S = (1/u0)*np.array([np.cross(Evec[i], Bvec[i]) for i in range(len(Bvec)) ])

fig,ax = plt.subplots()
ax.plot([i for i in range(len(Bvec))],S[:,2])

plt.show()

