# import matplotlib.pyplot as plt
#
# from ACESII_code.myImports import *
# from numpy.fft import rfft, fftfreq
#
#
# inputFile_deltaB = glob('C:\Data\ACESII\L3\deltaB\low\*Field_Aligned*')[0] # get the deltaB data
# inputFile_deltaE = glob('C:\Data\ACESII\L3\deltaE\low\*Field_Aligned*')[0] # get the deltaE data
# inputFile_B = 'C:\Data\ACESII\L2\low\ACESII_36364_l2_RingCore_Field_Aligned.cdf' # get the deltaB data
# inputFile_E = 'C:\Data\ACESII\L2\low\ACESII_36364_l2_E_Field_Field_Aligned_downsampled.cdf' # get the deltaE data
# inputFile_Langmuir = 'C:\Data\ACESII\L3\Langmuir\low\ACESII_36364_langmuir_fixed.cdf'
#
# ###############
# ### TOGGLES ###
# ###############
# wRegions = [0,1,2]
# wKeys = [0,1] # 0=> B_East, E_North, 1=> B_North,-E_East
# iLatRegions = [[71.25, 71.353205], # Quiet Time
#                [71.582, 71.68], # Dispersed Time
#                [71.86, 71.94]] # spatially aligned
# regionName = ['Quiet','Time_Aligned','Spatially_Aligned']
# freqLim = [0.01,15]
#
# testDat = []
#
# for wRegion in wRegions:
#     print(wRegion)
#     targetVar = [iLatRegions[wRegion], 'ILat']  # primary alfven region
#
#     data_dict_deltaB = deepcopy(loadDictFromFile(inputFile_deltaB, targetVar=targetVar))
#     data_dict_deltaE = deepcopy(loadDictFromFile(inputFile_deltaE, targetVar=targetVar))
#     data_dict_B = deepcopy(loadDictFromFile(inputFile_B,targetVar=targetVar))
#     data_dict_E = deepcopy(loadDictFromFile(inputFile_E,targetVar=targetVar))
#     data_dict_langmuir = deepcopy(loadDictFromFile(inputFile_Langmuir,targetVar=targetVar,wKeys_Load=['ni','Epoch','ILat']))
#
#     # downsample the langmuir data
#     indexVals = [np.abs(data_dict_langmuir['ILat'][0]-ILat).argmin() for ilt, ILat in enumerate(data_dict_B['ILat'][0])]
#     data_dict_langmuir['ni'][0] = deepcopy(data_dict_langmuir['ni'][0][indexVals])
#     data_dict_langmuir['Epoch'][0] = deepcopy(data_dict_langmuir['Epoch'][0][indexVals])
#     data_dict_langmuir['ILat'][0] = deepcopy(data_dict_langmuir['ILat'][0][indexVals])
#
#
#     ############################
#     # normalization temp plotting
#     ############################
#     dB_e = data_dict_deltaB['B_e'][0] # in tesla
#     dB_r = data_dict_deltaB['B_r'][0] # in tesla
#
#     dE_e = data_dict_deltaE['E_e'][0] # in mV/m
#     dE_r = data_dict_deltaE['E_r'][0] # in mV/m
#
#     B_e = 1E-9 * data_dict_B['B_e'][0] # in tesla
#     B_p = 1E-9 * data_dict_B['B_p'][0]
#     B_r = 1E-9 * data_dict_B['B_r'][0]
#     B0 =  1E-9 * data_dict_B['Bmag'][0]
#
#     E_e = data_dict_E['E_e'][0] # in V/m
#     E_p = data_dict_E['E_p'][0]
#     E_r = data_dict_E['E_r'][0]
#
#     ni = (cm_to_m**3) *data_dict_langmuir['ni'][0]
#
#     m_i_avg = 2.45E-26
#     VA_t = B0 /np.sqrt(u0*ni*m_i_avg)
#
#
#     for key in wKeys:
#
#         if key == 0:
#             Eperp_norm = E_r /(np.abs(VA_t)*np.abs(B0))
#             Bperp_norm =B_e / (np.abs(B0))
#             labelB = '$B_{e}(t)/|B_{0}(t)|$'
#             labelE = '$E_{r}(t)/|B_{0}(t)||V_{A}(t)|$'
#             pngTitle='ErBe'
#             deltaB = dB_e
#             deltaE = dE_r
#             labeldeltaB = '$\delta B_{e}(t)$ [nT]'
#             labeldeltaE = '$\delta E_{r}(t)$ [mV/m]'
#         else:
#             Eperp_norm = E_e / (np.abs(VA_t) * np.abs(B0))
#             Bperp_norm = B_r / (np.abs(B0))
#             labelB = '$B_{r}(t)/|B_{0}(t)|$'
#             labelE = '$E_{e}(t)/|B_{0}(t)||V_{A}(t)|$'
#             pngTitle = 'EeBr'
#             deltaB = dB_r
#             deltaE = -1*dE_e
#             labeldeltaB = '$ \delta B_{r}(t)$ [nT]'
#             labeldeltaE = '$-\delta E_{e}(t)$ [mV/m]'
#
#         # take the fourier transforms
#         N, T = len(data_dict_B['Epoch'][0]), 1 / 128
#         xf_B = fftfreq(N, T)[:N // 2]
#         N, T = len(data_dict_E['Epoch'][0]), 1 / 128
#         xf_E = fftfreq(N, T)[:N // 2]
#         yf_E = rfft(Eperp_norm)
#         FFT_E = 2.0 / N * np.abs(yf_E[0:N // 2])
#         yf_B = rfft(Bperp_norm)
#         FFT_B = 2.0 / N * np.abs(yf_B[0:N // 2])
#
#         V_A_ratio = np.abs(FFT_E / FFT_B)
#
#         ####################
#         ##### PLOTTING #####
#         ####################
#         fig, ax = plt.subplots(3)
#         fig.set_size_inches(6, 8)
#         fig.suptitle(f'{regionName[wRegion]} Region')
#
#         ax[0].plot(data_dict_deltaB['ILat'][0], deltaB, color='blue',label=labeldeltaB)
#         ax[0].plot(data_dict_deltaE['ILat'][0], deltaE, color='red', label=labeldeltaE)
#         ax[0].legend(loc='upper right')
#
#         ax[1].plot(xf_B, FFT_B,color='blue', label=labelB)
#         ax[1].plot(xf_E, FFT_E,color='red', label=labelE)
#         ax[1].set_ylim(1E-7, 1E-2)
#         ax[1].set_yscale('log')
#         # ax[0].set_xscale('log')
#         ax[1].set_xlim(freqLim[0],freqLim[1])
#         ax[1].legend(loc='upper right')
#
#         for i in range(5):
#             ax[1].axvline(0.55*(i+1),color='orange',alpha=0.5)
#             ax[2].axvline(0.55 * (i + 1), color='orange', alpha=0.5)
#
#         ax[2].plot(xf_B, V_A_ratio, color='black',linestyle='-',marker='o',markersize=1.3,linewidth=0.5)
#         # ax[1].scatter(xf_B, V_A_ratio, color='black', marker='o', s=40)
#         ax[2].set_ylabel('E\'/B\'')
#         ax[2].set_xlabel('Frequency [Hz]')
#         ax[2].axhline(1,color='red',alpha=0.5)
#         # ax[1].axhline(1.1, color='green', alpha=0.5)
#         # ax[1].axhline(0.9, color='green', alpha=0.5)
#         # ax[1].set_xscale('log')
#         ax[2].set_ylim(0,4)
#         ax[2].set_xlim(freqLim[0],freqLim[1])
#
#         plt.tight_layout()
#
#         plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\TEST\{regionName[wRegion]}_{pngTitle}.png')

import numpy as np
a = np.array([[1,2,3],[1,2,3],[1,2,3]])
b = np.array([[2,2,2],[5,5,5],[1,1,1]])

print(a/b)