# --- tmCDF_to_L0.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Read in CDF telemetry files and processes them to L0


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools
import os.path
import time
from ACESII_code.class_var_func import Done, prgMsg,setupPYCDF
start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files in
justPrintFileNames = False

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# select which files to convert
# [] --> all files
# [#1,#2,...etc] --> only specific files
wFiles = [3]

getMAGdata = True # get the mag data and the ESA data
justgetMAGdata = True # Only get the MAG data


#########################
# --- Special Toggles ---
#########################
plugHolesInData = False

# select how much of the back-end of the data not to include, given in decimal percents
wPercent = [0.04, 0.08]

# To fix the gSync Epoch dropouts, I get the previously non-fillval epoch values and add the median delta T to it
deltaT_median = 249984 # nanoseconds, in TT2000

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts,TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder,fliers
from ACESII_code.class_var_func import color, L0_TRICE_Quick
from glob import glob
from os.path import getsize
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)

def tmCDF_to_L0(wRocket, wFile, rocketFolderPath, justPrintFileNames,wflyer):

    # --- --- --- --- --- --- ---
    # --- FILE SPECIFIC DATA ---
    # --- --- --- --- --- --- ---

    # --- ACES II Flight/Integration Data ---
    if wRocket in [0,1,4,5]:
        rocketAttrs, data_dicts, deConvolveKeys = ACES_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
        globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L0'
        L0ModelData = L0_TRICE_Quick(wflyer)

    # --- TRICE II ---
    elif wRocket in [2,3]:
        globalAttrsMod = {}
        rocketAttrs = TRICE_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        L0ModelData = L0_TRICE_Quick(wflyer)

    # Set the paths for the file names
    tmCDFFiles = glob(f'{rocketFolderPath}tmCDF\{fliers[wflyer]}\*.cdf')
    L0Files = glob(f'{rocketFolderPath}L0\{fliers[wflyer]}\*.cdf')
    tmCDF_names = [ifile.replace(f'{rocketFolderPath}tmCDF\{fliers[wflyer]}\\', '') for ifile in tmCDFFiles]
    tmCDF_names_searchable = [cdfname.replace('.cdf','').replace('_','').replace('36359','').replace('36364','') for cdfname in tmCDF_names]
    L0_names = [ofile.replace(f'{rocketFolderPath}L0\{fliers[wflyer]}\\', '').replace(".cdf", '') for ofile in L0Files]
    L0_names_searchable = [ fname.replace('ACES_','').replace('l0_','').replace('eepaa_','').replace('leesa_','').replace('iepaa_','').replace('_v00','').replace('_','').replace('eepaa','').replace('leesa','').replace('iepaa','').replace('36359','').replace('36364','') for fname in L0_names]
    dataFile_name = tmCDFFiles[wFile].replace(f'{rocketFolderPath}tmCDF\{fliers[wflyer]}\\', '')

    # Output Naming Format: EEPAA, LEESA, IEPAA
    if wRocket in [0, 1]: # Integration Files
        dF_name = dataFile_name.replace('36_359_', '').replace('36_364_', '').replace('.cdf', '')
        fileoutName = [f'{rocketAttrs.missionNam}_{rocketID}_l0_{rocketAttrs.InstrNames_LC[i]}_{dF_name}_v00.cdf' for i in (range(len(rocketAttrs.InstrNames_LC)))]
    elif wRocket in [4,5]: #ACESII flight files
        dF_name = dataFile_name.replace('36359_', '').replace('36364_', '').replace('ACESII_','').replace('flight_',"").replace('_v00','').replace('.cdf', '')
        fileoutName = [f'{rocketAttrs.missionNam}_{rocketID}_l0_{rocketAttrs.InstrNames_LC[i]}_{dF_name}_v00.cdf' for i in (range(len(rocketAttrs.InstrNames_LC)))]
    else: #TRICE
        dF_name = dataFile_name.replace(rocketID, '').replace(rocketAttrs.missionNam, '').replace('_', '').replace('k0','').replace('.cdf', '')
        fileoutName = [f'{rocketAttrs.missionNam}_{rocketID}_l0_{rocketAttrs.InstrNames_LC[i]}_{dF_name}_v00.cdf' for i in (range(len(rocketAttrs.InstrNames_LC)))]

    # --- --- --- --- --- --- --- ------
    # --- PRODUCE TELEM DATA FILES ---
    # --- --- --- --- --- --- --- ------
    if justPrintFileNames:
        for i, file in enumerate(tmCDFFiles):
            anws = ["yes" if tmCDF_names_searchable[i].replace('.cdf', "") in L0_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made CDF: {:3s} '.format(i, tmCDF_names[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to L0 data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(tmCDFFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg('Getting data from tmCDFfile')
        with pycdf.CDF(tmCDFFiles[wFile]) as tmCDFDataFile:
            epoch_from_file = tmCDFDataFile['Epoch'][...]
            atmCDF_epoch = [pycdf.lib.datetime_to_tt2000(epoch_from_file[i]) for i in range(len(epoch_from_file))]
            atmCDF_sfid = tmCDFDataFile['sfid'][...]
            atmCDF_mf = tmCDFDataFile['minorframe'][...]
        Done(start_time)

        if plugHolesInData:  # if we decide to keep more the mf data before the interesting part
            # --- --- --- --- --- --- ---
            # --- REDUCE INPUT DATA ---
            # --- --- --- --- --- --- ---

            prgMsg('Resizing Data')

            # Don't use the last % of data in order to speed up the code
            if wPercent[wflyer] > 0:
                dontUse_amount = int(len(atmCDF_epoch)*wPercent[wflyer])
                atmCDF_epoch = atmCDF_epoch[0:-1*dontUse_amount]
                atmCDF_sfid = atmCDF_sfid[0:-1*dontUse_amount]
                atmCDF_mf = atmCDF_mf[0:-1*dontUse_amount]

            # Strip away the ends of the data so it starts/ends with a complete majorframe
            mf_startpoint, mf_endpoint = np.where(atmCDF_sfid == 0), np.where(atmCDF_sfid == 39)
            mf_start, mf_end = mf_startpoint[0][0], mf_endpoint[0][-1] + 1 # pythonic + 1

            tmCDF_epoch = np.array([atmCDF_epoch[i] for i in (range(mf_start, mf_end))])
            tmCDF_sfid = atmCDF_sfid[mf_start:mf_end]
            tmCDF_mf = atmCDF_mf[mf_start:mf_end]
            Done(start_time)

            prgMsg('Plugging Holes in Data')
            # --- --- --- --- --- --- --- --- --- ---
            # --- CORRECT GSYNCS AND SFID'S >= 40 ---
            # --- --- --- --- --- --- --- --- --- ---

            if wflyer == 0: # special case to clean up the high flyer data. Remove the single data point where sfid == 40
                index40 = np.where(tmCDF_sfid == 40) # should be at 3764720
                tmCDF_sfid = np.delete(tmCDF_sfid, index40, axis=0)
                tmCDF_mf = np.delete(tmCDF_mf, index40, axis=0)
                tmCDF_epoch = np.delete(tmCDF_epoch, index40, axis=0)


            fix_indicies = []
            counter = tmCDF_sfid[0]
            for index in tqdm(range(len(tmCDF_sfid))):

                if tmCDF_sfid[index] != counter:
                    fix_indicies.append([index, tmCDF_sfid[index], counter, tmCDF_epoch[index]])

                    if tmCDF_sfid[index] == -1:  # fabricate the gSync'd Epochs. No worry about starting with a -1 b/c of the code above
                        tmCDF_epoch[index] = tmCDF_epoch[index - 1] + deltaT_median
                        tmCDF_sfid[index] = counter

                counter += 1
                if counter >= 40:
                    counter = 0

            Done(start_time)
        else:
            tmCDF_epoch = atmCDF_epoch
            tmCDF_sfid = atmCDF_sfid
            tmCDF_mf = atmCDF_mf


        # --- --- --- --- --- ----
        # --- COLLECTING DATA ---
        # --- --- --- --- --- ----

        # --- pre-allocate memory for ESA data---
        No_of_esa_measurements = int(len(tmCDF_sfid)/4)
        for i in (range(rocketAttrs.NumOfInstr - 1)):
            data_dicts[i]['Epoch'][0] = np.zeros(shape=(No_of_esa_measurements), dtype='int64')
            data_dicts[i]['sfid'][0] = np.zeros(shape=(No_of_esa_measurements), dtype='uint8')
            data_dicts[i]['Sector_Counts'][0] = np.zeros(shape=(No_of_esa_measurements, rocketAttrs.num_of_sector_counts[i]), dtype='int32')

        ######################
        # --- GET MAG DATA ---
        ######################

        if getMAGdata:
            prgMsg('Collecting MAG data')
            RingCore_data = np.array([tmCDF_mf[iv][120 - 1] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) % 2 == 0])
            RingCore_sfid = np.array([tmCDF_sfid[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) % 2 == 0])
            RingCore_epoch = np.array([tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) % 2 == 0])

            # package all the MAG data into a data_dict
            data_dict_mag = {
                'RingCore_Data': [RingCore_data, {'DEPEND_0': 'RingCore_epoch', 'DEPEND_1': None, 'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','SCALETYP':'linear'}],
                'RingCore_sfid': [RingCore_sfid, {'DEPEND_0': 'RingCore_epoch', 'DEPEND_1': None, 'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','SCALETYP':'linear'}],
                'RingCore_epoch': [RingCore_epoch, {'DEPEND_0': 'RingCore_epoch', 'DEPEND_1': None, 'DEPEND_2': None,'FILLVAL': rocketAttrs.epoch_fillVal,'FORMAT': 'I5','UNITS': 'ns','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','MONOTON':'INCREASE','TIME_BASE':'J2000','TIME_SCALE':'Terrestrial Time','REFERENCE_POSITION':'Rotating Earth Geoid','SCALETYP':'linear'}]
            }

            if wRocket == 5: # Low Flyer has Ring Core and Tesseract
                Tesseract_data = np.array([tmCDF_mf[iv][16 - 1] for iv in (range(len(tmCDF_sfid)))])
                Tesseract_sfid = np.array(tmCDF_sfid)
                Tesseract_epoch = np.array(tmCDF_epoch)

                data_dict_mag = {**data_dict_mag, **{'Tesseract_data':[Tesseract_data, {'DEPEND_0': 'Tesseract_epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','SCALETYP':'linear'}]}}
                data_dict_mag = {**data_dict_mag, **{'Tesseract_sfid': [Tesseract_sfid, {'DEPEND_0': 'Tesseract_epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','SCALETYP':'linear'}]}}
                data_dict_mag = {**data_dict_mag, **{'Tesseract_epoch': [Tesseract_epoch, {'DEPEND_0': 'Tesseract_epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': rocketAttrs.epoch_fillVal,'FORMAT': 'I5','UNITS': 'ns','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','MONOTON':'INCREASE','TIME_BASE':'J2000','TIME_SCALE':'Terrestrial Time','REFERENCE_POSITION':'Rotating Earth Geoid','SCALETYP':'linear'}]}}

            Done(start_time)

        if not justgetMAGdata:
            # --- variables with special conditions ---
            prgMsg('Collecting variables with special conditions')
            for j in tqdm(range(rocketAttrs.NumOfInstr)):
                if j == 0:
                    data_dicts[j]['Boom_Monitor'][0] = [tmCDF_mf[iv][148 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 1]
                    data_dicts[j]['Epoch_monitors'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 17]
                    data_dicts[j]['28V_Monitor'][0] = [tmCDF_mf[iv][141 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 17]
                    data_dicts[j]['EXP_Current'][0] = [tmCDF_mf[iv][121 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1)%4 == 1]
                elif j == 1:
                    data_dicts[j]['Boom_Monitor'][0] = [tmCDF_mf[iv][148 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 6]
                    data_dicts[j]['Epoch_monitors'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 18]
                    data_dicts[j]['28V_Monitor'][0] = [tmCDF_mf[iv][141 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 18]
                    data_dicts[j]['EXP_Current'][0] = [tmCDF_mf[iv][121 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) % 4 == 0]
                elif j == 2:
                    data_dicts[j]['Boom_Monitor'][0] = [tmCDF_mf[iv][148 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 2]
                    data_dicts[j]['Epoch_monitors'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 33]
                    data_dicts[j]['28V_Monitor'][0] = [tmCDF_mf[iv][141 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 33]
                    data_dicts[j]['EXP_Current'][0] = [tmCDF_mf[iv][121 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) % 4 == 2]
                elif j == 3:

                    data_dicts[j]['Boom_Monitor_1'][0] = [tmCDF_mf[iv][148 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 3]
                    data_dicts[j]['Epoch_monitor_1'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 3]

                    data_dicts[j]['Boom_Monitor_2'][0] = [tmCDF_mf[iv][148 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 4]
                    data_dicts[j]['Epoch_monitor_2'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 4]

                    data_dicts[j]['28V_Monitor'][0] = [tmCDF_mf[iv][141 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 19]
                    data_dicts[j]['Epoch_monitors'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 19]

                    data_dicts[j]['EXP_Current'][0] = [tmCDF_mf[iv][121 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) % 4 == 3]

            Done(start_time)

            # --- Major Frame counter ---
            prgMsg('Counting major frames')
            for ii in (range(len(tmCDF_sfid))): # NOTE: Major frame counter 1 is the first 16 binary digits of the majorframe number and major frame counter 2 is the next 16 binary digits, which totals a 32 binary word.
                if (tmCDF_sfid[ii] + 1) not in [4, 12, 20, 28, 36]:
                    for j in (range(rocketAttrs.NumOfInstr)):  # MajF counter 1
                        data_dicts[j]['major_frame_counter'][0].append([tmCDF_mf[ii][72 - 1]])
                else:
                    for j in (range(rocketAttrs.NumOfInstr)):  # MajF counter 2
                        data_dicts[j]['major_frame_counter'][0].append([tmCDF_mf[ii][136 - 1], tmCDF_mf[ii][72 - 1]])
            Done(start_time)

            # --- LP Data ---
            prgMsg('Collecting LP Data')
            data_dicts[3]['Channel_Counts'][0] = [[tmCDF_mf[iv][word - 1] for word in rocketAttrs.LP_words] for iv in (range(len(tmCDF_sfid))) ]
            data_dicts[3]['Epoch'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid)))]
            data_dicts[3]['sfid'][0] = [tmCDF_sfid[iv] for iv in (range(len(tmCDF_sfid)))]
            Done(start_time)

            # --- --- --- --- --- --- --- --- --- --- --- ---
            # --- COLLECT MINORFRAME DATA FOR THE ESAs  ---
            # --- --- --- --- --- --- --- --- --- --- --- ---

            for iii in range(rocketAttrs.NumOfInstr - 1):
                words = rocketAttrs.ESA_words[iii]
                deConvolveKey = deConvolveKeys[iii]
                deConvolveItems = deConvolveKey.items()
                SectorDataCounter = 0

                for iv in tqdm(range(No_of_esa_measurements)):
                    esaMeasurement = tmCDF_mf[4 * iv: 4 * (iv + 1)]
                    data_dicts[iii]['sfid'][0][iv] = tmCDF_sfid[4 * iv]
                    data_dicts[iii]['Epoch'][0][iv] = tmCDF_epoch[4 * iv]

                    for index in range(len(esaMeasurement)):
                        mfdata = [esaMeasurement[index][k - 1] for k in words]
                        keyset = [val for key, val in deConvolveItems if key in [i for i in range(1 + index * len(words), 1 + (1 + index) * len(words))]]

                        for v in (range(len(mfdata))):  # place all the data
                            if keyset[v] == 'Sector_Counts':
                                data_dicts[iii]['Sector_Counts'][0][iv][SectorDataCounter] = mfdata[v]
                                SectorDataCounter += 1

                                if SectorDataCounter == rocketAttrs.num_of_sector_counts[iii]:
                                    SectorDataCounter = 0

                            elif keyset[v] == 'status_word_2':
                                data_dicts[iii]['HV_div16'][0].append((mfdata[v] >> 14) & 1)
                                data_dicts[iii]['HV_enable'][0].append((mfdata[v] >> 7) & 1)
                                data_dicts[iii]['sweep_step'][0].append((mfdata[v] >> 8) & 63)
                                data_dicts[iii]['TP5_enable'][0].append((mfdata[v] >> 5) & 1)

                            elif keyset[v] == '625kHz_Clock_Input':
                                data_dicts[iii][keyset[v]][0].append(mfdata[v])
                                if mfdata[v] / rocketAttrs.ESA_CLK_INPUT >= 2:
                                    data_dicts[iii]['Count_Interval'][0].append(65535)
                                else:
                                    data_dicts[iii]['Count_Interval'][0].append(mfdata[v] / rocketAttrs.ESA_CLK_INPUT)
                            elif keyset[v] != None:
                                data_dicts[iii][keyset[v]][0].append(mfdata[v])

            # --- --- --- --- ---
            # --- PROCESSING ---
            # --- --- --- --- ---
            prgMsg('Processing Data')
            for jjj in range(rocketAttrs.NumOfInstr):

                # --- Major Frame Counter ---
                MF_counter, get_first_major_frame = [], []
                MFDat = data_dicts[jjj]['major_frame_counter'][0]
                topNum, botNum = 0, 0

                for data in MFDat:
                    if len(data) == 1:
                        botNum = data[0]
                    elif len(data) == 2:
                        botNum = data[1]
                        topNum = data[0]

                        if get_first_major_frame == []: # used to start the minor frame counter
                            if topNum != 0:
                                get_first_major_frame.append( (topNum << 16) - 1 ) # -1 to remove the MF that I start on
                            else:
                                get_first_major_frame.append(botNum - 1)

                    num = (topNum << 16) + botNum
                    MF_counter.append(num)

                data_dicts[jjj]['major_frame_counter'][0] = MF_counter

                # --- minor Frame Counter ---
                data_dicts[jjj]['minor_frame_counter'][0] = np.array([i + (get_first_major_frame[0]*40) for i in (range(1,len(tmCDF_sfid)+1))])

                # --- Convert all data to numpy arrays ---
                for key, dataVal in data_dicts[jjj].items():
                    data_dicts[jjj][key][0] = np.array(dataVal[0])


            Done(start_time)

            # --- --- --- --- --- --- ---
            # --- WRITE OUT THE DATA ---
            # --- --- --- --- --- --- ---
            writeOut_path = [f'{rocketFolderPath}L0\{fliers[wflyer]}\\{fileoutName[i]}' for i in range(rocketAttrs.NumOfInstr)]

            # --- loop through instruments ---
            for v in range(rocketAttrs.NumOfInstr):
                prgMsg(f'Writing out {rocketAttrs.InstrNames[v]} data')
                outputPath = writeOut_path[v]

                # --- delete output file if it already exists ---
                if os.path.exists(outputPath):
                    os.remove(outputPath)

                # --- open the output file ---
                with pycdf.CDF(outputPath,'') as L0File:
                    L0File.readonly(False)

                    # --- write out global attributes ---
                    inputGlobDic = L0ModelData.cdfFile.globalattsget()
                    for key, val in inputGlobDic.items():
                        if key == 'Descriptor':
                            globalAttrsMod[key] = rocketAttrs.InstrNames_Full[v]

                        if key in globalAttrsMod:
                            L0File.attrs[key] = globalAttrsMod[key]
                        else:
                            L0File.attrs[key] = val

                    # --- WRITE OUT DATA ---
                    for varKey, varVal in data_dicts[v].items():
                        if varKey == 'Epoch' or varKey == 'Epoch_monitors':
                            L0File.new(varKey, data=varVal[0], type = 33)
                        elif varKey == 'Sector_Counts':
                            L0File.new(varKey, data=varVal[0], type=pycdf.const.CDF_UINT2)
                        else:
                            L0File.new(varKey, data=varVal[0])

                        # --- Write out the attributes and variable info ---
                        for attrKey, attrVal in data_dicts[v][varKey][1].items():
                            if attrKey == 'VALIDMIN':
                                L0File[varKey].attrs[attrKey] = varVal[0].min()
                            elif attrKey == 'VALIDMAX':
                                L0File[varKey].attrs[attrKey] = varVal[0].max()
                            elif attrVal != None:
                                L0File[varKey].attrs[attrKey] = attrVal

                Done(start_time)



        if getMAGdata:
            prgMsg('Writing out MAG data')

            if wRocket == 4:
                fileoutName = fileoutName[0].replace('eepaa_','mag_').replace('iepaa_','mag_').replace('leesa_','mag_')
            elif wRocket == 5:
                fileoutName = fileoutName[0].replace('eepaa_', 'mag&tesseract_').replace('iepaa_', 'mag&tesseract_').replace('leesa_','mag&tesseract_')

            outputPath = f'{rocketFolderPath}mag\{fliers[wflyer]}\\{fileoutName}'


            # --- delete output file if it already exists ---
            if os.path.exists(outputPath):
                os.remove(outputPath)

            # --- open the output file ---
            with pycdf.CDF(outputPath, '') as magFile:
                magFile.readonly(False)

                # --- write out global attributes ---
                inputGlobDic = L0ModelData.cdfFile.globalattsget()
                for key, val in inputGlobDic.items():
                    if key in globalAttrsMod:
                        magFile.attrs[key] = globalAttrsMod[key]
                    else:
                        magFile.attrs[key] = val

                # --- WRITE OUT DATA ---
                for varKey, varVal in data_dict_mag.items():
                    if 'epoch' in varKey.lower():  # epoch data
                        magFile.new(varKey, data=varVal[0], type=33)
                    else:  # other data
                        magFile.new(varKey, data=varVal[0], type=pycdf.const.CDF_REAL8)

                    # --- Write out the attributes and variable info ---
                    for attrKey, attrVal in data_dict_mag[varKey][1].items():
                        if attrKey == 'VALIDMIN':
                            magFile[varKey].attrs[attrKey] = varVal[0].min()
                        elif attrKey == 'VALIDMAX':
                            magFile[varKey].attrs[attrKey] = varVal[0].max()
                        elif attrVal != None:
                            magFile[varKey].attrs[attrKey] = attrVal

            Done(start_time)




# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 0:  # ACES II Integration High
    rocketFolderPath = Integration_data_folder
    wflyer = 0
elif wRocket == 1: # ACES II Integration Low
    rocketFolderPath = Integration_data_folder
    wflyer = 1
elif wRocket == 2:  # TRICE II High
    rocketFolderPath = TRICE_data_folder
    wflyer = 0
elif wRocket == 3: # TRICE II Low
    rocketFolderPath = TRICE_data_folder
    wflyer = 1
elif wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}tmCDF\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        tmCDF_to_L0(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}tmCDF\{fliers[wflyer]}\*.cdf')))):
            tmCDF_to_L0(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            tmCDF_to_L0(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)