# --- missionAttributes.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to store the attributes specific to certain rockets


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---

import numpy as np
import datetime as dt

# --- --- --- --- ---
# --- ATTRIBUTES ---
# --- --- --- --- ---

class makeRocketAttrs:
    def __init__(self, missionAttrs):
        self.missionNam = missionAttrs['missionNam']
        self.rocketID = missionAttrs['rocketID']
        self.globalAttributes = missionAttrs['globalAttributes']
        self.g_nWords = missionAttrs['g_nWords']
        self.InstrNames = missionAttrs['InstrNames']
        self.InstrNames_LC = missionAttrs['InstrNames_LC']
        self.InstrNames_Full = missionAttrs['InstrNames_Full']
        self.ESA_words = missionAttrs['ESA_words']
        self.nDataLen = missionAttrs['nDataLen']
        self.gSync = missionAttrs['gSync']
        self.ESA_CLK_INPUT = missionAttrs['ESA_CLK_INPUT']
        self.LP_words = missionAttrs['LP_words']
        self.words_per_mf = missionAttrs['words_per_mf']
        self.num_of_sector_counts = missionAttrs['num_of_sector_counts']
        self.g_nNomRecLen = missionAttrs['g_nNomRecLen']
        self.nNomDataLen = missionAttrs['nNomDataLen']
        self.NumOfInstr = missionAttrs['NumOfInstr']
        self.epoch_fillVal = missionAttrs['epoch_fillVal']
        self.mfPerFullSweep = missionAttrs['mfPerFullSweep']
        self.Instr_Energy = missionAttrs['Instr_Energy']
        self.geometric_factor = missionAttrs['geometric_factor']
        self.deadtime = missionAttrs['deadtime']
        self.Instr_sector_to_pitch = missionAttrs['Instr_sector_to_pitch']
        self.Instr_words_per_sweep = missionAttrs['Instr_words_per_sweep']
        self.Instr_Energy = missionAttrs['Instr_Energy']
        self.esaMaxCounts = missionAttrs['esaMaxCounts']
        self.LP_Variables = missionAttrs['LP_Variables']
        self.MinorFrameTime = missionAttrs['MinorFrameTime']
        self.LPSamplePeriod = missionAttrs['LPSamplePeriod']
        self.Launch_Times = missionAttrs['Launch_Times']
        self.Andoya_magnetic_inclination = missionAttrs['Andoya_magnetic_inclination']
        self.LPFixed_calResistances = missionAttrs['LPFixed_calResistances']
        self.LPswept_voltage_range = missionAttrs['LPswept_voltage_range']
        self.LPFixedProbeBias = missionAttrs['LPFixedProbeBias']
        self.Epoch_range_to_determine_stepDAC = missionAttrs['Epoch_range_to_determine_stepDAC']



# --- --- --- --- --- ---
# --- ACES II Mission ---
# --- --- --- --- --- ---
def ACES_mission_dicts():

    from copy import deepcopy

    ACESII_attrs_dict = {
        'Andoya_magnetic_inclination':78.1300,
        'missionNam': 'ACESII',
        'rocketID': ['36359', '36364'],
        'globalAttributes': [
            {'Source_name': f'ACESII_36359>ACES II 36.359',
             'Data_type': 'K0>Key Parameter',
             'PI_name': 'Bounds',
             'Logical_source': f'aces_36359_',
             'Logical_file_id': f'aces_36359_00000000_v01',
             'Logical_source_description': 'Raw Data from the ACESII mission organized by minorframe.150 words per minor frame.40 minor frames to a major frame.',
             'TEXT': 'Raw Data from the ACESII mission organized by minorframe.150 words per minor frame.40 minor frames to a major frame.'
             }, {
                'Source_name': f'ACESII_36364>ACES II 36.364',
                'Data_type': 'K0>Key Parameter',
                'PI_name': 'Bounds',
                'Logical_source': f'aces_36364_',
                'Logical_file_id': f'aces_36364_00000000_v01',
                'Logical_source_description': 'Raw Data from the ACESII mission organized by minorframe.150 words per minor frame.40 minor frames to a major frame.',
                'TEXT': 'Raw Data from the ACESII mission organized by minorframe.150 words per minor frame.40 minor frames to a major frame.'
            }],
        'Launch_Times': [722236869184000000, 722236969184000000], # TT2000 values corresponding to 17:20:00 and 17:21:40 for high/low flyer, respectively.
        'g_nWords': 150,
        'g_nNomRecLen': 150 * 2 + 12,
        'nNomDataLen': 150 - 2,
        'gSync': b'\x40\x28\x6b\xfe',
        # g_sSync = "\xfe\x6b\x28\x40"  is the real frame sync but the data is little endian
        'nDataLen': "74I",
        'InstrNames': ['EEPAA', 'LEESA', 'IEPAA', 'LP'],
        'InstrNames_LC': ['eepaa', 'leesa', 'iepaa', 'lp'],
        'NumOfInstr': 4,
        'mfPerFullSweep': 4,
        'num_of_sector_counts': [21, 21, 7, 8],
        'words_per_mf': [9, 9, 7, 9],
        'ESA_words': [[10, 26, 43, 60, 85, 101, 114, 126, 142],
                      [14, 28, 45, 61, 89, 103, 116, 132, 144],
                      [11, 27, 44, 86, 102, 115, 143]],
        'epoch_fillVal': -9223372036854775808,
        'esaMaxCounts': 4095,
        'ESA_CLK_INPUT': 625,
        'LP_words': [15, 29, 46, 66, 90, 104, 117, 133, 145],
        'LP_Variables': ["deltaNdivN","step", "ne_swept", "ni_swept","ni"],
        'MinorFrameTime': 250000,
        'LPSamplePeriod': 31250,
        'LPFixedProbeBias':[-5.05,-4.96],
        'LPFixed_calResistances':[ {'Open':165,'500M':2170,'250M':2308,'100M':2540,'50M':2710,'10M':3122,'5M':3299},
                                   {'Open':165,'500M':2063,'250M':2260,'100M':2500,'50M':2660,'10M':3060,'5M':3194}],
        'LPswept_voltage_range':[[-4.72,2.12],[-4.68,2.08]],
        'Epoch_range_to_determine_stepDAC': [[dt.datetime(2022, 11, 20,17,24,50,300),dt.datetime(2022, 11, 20,17,25,10,310)],
                                             [dt.datetime(2022, 11, 20,17,24,50,400),dt.datetime(2022, 11, 20,17,25,10,450)]],
        'InstrNames_Full': ['EEPAA>Electron Energy Pitch Angle Analyzer',
                            'LEESA>Low Energy Electrostatic Analyzer',
                            'IEPAA>Ion Energy Pitch Angle Analyzer',
                            'Langmuir_Probe>Langmuir Probe'],
        'geometric_factor': [
            [0.000174 for i in range(21)],
            [0.000174 for i in range(21)],
            [0.000174 for i in range(7)]],
        'deadtime':[2.e-07,2.e-07],
        'Instr_sector_to_pitch': [
            [-10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190],
            [-10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190],
            [180, 150, 120, 90, 60, 30, 0]
        ],
        'Instr_words_per_sweep': [49, 49, 49],
        'Instr_Energy': [ #eepaa
                        np.array(
                        [12750.00, 10905.65, 9328.10, 7978.75, 6824.58, 5837.38, 4992.97, 4270.72,
                          3652.94,  3124.52, 2672.55, 2285.95, 1955.28, 1672.44, 1430.51, 1223.58,
                          1046.58,   895.19,  765.70,  654.94,  560.20,  479.16,  409.85,  350.56,
                           299.85,   256.48,  219.38,  187.64,  160.50,  137.28,  117.42,  100.44,
                            85.91,    73.48,   62.85,   53.76,   45.98,   39.33,   33.64,   28.78,
                            24.61,    21.05,   18.01,   15.40,   13.17,   11.27,    9.64,    8.24, 7.05]),
                        # leesa
                        np.array(
                        [12750.00, 10905.65, 9328.10, 7978.75, 6824.58, 5837.38, 4992.97, 4270.72,
                          3652.94,  3124.52, 2672.55, 2285.95, 1955.28, 1672.44, 1430.51, 1223.58,
                          1046.58,   895.19,  765.70,  654.94,  560.20,  479.16,  409.85,  350.56,
                           299.85,   256.48,  219.38,  187.64,  160.50,  137.28,  117.42,  100.44,
                            85.91,    73.48,   62.85,   53.76,   45.98,   39.33,   33.64,   28.78,
                            24.61,    21.05,   18.01,   15.40,   13.17,   11.27,    9.64,    8.24, 7.05])/1000,
                        # iepaa
                        np.array(
                        [12750.00, 10905.65, 9328.10, 7978.75, 6824.58, 5837.38, 4992.97, 4270.72,
                         3652.94,   3124.52, 2672.55, 2285.95, 1955.28, 1672.44, 1430.51, 1223.58,
                         1046.58,    895.19,  765.70,  654.94,  560.20,  479.16,  409.85,  350.56,
                          299.85,    256.48,  219.38,  187.64,  160.50,  137.28,  117.42,  100.44,
                           85.91,     73.48,   62.85,   53.76,   45.98,   39.33,   33.64,   28.78,
                           24.61,     21.05,   18.01,   15.40,   13.17,   11.27,    9.64,    8.24, 7.05])
                        ]
    }

    ACESII_attrs = makeRocketAttrs(ACESII_attrs_dict)


    # Initialization Information for L0 Variables
    EEPAA_data_dict = {
        'Epoch':                    [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': ACESII_attrs.epoch_fillVal,'FORMAT': 'I5','UNITS': 'ns','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','MONOTON':'INCREASE','TIME_BASE':'J2000','TIME_SCALE':'Terrestrial Time','REFERENCE_POSITION':'Rotating Earth Geoid','SCALETYP':'linear'}],
        'STEPPER_Voltage_ADC':      [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'MCP_Voltage_ADC':          [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'MCP_Current_ADC':          [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'STACK_Voltage_ADC':        [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'STACK_Current_ADC':        [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'minor_frame_counter':      [[],                                             {'DEPEND_0':  None,  'DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'Count_Interval':           [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'major_frame_counter':      [[],                                             {'DEPEND_0':  None,  'DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'Sector_Counts':            [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': 'Sector_Number','DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','SCALETYP':'linear'}],
        'Sector_Number':            [[[i + 1 for i in (range(21))]],                 {'DEPEND_0':  None,  'DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'Pitch_Angle':              [[[(-10 + i*10) for i in (range(21))]],          {'DEPEND_0':  None,  'DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': 'deg','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear','LABLAXIS': 'Pitch_Angle'}],
        '625kHz_Clock_Input':       [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        '3p3V_Voltage_Monitor_ADC': [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        '5V_Voltage_Monitor_ADC':   [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'Temperature_Monitor_ADC':  [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'sfid':                     [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'EXP_Current':              [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        '28V_Monitor':              [[],                                             {'DEPEND_0': 'Epoch_monitors','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'Boom_Monitor':             [[],                                             {'DEPEND_0': 'Epoch_monitors','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'Epoch_monitors':           [[],                                             {'DEPEND_0': 'Epoch_monitors','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': ACESII_attrs.epoch_fillVal,'FORMAT': 'I5','UNITS': 'ns','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','MONOTON':'INCREASE','TIME_BASE':'J2000','TIME_SCALE':'Terrestrial Time','REFERENCE_POSITION':'Rotating Earth Geoid','SCALETYP':'linear'}],
        'sync_word':                [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
        'status_word_1':            [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
        'status_word_3':            [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
        'HV_div16':                 [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
        'HV_enable':                [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
        'sweep_step':               [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
        'TP5_enable':               [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]
    }
    for key, val in EEPAA_data_dict.items():
        EEPAA_data_dict[key][1]['LABLAXIS'] = key
        EEPAA_data_dict[key][1]['FIELDNAM'] = key
        EEPAA_data_dict[key][1]['CATDESC'] = key

    # -- LEESA ---
    LEESA_data_dict = deepcopy(EEPAA_data_dict)

    # -- IEPAA ---
    IEPAA_data_dict = deepcopy(EEPAA_data_dict)
    del IEPAA_data_dict['Sector_Counts'], IEPAA_data_dict['Sector_Number'], IEPAA_data_dict['Pitch_Angle']
    IEPAA_data_dict['Sector_Counts'] = [[], EEPAA_data_dict['Sector_Counts'][1]]  # IEPAA has 7 anodes with 30deg/anode totally 210deg coverage
    IEPAA_data_dict['Sector_Number'] = [[j + 1 for j in (range(ACESII_attrs.words_per_mf[2]))], EEPAA_data_dict['Sector_Number'][1]]  # IEPAA has 7 anodes with 30deg/anode totally 210deg coverage
    IEPAA_data_dict['Pitch_Angle'] = [[j*30 for j in range(ACESII_attrs.words_per_mf[2])], EEPAA_data_dict['Pitch_Angle'][1]]

    # -- LPs ---
    LP_data_dict = deepcopy(EEPAA_data_dict)
    LPremove = ['STEPPER_Voltage_ADC',    'MCP_Voltage_ADC',         'MCP_Current_ADC',
                'STACK_Voltage_ADC',      'STACK_Current_ADC',       'Sector_Counts',
                'Sector_Number',          'Pitch_Angle',             '3p3V_Voltage_Monitor_ADC',
                '5V_Voltage_Monitor_ADC', 'Temperature_Monitor_ADC',
                'Boom_Monitor',           'Count_Interval',          '625kHz_Clock_Input',
                'sync_word',              'status_word_1',
                'TP5_enable',             'status_word_3',           'HV_div16',
                'HV_enable',              'sweep_step']

    for thing in LPremove:
        del LP_data_dict[thing]

    LP_data_dict['Channel_Number'] =   [[[i for i in (range(len(ACESII_attrs_dict['LP_words'])))]], {'CATDESC': 'Channel_Number', 'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FIELDNAM': 'Channel_Number', 'FILLVAL': -1, 'FORMAT': 'I5', 'LABLAXIS': 'Channel_Number', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]
    LP_data_dict['Channel_Counts'] =   [[],                            {'CATDESC': 'Channel_Counts', 'DEPEND_0': 'Epoch', 'DEPEND_1': 'Channel_Number', 'DEPEND_2': None, 'FIELDNAM': 'Channel_Counts', 'FILLVAL': -1, 'FORMAT': 'I5','LABLAXIS': 'Channel_Counts', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None,'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]
    LP_data_dict['Boom_Monitor_1'] =   [[],                            {'CATDESC': 'Boom_Monitor_1', 'DEPEND_0': 'Epoch_monitor_1', 'DEPEND_1': None, 'DEPEND_2': None, 'FIELDNAM': 'Boom_Monitor_1', 'FILLVAL': -1, 'FORMAT': 'I5', 'LABLAXIS': 'Boom_Monitor_1', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]
    LP_data_dict['Boom_Monitor_2'] =   [[],                            {'CATDESC': 'Boom_Monitor_2', 'DEPEND_0': 'Epoch_monitor_2', 'DEPEND_1': None, 'DEPEND_2': None, 'FIELDNAM': 'Boom_Monitor_2', 'FILLVAL': -1, 'FORMAT': 'I5', 'LABLAXIS': 'Boom_Monitor_2', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]
    LP_data_dict['Epoch_monitor_1'] =  [[],                            {'CATDESC': 'Epoch_monitor_1','DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII_attrs.epoch_fillVal, 'FORMAT': 'I5', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time', 'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear'}]
    LP_data_dict['Epoch_monitor_2'] =  [[],                            {'CATDESC': 'Epoch_monitor_2','DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII_attrs.epoch_fillVal, 'FORMAT': 'I5', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time', 'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear'}]

    data_dicts = [EEPAA_data_dict, LEESA_data_dict, IEPAA_data_dict,LP_data_dict,LP_data_dict]



    # --- Deconvolution Keys for the ESAs ---
    deConvolveKey_EEPAA = {
        1: None,
        2: 'sync_word',
        3: 'status_word_1',
        4: 'status_word_2',
        5: 'status_word_3',
        6: 'STEPPER_Voltage_ADC',
        7: 'MCP_Voltage_ADC',
        8: 'MCP_Current_ADC',
        9: 'STACK_Voltage_ADC',
        10: 'STACK_Current_ADC',
        11: '3p3V_Voltage_Monitor_ADC',
        12: '5V_Voltage_Monitor_ADC',
        13: 'Temperature_Monitor_ADC',
        14: 'Sector_Counts',
        15: 'Sector_Counts',
        16: 'Sector_Counts',
        17: 'Sector_Counts',
        18: 'Sector_Counts',
        19: 'Sector_Counts',
        20: 'Sector_Counts',
        21: 'Sector_Counts',
        22: 'Sector_Counts',
        23: 'Sector_Counts',
        24: 'Sector_Counts',
        25: 'Sector_Counts',
        26: 'Sector_Counts',
        27: 'Sector_Counts',
        28: 'Sector_Counts',
        29: 'Sector_Counts',
        30: 'Sector_Counts',
        31: 'Sector_Counts',
        32: 'Sector_Counts',
        33: 'Sector_Counts',
        34: 'Sector_Counts',
        35: '625kHz_Clock_Input',
        36: None
    }
    deConvolveKey_LEESA = deepcopy(deConvolveKey_EEPAA)
    deConvolveKey_IEPAA = deepcopy(deConvolveKey_EEPAA)
    deConvoleKeys = [deConvolveKey_EEPAA,deConvolveKey_LEESA,deConvolveKey_IEPAA]
    for i in [21 + i for i in (range(16))]:
        del deConvolveKey_IEPAA[i]
    deConvolveKey_IEPAA[21] = '625kHz_Clock_Input'
    for i in [22,23,24,25,26,27,28]:
        deConvolveKey_IEPAA[i] = None

    return ACESII_attrs,data_dicts,deConvoleKeys



# --- --- --- --- --- ---
# --- TRICEII Mission ---
# --- --- --- --- --- ---
def TRICE_mission_dicts():
    TRICE_attrs_dict = {
        'missionNam' : 'TRICE',
        'rocketID': ['52003', '52004'],
        'g_nWords': 150,
        'g_nNomRecLen': 150 * 2 + 12,
        'nNomDataLen': 150 - 2,
        'gSync' : b'\x40\x28\x6b\xfe',
        'nDataLen': "59I",
        'InstrNames_LC': ['eepaa']
    }

    TRICEII_attrs = makeRocketAttrs(TRICE_attrs_dict)

    return TRICEII_attrs



