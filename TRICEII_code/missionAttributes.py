# --- missionAttributes.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to store the attributes specific to certain rockets

# --- --- --- --- --- ---
# --- ACES II Mission ---
# --- --- --- --- --- ---

class makeRocketAttrs:
    def __init__(self, missionAttrs):
        self.missionNam = missionAttrs['missionNam']
        self.rocketID = missionAttrs['rocketID']
        self.globalAttributes = missionAttrs['globalAttributes']
        self.InstrNames_Full = missionAttrs['InstrNames_Full']
        self.InstrNames_LC = missionAttrs['InstrNames_LC']
        self.fillVal = missionAttrs['fillVal']
        self.epoch_starts = missionAttrs['Epoch_Starts']
        self.geometric_factor = missionAttrs['geometric_factor']
        self.deadtime = missionAttrs['deadtime']


def TRICE_mission_dicts():

    missionAttributes = {
        'missionNam': 'TRICE',
        'rocketID': ['52003', '52004'],
        'fillVal':-1,
        'fillVal_epoch': -9223372036854775808,
        'Epoch_Starts': [[2018, 12, 8, 8, 26, 0, 000, 000, 0], [2018, 12, 8, 8, 28, 0, 000, 000, 0]],
        'geometric_factor':[
        [0.0001914 , 0.00014442, 0.000174  , 0.00016356, 0.00019314,
        0.00016182, 0.000174  , 0.00017748, 0.00018096, 0.00015834,
        0.00018792, 0.00017226, 0.00017748, 0.00017574, 0.00016704,
        0.00015834, 0.00018618, 0.00019488, 0.0001479 , 0.00016182,
        0.0001914 ],
        [0.000174, 0.000174, 0.000174, 0.000174, 0.000174, 0.000174,
         0.000174, 0.000174, 0.000174, 0.000174, 0.000174, 0.000174,
         0.000174, 0.000174, 0.000174, 0.000174, 0.000174, 0.000174,
         0.000174, 0.000174, 0.000174]
        ],
        'deadtime':[2.e-07,2.e-07],
        'g_nWords': 150,
        'g_nNomRecLen': 150 * 2 + 12,
        'nNomDataLen': 150 - 2,
        'gSync': b'\x40\x28\x6b\xfe',
        'nDataLen': "59I",
        'InstrNames_LC': 'eepaa',
        'InstrNames_Full':'EEPAA>Electron Energy Pitch Angle Analyzer',
        'globalAttributes': [
            {'Source_name': f'TRICE_53003>TRICE II High 52.003',
             'Data_type': 'K0>Key Parameter',
             'PI_name': 'Kletzing',
             'Logical_source': f'trice52003_',
             'Logical_file_id': f'trice52003_00000000_v01',
             'Logical_source_description': 'Data from the TRICEII mission',
             'TEXT': 'Data from the TRICEII mission'
             },
            {'Source_name': f'TRICE_53004>TRICE II Low 52.004',
             'Data_type': 'K0>Key Parameter',
             'PI_name': 'Kletzing',
             'Logical_source': f'trice52004_',
             'Logical_file_id': f'trice52004_00000000_v01',
             'Logical_source_description': 'Data from the TRICEII mission',
             'TEXT': 'Data from the TRICEII mission'
             }],
        'Launcher_Settings':
            [{'T0': '12/8/2018 8:26:00 GMT', 'Launcher': 'Pad 7, U3','Latitude': '69.2941 Deg',
             'Longitude': '16.0189 Deg', 'Altitude': '16.4ft', 'Elevation': '81.8Deg',
             'Azimuth': '12.3Deg', 'Bank': '0 deg'},
             {'T0': '12/8/2018 8:28:00 GMT', 'Launcher': 'Athena', 'Latitude': '69.2942 Deg',
             'Longitude': '16.0193 Deg', 'Altitude': '16.4ft', 'Elevation': '81.5Deg',
             'Azimuth': '12.8Deg', 'Bank': '0 deg'}
            ]
    }


    # Initialization Information for L0 Variables
    data_dict_temp = {
        'Epoch': [[], {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': missionAttributes['fillVal_epoch'],'FORMAT': 'I5', 'UNITS': 'ns','VALIDMIN':None, 'VALIDMAX': None, 'VAR_TYPE':'support_data','MONOTON':'INCREASE','TIME_BASE':'J2000','TIME_SCALE':'Terrestrial Time','REFERENCE_POSITION':'Rotating Earth Geoid','SCALETYP':'linear'}],
        'data':  [[], {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': missionAttributes['fillVal'], 'FORMAT': 'I5', 'UNITS': None ,'VALIDMIN':None, 'VALIDMAX': None, 'VAR_TYPE':'support_data', 'SCALETYP':'linear'}],
    }

    TRICE_attrs = makeRocketAttrs(missionAttributes)



    return TRICE_attrs,missionAttributes,data_dict_temp



