import cdflib.epochs

print('Importing Variables: ',end='')
import numpy as np
from cdflib import cdfread,cdfwrite
import csv
from Variables import EEPAA_file_high_info,EEPAA_file_low_info,atts_high,atts_low,EEPAA_file_high,EEPAA_file_low
from functions import write_var_to_file_attitude
from files import root,output
print('Done')

atti_filenam_high = '52003Kletzing_timeSynchronized.csv'
atti_filenam_low = '52004Kletzing_timeSynchronized.csv'

csvfiles = [root + atti_filenam_low, root + atti_filenam_high]

atti_control_dat = [[],[]] #Keep the [low_dat,high_dat] structure thoughout all lists
varnames = [[],[]]
listatti_control_dat = []

#-----------------------
#Store the data in lists
#-----------------------
for l in range(2):
    with open(csvfiles[l]) as file:
        reader = csv.reader(file)
        for row in reader:
            atti_control_dat[l].append(row)

# Prepare the data in separated column formats
for i in range(2):
    listdata = atti_control_dat[i]
    varnames[i] = listdata[4]
    varnames[i].append('Epoch')

    listatti_control_dat.append([[] for j in range(27)])

    if i == 0:
        row_counter = 0
        for row in atti_control_dat[i]:
            if row_counter >=5:
                for k in range(27):
                    listatti_control_dat[i][k].append(row[k])

            row_counter +=1

    if i == 1:
        row_counter = 0
        for row in atti_control_dat[i]:
            if row_counter >= 5:
                for k in range(27):
                    listatti_control_dat[i][k].append(row[k])

            row_counter += 1

#-------------------------------------
# Convert time-since-launch into epoch
#-------------------------------------

# -----------------------
# OUTPUT DATA TO CDF FILE
# -----------------------

Launcher_Settings = [ {'T0':'12/8/2018 8:28:00 GMT','Launcher': 'Athena', 'Latitude': '69.2942 Deg', 'Longitude':'16.0193 Deg','Altitude':'16.4ft','Elevation':'81.5Deg','Azimuth':'12.8Deg','Bank':'0 deg'},{'T0':'12/8/2018 8:26:00 GMT','Launcher': 'Pad 7, U3', 'Latitude': '69.2941 Deg', 'Longitude':'16.0189 Deg','Altitude':'16.4ft','Elevation':'81.8Deg','Azimuth':'12.3Deg','Bank':'0 deg'}  ]

#----------------------------------
#File information and file creation
#----------------------------------

#HIGH
output_paths = [root+ output + 'TRICE_52004_20181208T082243_attitude_control&_v1.1.1',root+ output + 'TRICE_52003_20181208T082239_attitude_control&_v1.1.1']
file_info = EEPAA_file_high_info
file_info['CDF'],file_info['zVariables'],file_info['Version'] = output_paths[1],varnames[1],'1.1.1'

#Output the launcher settings to attributes
file_attrs = file_info['Attributes']
file_attrs.append(Launcher_Settings[1])
file_info['Attributes'] = file_attrs

# print(file_info)
attitude_control_high_file = cdfwrite.CDF(root+ output + 'TRICE_52003_20181208T082239_attitude_control&_v1.1.1', delete=True,cdf_spec=file_info) #high file


#LOW
file_info = EEPAA_file_low_info
file_info['CDF'],file_info['zVariables'],file_info['Version'] = output_paths[0],varnames[0],'1.1.1'

#Output the launcher settings to attributes
file_attrs = file_info['Attributes']
file_attrs.append(Launcher_Settings[0])
file_info['Attributes'] = file_attrs

attitude_control_low_file = cdfwrite.CDF(root+ output + 'TRICE_52004_20181208T082243_attitude_control&_v1.1.1', delete=True,cdf_spec=file_info)#low file

output_files = [attitude_control_low_file,attitude_control_high_file]
Epochs_starts = [[2018,12,8,8,28,0,000,000,0],[2018,12,8,8,26,0,000,000,0]]



# ------------------------------
# WRITE THE VARIABLE INFORMATION
# ------------------------------


# attrs = [namestr,fieldnam_str (should be same as namestr), fillval, data_min(numpy dtype float32),data_max (numpy dtype float32),scaletype_str]
# infos = [which_variable_number,datatype,length of data,data description,pads]


x_axis_depend = 'Epoch'
for i in range(2): #Time

    #%%%%%time%%%%%
    varnum = 0
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])

    attrs = ['Time_since_launch',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','seconds','series']
    infos = [varnum, 44, len(vardata), 'Time_Since_Launch',[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    varattributes['VAR_TYPE'] = 'Support_data'
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%Roll%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Roll (Euler)',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%pitch%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Pitch (Euler)',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%RollRate%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Roll_Rate',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','Hz','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%AoA%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Angle_of_Attack',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #%%%%%a11%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Row1Col1_cosinemat_ENU',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    varattributes['VAR_TYPE'] = 'Support_data'
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #%%%%%a12%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Row1Col2_cosinemat_ENU',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    varattributes['VAR_TYPE'] = 'Support_data'
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #%%%%%a13%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Row1Col3_cosinemat_ENU',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    varattributes['VAR_TYPE'] = 'Support_data'
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #%%%%%a21%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Row2Col1_cosinemat_ENU',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    varattributes['VAR_TYPE'] = 'Support_data'
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #%%%%%a22%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Row2Col2_cosinemat_ENU',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    varattributes['VAR_TYPE'] = 'Support_data'
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #%%%%%a23%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Row2Col3_cosinemat_ENU',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    varattributes['VAR_TYPE'] = 'Support_data'
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%a31%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Row3Col1_cosinemat_ENU',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    varattributes['VAR_TYPE'] = 'Support_data'
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%a32%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Row3Col2_cosinemat_ENU',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    varattributes['VAR_TYPE'] = 'Support_data'
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #%%%%%a33%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Row3Col3_cosinemat_ENU',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    varattributes['VAR_TYPE'] = 'Support_data'
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #%%%%%X_Az%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['X_Azimuth_ENU_fixedframe',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%X_EL%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['X_Elevation_ENU_fixedframe',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%X_Az%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Y_Azimuth_ENU_fixedframe',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%Y_EL%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Y_Elevation_ENU_fixedframe',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #%%%%%Z_Az%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Z_Azimuth_ENU_fixedframe',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%Z_EL%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Z_Elevation_ENU_fixedframe',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%F_az%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Payload_VelocityVector_Azimuth_ENU_fixedframe',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #%%%%%F_el%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Payload_VelocityVector_Elevation_ENU_fixedframe',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%latg%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Latitude_Geodetic',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)


    #%%%%%longg%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Longitude_Geodetic',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','deg','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #%%%%%alt%%%%%
    varnum += 1
    vardata = np.array([float(dat) for dat in listatti_control_dat[i][varnum]])
    attrs = ['Altitude',[-1.e+31],[vardata.min()],[vardata.max()], 'linear','meters','series']
    infos = [varnum, 44, len(vardata), attrs[0],[-9223372036854775807]]

    varattributes = {'CATDESC': attrs[0], 'DEPEND_0': x_axis_depend, 'DISPLAY_TYPE': attrs[6], 'FIELDNAM': attrs[0], 'FILLVAL':np.array(attrs[1],dtype='float32'), 'FORMAT': 'E12.2', 'LABLAXIS': attrs[0], 'LABL_PTR_1': attrs[5], 'LABL_PTR_2': 'eepaa_LABL_2', 'UNITS': attrs[5], 'VALIDMIN': attrs[2], 'VALIDMAX': attrs[3], 'VAR_TYPE': 'data', 'SCALETYP':attrs[4] }
    varinfo = {'Variable': attrs[0], 'Num': infos[0], 'Var_Type': 'zVariable', 'Data_Type': infos[1], 'Data_Type_Description': infos[3], 'Num_Elements': 1, 'Num_Dims': 0, 'Dim_Sizes': [], 'Sparse': 'No_sparse', 'Last_Rec': infos[2], 'Rec_Vary': True, 'Dim_Vary': [], 'Pad': np.array(infos[4],dtype='float32'), 'Compress': 0, 'Block_Factor': 0}
    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)

    #%%%%%Epoch%%%%%
    varnum += 1
    #Produce the epoch using the given excel file's time tags
    startval = cdflib.epochs.CDFepoch.compute_tt2000(Epochs_starts[i])
    vardata = np.array([(float(value)*(10**9) + startval) for value in listatti_control_dat[i][0]],dtype='float64')
    attrs = ['Epoch', [-1.e+31], [vardata.min()], [vardata.max()], 'linear', 'ns', 'series']
    infos = [varnum, 33, len(vardata), attrs[0], [-9223372036854775807]]

    epochattrs = EEPAA_file_high.varattsget('Epoch')
    epochattrs['VALIDMIN'] = attrs[2]
    epochattrs['VALIDMAX'] = attrs[3]
    epochattrs['VAR_TYPE'] = 'data'

    varinfoa = EEPAA_file_high.varinq('Epoch')
    varinfoa['Last_Rec'] =infos[2]

    varattributes = epochattrs
    varattributes['DEPEND_0'] = 'Epoch'
    varattributes['CATDESC'] = 'Epoch'
    varinfo = varinfoa


    print(varattributes)
    print(varinfo)

    output_files[i].write_var(varinfo, var_attrs=varattributes, var_data=vardata)



print('Done')



