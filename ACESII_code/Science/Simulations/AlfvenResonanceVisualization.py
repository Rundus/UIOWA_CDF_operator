# --- imports ---
import matplotlib.pyplot as plt
import numpy as np
from numpy import linspace
from numpy import array
from myspaceToolsLib.physicsVariables import m_to_km,q0,m_e
from myspaceToolsLib.CDF_load import outputCDFdata


# --- SIMULATION TOGGLES ---

simLen = 550 # how many delta T steps to simulate
deltaT = 0.1 # in seconds
alt_Rez = 2000 # number of points in the altitude grid
simAltLow = 200*m_to_km # low altitude (in meters)
simAltHigh = 20000*m_to_km # high altitude (in meters)

simAlt = linspace(simAltLow, simAltHigh, alt_Rez)  # in METERS
simTime = linspace(0, deltaT * simLen, simLen + 1)  # in seconds

# Wave
wavefreq = 0.5 # in Hz
tau0 = 1/wavefreq
waveSpeed = 1000*m_to_km
lambdaPara = waveSpeed/wavefreq  # in km
Z0Start = 8000*m_to_km - (lambdaPara/2) # initial position of parallel E-Field
EparaAmp = 0.0001

# PARTICLES
elec1_Z0 = 5000*m_to_km  # in km
elect1_V0 = waveSpeed*0.3 # ABSOULUTE VALUE (sign doesn't matteR)
elec2_Z0 = 12500*m_to_km # in km
elect2_V0 = waveSpeed*2 # ABSOULUTE VALUE (sign doesn't matteR)




# --- OUTPUT ---
outputPath = r'C:\Data\Simulations\AlfvenResonanceVisualization.cdf'
data_dict_sim = {'EField_wf':[[],{'DEPEND_0': 'simTime', 'DEPEND_1': 'simAlt', 'UNITS': 'V/m', 'LABLAXIS': 'Epara'}],
                 'Phi_wf':[[],{'DEPEND_0': 'simTime', 'DEPEND_1': 'simAlt', 'UNITS': 'eV', 'LABLAXIS': 'Wave Potential'}],
                 'elec_wf_x':[[],{'DEPEND_0': 'simTime', 'DEPEND_1': 'simAlt', 'UNITS': 'm', 'LABLAXIS': 'Position'}],
                 'elec_wf_v':[[],{'DEPEND_0': 'simTime', 'DEPEND_1': 'simAlt', 'UNITS': 'm/s', 'LABLAXIS': 'Velocity'}],
                 'accel_wf':[[],{'DEPEND_0': 'simTime', 'DEPEND_1': 'simAlt', 'UNITS': 'm/s^2', 'LABLAXIS': 'Accelertion_LabFrame'}],
                 'EField_lf':[[],{'DEPEND_0': 'simTime', 'DEPEND_1': 'simAlt', 'UNITS': 'V/m', 'LABLAXIS': 'Epara'}],
                 'Phi_lf':[[],{'DEPEND_0': 'simTime', 'DEPEND_1': 'simAlt', 'UNITS': 'eV', 'LABLAXIS': 'Wave Potential'}],
                 'elec_lf_x':[[],{'DEPEND_0': 'simTime', 'DEPEND_1': 'simAlt', 'UNITS': 'x', 'LABLAXIS': 'Position'}],
                 'elec_lf_v':[[],{'DEPEND_0': 'simTime', 'DEPEND_1': 'simAlt', 'UNITS': 'm/s', 'LABLAXIS': 'Velocity'}],
                 'accel_lf':[[],{'DEPEND_0': 'simTime', 'DEPEND_1': 'simAlt', 'UNITS': 'm/s^2', 'LABLAXIS': 'acceleration_LabFrame'}],
                 'simAlt':[simAlt,{'DEPEND_0': 'simAlt', 'UNITS': 'm', 'LABLAXIS': 'simAlt'}],
                 'simTime':[simTime,{'DEPEND_0': 'simAlt', 'UNITS': 'seconds', 'LABLAXIS': 'simTime'}]}


# --- Plot toggles ---
Plot_Initial_E = True

def Epara_generator(z, t, Vel, initialZ):
    # the middle of wave
    if initialZ - Vel * (t - tau0) > z > initialZ - Vel * t:
        EperpVal = EparaAmp * (np.sin(((z - initialZ) + Vel * t) * (2 * np.pi / (Vel * tau0))))
    else:
        EperpVal = 0

    return EperpVal

def Euler(yn, funcVal_n):
    yn1 = yn + deltaT * funcVal_n
    return yn1

def AdamsBashforth(yn1, funcVal_n, funcVal_n1):
    yn2 = yn1 + (3 / 2) * deltaT * funcVal_n1 - (1 / 2) * deltaT * funcVal_n
    return yn2

#%%%%%%%%%%%%%%%%%%%%%%%%
#### LABRATORY FRAME ####
#%%%%%%%%%%%%%%%%%%%%%%%%
# --- Initialize Parallel E-Field ---

# populate the Epara Data
for t in simTime:
    data_dict_sim['EField_lf'][0].append([Epara_generator(z=val,t=t, initialZ=Z0Start,Vel=waveSpeed) for val in simAlt])
    data_dict_sim['EField_wf'][0].append([Epara_generator(z=val, t=t, initialZ=Z0Start, Vel=0) for val in simAlt])


# --- populate the electrons ---
# INITALIZE - Lab Frame
data_dict_sim['elec_lf_x'][0].append([elec1_Z0, elec2_Z0])
data_dict_sim['elec_lf_v'][0].append([-1*np.abs(elect1_V0), -1*np.abs(elect2_V0)])
EparaInit = data_dict_sim['EField_lf'][0][0]
data_dict_sim['accel_lf'][0].append([(-q0/m_e)*EparaInit[np.abs(simAlt - data_dict_sim['elec_lf_x'][0][0][0]).argmin()],
                                     (-q0/m_e)*EparaInit[np.abs(simAlt - data_dict_sim['elec_lf_x'][0][0][1]).argmin()]])

# INITALIZE - Wave Frame
data_dict_sim['elec_wf_x'][0].append([elec1_Z0, elec2_Z0])
data_dict_sim['elec_wf_v'][0].append([-1*(elect1_V0 - waveSpeed), -1*(elect2_V0 - waveSpeed)])
EparaInit = data_dict_sim['EField_wf'][0][0]
data_dict_sim['accel_wf'][0].append([(-q0/m_e)*EparaInit[np.abs(simAlt - data_dict_sim['elec_wf_x'][0][0][0]).argmin()],
                                     (-q0/m_e)*EparaInit[np.abs(simAlt - data_dict_sim['elec_wf_x'][0][0][1]).argmin()]])

# EULER - Lab Frame
data_dict_sim['elec_lf_v'][0].append([Euler(yn=data_dict_sim['elec_lf_v'][0][0][0],funcVal_n=data_dict_sim['accel_lf'][0][0][0]), Euler(yn=data_dict_sim['elec_lf_v'][0][0][1],funcVal_n=data_dict_sim['accel_lf'][0][0][1])])
data_dict_sim['elec_lf_x'][0].append([Euler(yn=data_dict_sim['elec_lf_x'][0][0][0],funcVal_n=data_dict_sim['elec_lf_v'][0][0][0]), Euler(yn=data_dict_sim['elec_lf_x'][0][0][1],funcVal_n=data_dict_sim['elec_lf_v'][0][0][1])])
Epara_Euler = data_dict_sim['EField_lf'][0][1]
data_dict_sim['accel_lf'][0].append([(-q0/m_e)*Epara_Euler[np.abs(simAlt -  data_dict_sim['elec_lf_x'][0][1][ptclIdx]).argmin()] for ptclIdx in range(2)])

# EULER - Wave Frame
data_dict_sim['elec_wf_v'][0].append(  [Euler(yn=data_dict_sim['elec_wf_v'][0][0][0],funcVal_n=data_dict_sim['accel_wf'][0][0][0]), Euler(yn=data_dict_sim['elec_wf_v'][0][0][1],funcVal_n=data_dict_sim['accel_wf'][0][0][1])])
data_dict_sim['elec_wf_x'][0].append([Euler(yn=data_dict_sim['elec_wf_x'][0][0][0],funcVal_n=data_dict_sim['elec_wf_v'][0][0][0]), Euler(yn=data_dict_sim['elec_wf_x'][0][0][1],funcVal_n=data_dict_sim['elec_wf_v'][0][0][1])])
Epara_Euler = data_dict_sim['EField_wf'][0][1]
data_dict_sim['accel_wf'][0].append([(-q0/m_e)*Epara_Euler[np.abs(simAlt -  data_dict_sim['elec_wf_x'][0][1][ptclIdx]).argmin()] for ptclIdx in range(2)])

for idx, t in enumerate(simTime):

    if idx > 1:
        # ---------------------------
        # Adams Bashforth - Lab Frame
        # ---------------------------
        Epara_n = data_dict_sim['EField_lf'][0][idx-2]
        elec_x_n = data_dict_sim['elec_lf_x'][0][idx-2]
        elec_v_n = data_dict_sim['elec_lf_v'][0][idx-2]
        accel_n = data_dict_sim['accel_lf'][0][idx-2]
        Epara_n1 = data_dict_sim['EField_lf'][0][idx - 1]
        elec_x_n1 = data_dict_sim['elec_lf_x'][0][idx - 1]
        elec_v_n1 = data_dict_sim['elec_lf_v'][0][idx - 1]
        accel_n1 = data_dict_sim['accel_lf'][0][idx - 1]
        Epara_n2 = data_dict_sim['EField_lf'][0][idx]

        # determine new velocity/position
        elec_v_n2 = [AdamsBashforth(yn1=elec_v_n1[ptclIdx],funcVal_n=accel_n[ptclIdx],funcVal_n1=accel_n1[ptclIdx]) for ptclIdx in range(len(elec_x_n))]
        elec_x_n2 = [AdamsBashforth(yn1=elec_x_n1[ptclIdx],funcVal_n=elec_v_n[ptclIdx],funcVal_n1=elec_v_n1[ptclIdx]) for ptclIdx in range(len(elec_x_n))]
        data_dict_sim['elec_lf_v'][0].append(elec_v_n2)
        data_dict_sim['elec_lf_x'][0].append(elec_x_n2)

        # determine new acceleration
        data_dict_sim['accel_lf'][0].append([(-q0/m_e)*Epara_n2[np.abs(simAlt -  elec_x_n2[ptclIdx]).argmin()] for ptclIdx in range(len(elec_x_n))])

        # ----------------------------
        # Adams Bashforth - Wave Frame
        # ----------------------------
        # get previous data
        Epara_n = data_dict_sim['EField_lf'][0][idx - 2]
        elec_x_n = data_dict_sim['elec_lf_x'][0][idx - 2]
        elec_v_n = data_dict_sim['elec_lf_v'][0][idx - 2]
        accel_n = data_dict_sim['accel_lf'][0][idx - 2]
        Epara_n1 = data_dict_sim['EField_lf'][0][idx - 1]
        elec_x_n1 = data_dict_sim['elec_lf_x'][0][idx - 1]
        elec_v_n1 = data_dict_sim['elec_lf_v'][0][idx - 1]
        accel_n1 = data_dict_sim['accel_lf'][0][idx - 1]
        Epara_n2 = data_dict_sim['EField_lf'][0][idx]

        # determine new velocity/position
        elec_v_n2 = [AdamsBashforth(yn1=elec_v_n1[ptclIdx], funcVal_n=accel_n[ptclIdx], funcVal_n1=accel_n1[ptclIdx]) for ptclIdx in range(len(elec_x_n))]
        elec_x_n2 = [AdamsBashforth(yn1=elec_x_n1[ptclIdx], funcVal_n=elec_v_n[ptclIdx], funcVal_n1=elec_v_n1[ptclIdx]) for ptclIdx in range(len(elec_x_n))]
        data_dict_sim['elec_lf_v'][0].append(elec_v_n2)
        data_dict_sim['elec_lf_x'][0].append(elec_x_n2)

        # determine new acceleration
        data_dict_sim['accel_lf'][0].append([(-q0 / m_e) * Epara_n2[np.abs(simAlt - elec_x_n2[ptclIdx]).argmin()] for ptclIdx in range(len(elec_x_n))])



for t in range(len(simTime)):
    fig, ax = plt.subplots()
    ax.plot(simAlt, data_dict_sim['EField_wf'][0][t])
    ax.scatter(data_dict_sim['elec_wf_x'][0][t][0],0,color='red',s=40)
    ax.scatter(data_dict_sim['elec_wf_x'][0][t][1], 0, color='green',s=40)
    ax.set_xlim(0,simAltHigh)
    plt.show()

#####################
# --- OUTPUT DATA ---
#####################
outputCDFdata(outputPath=outputPath,data_dict=data_dict_sim)



