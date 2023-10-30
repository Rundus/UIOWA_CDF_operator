# --- TestParticle.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:
import matplotlib.pyplot as plt
import numpy as np

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.myImports import *
from scipy.stats import linregress
start_time = time.time()


colors = ['tab:purple','tab:orange','tab:red','tab:green','tab:blue']

# --- --- --- --- --- --- --
# --- SIMULATION TOGGLES ---
# --- --- --- --- --- --- --
Z0 = 6500 # initial altitude of particles (in km)
Z0_wave = 7500 # initial altitude of the wave (in km)
N = 64 # number of particles in each energy
ptcl_mass = m_e # mass of the particle
ptcl_charge = q0 # charge of the particle
Energies = [10, 20, 30, 40, 50] # choose some energies
simLen = 10 # how many delta T steps to simulate
deltaT = 0.1  # in seconds
Amplitude = 0.0012  # V/m
lamda = (7500 - 3500)  # in km
WaveSpeed = 2000  # in km/s
inverse = False

def testParticle_Sim():

    def E_Field(simTime, Alt):

        sign = -1 if inverse else 1
        deltaPos = simTime * WaveSpeed

        if Alt > Z0_wave - deltaPos:  # if we're above the wave
            E = 0
        elif (Z0_wave - deltaPos - 0.05 * lamda) <= Alt < (Z0_wave - deltaPos):  # if we're in the first risetime
            slope, intercept, r, p, se = linregress(x=[Z0_wave - deltaPos, Z0_wave - deltaPos - 0.05 * lamda], y=[0, Amplitude])
            E = sign * (slope * Alt + intercept)
        elif (Z0_wave - deltaPos - 0.45 * lamda) <= Alt < (Z0_wave - deltaPos - 0.05 * lamda):  # if we're in the first max amplitude part
            E = sign * Amplitude
        elif (Z0_wave - deltaPos - 0.55 * lamda) <= Alt < (Z0_wave - deltaPos - 0.45 * lamda):  # if we're in the middle transition region
            slope, intercept, r, p, se = linregress(x=[Z0_wave - deltaPos - 0.45 * lamda, Z0_wave - deltaPos - 0.55 * lamda], y=[Amplitude, -1 * Amplitude])
            E = sign * (slope * Alt + intercept)
        elif (Z0_wave - deltaPos - 0.95 * lamda) <= Alt < (Z0_wave - deltaPos - 0.55 * lamda):  # if we're in the lowest max amplitude part
            E = -1 * sign * Amplitude
        elif (Z0_wave - deltaPos - lamda) <= Alt < (Z0_wave - deltaPos - 0.95 * lamda):  # if we're in the last risetime
            slope, intercept, r, p, se = linregress(x=[Z0_wave - deltaPos - 0.95 * lamda, Z0_wave - deltaPos - lamda], y=[-1 * Amplitude, 0])
            E = sign * (slope * Alt + intercept)
        elif Alt < (Z0_wave - deltaPos - lamda):  # if we're outside the wave on the bottom edge
            E = 0

        return E

    def forceFunc(simTime, Alt):
        return (ptcl_charge / ptcl_mass) * E_Field(simTime, Alt)

    def AdamsBashforth(yn1, funcVal_n, funcVal_n1):
        yn2 = yn1 + (3 / 2) * deltaT * funcVal_n1 - (1 / 2) * deltaT * funcVal_n
        return yn2

    def Euler(yn, funcVal_n):
        yn1 = yn + deltaT * funcVal_n
        return yn1


    # --- --- --- --- --- --- --- ----
    # --- INITIALIZE THE PARTICLES ---
    # --- --- --- --- --- --- --- ----

    # Generate initial position and velocities
    prgMsg('Generating Particles')
    ptcl_pos = np.array([[Z0 for i in range(N)] for engy in Energies ])
    ptcl_vel_par = []
    ptcl_vel_perp = []

    for engy in Energies:
        V_mag = np.sqrt(2*engy*ptcl_charge/ptcl_mass)
        pitches = np.radians(np.linspace(0, 360, N))
        ptcl_vel_par.append([V_mag*np.cos(ptch)/1000 for ptch in pitches])
        ptcl_vel_perp.append([V_mag*np.sin(ptch)/1000 for ptch in pitches])

    ptcl_vel_perp = np.array(ptcl_vel_perp)
    ptcl_vel_par = np.array(ptcl_vel_par)

    # create the data dict
    data_dict = {}
    for i,engy in enumerate(Energies):
        data_dict = {**data_dict, **{f"{engy}_vpar":[ptcl_vel_par[i]]}}
        data_dict = {**data_dict, **{f"{engy}_vperp": [ptcl_vel_perp[i]]}}
        data_dict = {**data_dict, **{f"{engy}_zpos":[ptcl_pos[i]]}}

    Done(start_time)


    # # Show the Energies
    # fig, ax = plt.subplots()
    # for i in range(len(Energies)):
    #     ax.scatter(data_dict[f'{Energies[i]}_vperp'],data_dict[f'{Energies[i]}_vpar'],color=colors[i])
    # plt.show()

    # --- --- --- --- --- --- -
    # --- IMPLEMENT THE SIM ---
    # --- --- --- --- --- --- -
    Alt = np.linspace(1000,8000,1000)
    simTime = [deltaT*i for i in range(simLen)]



    # --- IMPLEMENT 1 STEP OF EULER ---
    for engy in Energies: # number of energies

        vPars = data_dict[f"{engy}_vpar"][0]
        zPos = data_dict[f"{engy}_zpos"][0]

        newVpars = []
        newZpos = []

        # Determine the new velocities
        for i in range(N): # numer of particles

            newVel = Euler(yn=vPars[i], funcVal_n=forceFunc(simTime[0], zPos[i])/1000 )
            newVpars.append(newVel)
            newPos = Euler(yn= zPos[i], funcVal_n=newVel)
            newZpos.append(newPos)

        # store the new data
        data_dict[f"{engy}_vperp"].append(data_dict[f"{engy}_vperp"][0])
        data_dict[f"{engy}_vpar"].append(newVpars)
        data_dict[f"{engy}_zpos"].append(newZpos)

    fig, ax = plt.subplots(2)
    for i in range(2):
        for j,engy in enumerate(Energies):
            ax[i].scatter(data_dict[f"{engy}_vperp"][i], data_dict[f"{engy}_vpar"][i],color=colors[j])
    plt.show()

    # --- IMPLEMENT ADAMS BASHFORTH---
    # for j,tme in enumerate(simTime[2:]):
    #
    #     for engy in Energies: # number of energies
    #
    #         vPars_n = data_dict[f"{engy}_vpar"][0]
    #         zPos_n = data_dict[f"{engy}_zpos"][0]
    #
    #         vPars_n1 = data_dict[f"{engy}_vpar"][0]
    #         zPos_n1 = data_dict[f"{engy}_zpos"][0]
    #
    #         newVpars = []
    #         newZpos = []
    #
    #         # Determine the new velocities
    #         for i in range(N): # numer of particles
    #
    #             newVel = Euler(yn=vPars[i], funcVal_n=forceFunc(simTime[0], zPos[i]) )
    #             newVpars.append(newVel)
    #             newPos = Euler(yn= zPos[i], funcVal_n=newVel)
    #             newZpos.append(newPos)
    #
    #         # store the new data
    #         data_dict[f"{engy}_vperp"].append(data_dict[f"{engy}_vperp"][0])
    #         data_dict[f"{engy}_vpar"].append(newVpars)
    #         data_dict[f"{engy}_zpos"].append(newZpos)









    def ForceFunc(E, mass, charge):
        return (charge / mass) * E

    def AdamsBashforth(yn1, yn, deltaT, func):
        yn2 = yn1 + (3 / 2) * deltaT * func(yn1) - (1 / 2) * deltaT * func(yn)
        return yn2

    def Euler(yn, deltaT, func):
        yn1 = yn + deltaT * func(yn)
        return yn1


    return




# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

testParticle_Sim()