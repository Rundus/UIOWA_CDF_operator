# --- TestParticle.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Test particle simulaton that was recreated using ideas from a Knudsen and Wu paper:
# https://doi. org/10.1029/2020JA028005

# TODO: Use ion concentrations to better estimate rho_m in alfven velocity

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.myImports import *
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from scipy.stats import linregress
from ACESII_code.class_var_func import CHAOS, lightSpeed, u0
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import matplotlib.patches as mpatches
start_time = time.time()


############################
# --- --- --- --- --- --- --
# --- SIMULATION TOGGLES ---
# --- --- --- --- --- --- --
############################
m_to_km = 1E3

# --- simToggles ---
simLen = 330 # how many delta T steps to simulate
deltaT = 0.01 # in seconds
simAltLow = 500*m_to_km # low altitude in meters
simAltHigh = 12000*m_to_km # high altitude in meters

# --- particle toggles ---
Z0 = (1*6371)*m_to_km # initial altitude of particles (in meters)
N = 64 # number of particles in each energy
ptcl_mass = m_e # mass of the particle
ptcl_charge = q0 # charge of the particle

# The range of energies for each color (the rightmost linspace point is ignored)
simEnergyRanges = [[0.01, 1],
                   [1, 5],
                   [5, 10],
                   [10, 15],
                   [15, 20],
                   [20, 25],
                  [25, 30]
                   ]

# the color choice for each set of particles
simColors = [['tab:purple'],
             ['tab:orange'],
             ['tab:red'],
             ['tab:blue'],
             ['tab:green'],
             ['tab:brown'],
             ['tab:pink']]

# the number of Energies to linspace the energy ranges
simPtclsPerEngyRange = [5 for i in range(len(simEnergyRanges))]

# which pitch angles should the initial partcles be at
customPitches = True
simPtcls_Initial_Pitch_Range = (list(np.linspace(-25, 25, int(1*N/8))) +
                                list(np.linspace(155, 205, int(5*N/8)))+
                                list(np.linspace(25, 155, int(N/8)))+
                                list(np.linspace(205, 335, int(N/8))))


# --- WAVE E-Field Toggles ---
flipEField = True
dynamicWaveSpeed = False
Z0_wave = (1*6371+3000)*m_to_km # initial altitude of the wave (in meters)

# static
Amplitude_static = 0.15/1000  # V/m
lamda_static = (7500 - 4500)*m_to_km  # in meters
WaveSpeed_Static = 2500*m_to_km  # in meters/s
duty = 40 # duty cycle percent (Max is 50)

# dynamic
Eperp = 12E-3 # V/m
alpha = 1/10 # defines the ratio of the wave's wavelength to the total simulated altitude size
kperp_obs = 1/alpha

# --- B-Field Toggles ---
Lshell = 9

# --- observation toggles ---
obsHeight = 400*m_to_km # in meters

# --- diagnostic plotting ---
plotEField = False
plotMirror_vs_electrostatic = False
plotInputEnergies = False
plotEuler_1step = False

#### MAJOR TOGGLES ###
# --- Particle Trajectory Movies ---
outputAnimation = True
fps = 30
# -------------------------
simulated_EEPAA_plots = True
pitches_to_plot = [i for i in range(21)] # indicies of pitch bins that will be plotted from [-10, 0 , 10, 20, 30, ...]
# -------------------------
simulated_pitchAnglePlot_EEPAA_animation = True
fps_ptichAnglePlot = 15
# -------------------------
outputObservedParticleData = True




outputTitle = 'TestParticle\TestParticle_above_invTrue.mp4' if flipEField else 'TestParticle\TestParticle_above_invFalse.mp4'

def testParticle_Sim():

    ##########################
    # --- SIMULATION START ---
    ##########################
    Alt = np.linspace(simAltLow, simAltHigh, 2000)  # in METERS
    simTime = [deltaT * i for i in range(simLen)]

    # determine the Energies and colors used in the simulation
    Energies = []
    colors = []
    for i in range(len(simEnergyRanges)):
        temp = np.linspace(simEnergyRanges[i][0], simEnergyRanges[i][1],simPtclsPerEngyRange[i])
        for j,val in enumerate(temp):
            if j != len(temp)-1:
                Energies.append(val)
                colors.append(simColors[i][0])

    # --- --- --- --- --- --- --- --- --- --
    # --- GENERATE THE GEOMAGNETIC FIELD ---
    # --- --- --- --- --- --- --- --- --- --
    prgMsg('Getting Geomagnetic Field')

    R_REF = 6371.2 * m_to_km  # in meters
    geomagAlts = [((alt + R_REF) / R_REF) for alt in Alt]
    geomagLats = np.array([np.degrees(np.arccos(radi / Lshell)) for radi in geomagAlts])
    geomagLongs = np.array([111.83 for i in range(len(Alt))])
    times = [dt.datetime(2022, 11, 20, 17, 20, 00, 000) for i in range(len(Alt))]
    Pos = np.array([geomagAlts, geomagLats, geomagLongs]).transpose()
    ISOtime = [times[i].isoformat() for i in range(len(times))]
    cvals_MAG = coord.Coords(Pos, 'MAG', 'sph')
    cvals_MAG.ticks = Ticktock(ISOtime, 'ISO')
    cvals_GDZ = cvals_MAG.convert('GEO', 'sph')
    Lat_geo = cvals_GDZ.lati

    # Get the Chaos model
    B = CHAOS(Lat_geo, [15.25 for i in range(len(Alt))], np.array(Alt) / m_to_km, times)
    Bmag = (1E-9) * np.array([np.linalg.norm(Bvec) for Bvec in B])
    Bgrad = [(Bmag[i + 1] - Bmag[i]) / (Alt[i + 1] - Alt[i]) for i in range(len(Bmag) - 1)] + [0]
    Done(start_time)

    # --- --- --- --- --- --- --- --- ---
    # --- GENERATE THE ELECTRIC FIELD ---
    # --- --- --- --- --- --- --- --- ---

    # --- determine the density over all altitudes ---
    # Description: returns density for altitude "z [km]" in m^-3
    h = 0.06 * (R_REF / m_to_km)  # in km from E's surface
    n0 = 6E4
    n1 = 1.34E7
    z0 = 0.05 * (R_REF / m_to_km)  # in km from E's surface
    n = [(cm_to_m ** 3) * (n0 * np.exp(-1 * ((alt / m_to_km) - z0) / h) + n1 * ((alt / m_to_km) ** (-1.55))) for alt in Alt]  # calculated density (in cm^-3)

    # --- determine the electron plasma density and skin depth ---
    plasmaFreq = [np.sqrt((n[j] * q0 * q0) / (m_e * ep0)) for j in range(len(n))]
    skinDepth = [lightSpeed / freq for freq in plasmaFreq]

    # --- determine kperp throughout all altitudes ---
    kperp = [kperp_obs * np.sqrt(Bmag[j] / Bmag[0]) for j in range(len(Alt))]

    # --- determine MHD alfven speed ---
    VA_MHD = [Bmag[j]/np.sqrt(u0*IonMasses[0]*n[j]) for j in range(len(Alt))]

    # Desciption: the output variable is a[ [len(Alt)], ... len(simtime)] object that describes the electric field at
    # every altitude for each time in the simulation

    E_Field = []

    for k, tme in enumerate(simTime):

        if duty > 50:
            raise Exception('Duty cannot be more than 50%')

        E_temp= []
        for j, alt in enumerate(Alt):

            # determine the wave parameters
            sign = -1 if flipEField else 1 # flip the sign of the amplitude
            WaveSpeed = 5 if dynamicWaveSpeed else WaveSpeed_Static
            lamda = 5 if dynamicWaveSpeed else lamda_static
            deltaPos = 5 if dynamicWaveSpeed else tme * WaveSpeed_Static
            Amplitude = 5 if dynamicWaveSpeed else Amplitude_static

            riseTime = lamda*(50 - duty)/(2*100)
            dutyTime = lamda*(duty/100)

            if alt > Z0_wave - deltaPos:  # if we're above the wave
                E = 0
            elif (Z0_wave - deltaPos - riseTime) <= alt < (Z0_wave - deltaPos):  # if we're in the first risetime
                slope, intercept, r, p, se = linregress(x=[Z0_wave - deltaPos, Z0_wave - deltaPos - riseTime], y=[0, Amplitude])
                E = sign * (slope * alt + intercept)
            elif (Z0_wave - deltaPos - (dutyTime+riseTime)) <= alt < (Z0_wave - deltaPos - riseTime):  # if we're in the first max amplitude part
                E = sign * Amplitude
            elif (Z0_wave - deltaPos - (dutyTime+3*riseTime)) <= alt < (Z0_wave - deltaPos - (dutyTime+riseTime)):  # if we're in the middle transition region
                slope, intercept, r, p, se = linregress(x=[Z0_wave - deltaPos - (dutyTime+riseTime), Z0_wave - deltaPos - (dutyTime+3*riseTime)], y=[Amplitude, -1 * Amplitude])
                E = sign * (slope * alt + intercept)
            elif (Z0_wave - deltaPos - (2*dutyTime+3*riseTime)) <= alt < (Z0_wave - deltaPos - (dutyTime+3*riseTime)):  # if we're in the lowest max amplitude part
                E = -1 * sign * Amplitude
            elif (Z0_wave - deltaPos - lamda) <= alt < (Z0_wave - deltaPos - (2*dutyTime+3*riseTime)):  # if we're in the last risetime
                slope, intercept, r, p, se = linregress(x=[Z0_wave - deltaPos - (2*dutyTime+3*riseTime), Z0_wave - deltaPos - lamda], y=[-1 * Amplitude, 0])
                E = sign * (slope * alt + intercept)
            elif alt < (Z0_wave - deltaPos - lamda):  # if we're outside the wave on the bottom edge
                E = 0
            else:
                E = 0

            E_temp.append(E)
        E_Field.append(E_temp)

    #--- --- --- --- --- --- --- -
    # --- SIMULATION FUNCTIONS ---
    #--- --- --- --- --- --- --- -

    def forceFunc(simTime_Index, alt_Index, mu, deltaB): # gives meters/s^2
        # return (-1 * ptcl_charge / ptcl_mass) * E_Field(simTime, Alt)
        # return  - (mu / ptcl_mass) * deltaB
        return (-1*ptcl_charge / ptcl_mass) * E_Field[simTime_Index][alt_Index] - (mu/ptcl_mass)*deltaB

    def AdamsBashforth(yn1, funcVal_n, funcVal_n1):
        yn2 = yn1 + (3 / 2) * deltaT * funcVal_n1 - (1 / 2) * deltaT * funcVal_n
        return yn2

    def Euler(yn, funcVal_n):
        yn1 = yn + deltaT * funcVal_n
        return yn1

    def calcE(Vperp,Vpar,mass,charge):
        return 0.5*(mass/charge)*(Vperp**2 + Vpar**2)

        # return 0.5 * (mass / charge) * (Vpar ** 2)


    # --- --- --- --- --- --- --- ----
    # --- INITIALIZE THE PARTICLES ---
    # --- --- --- --- --- --- --- ----
    prgMsg('Generating Particles')

    # --- create an empty data_dict ---
    data_dict = {}
    varNames = ['vpar','vperp','zpos','force','moment','Bmag','Bgrad','observed','energies']
    for i in range(len(Energies)):
        data_dict = {**data_dict, **{f"{Energies[i]}_{vNam}": [] for vNam in varNames}}

    # --- Populate Initial Variables ---
    for i, engy in enumerate(Energies):
        V_mag = np.sqrt(2 * engy * ptcl_charge / ptcl_mass)

        pitches = np.radians(np.linspace(0, 360, N)) if not customPitches else np.radians(np.array(simPtcls_Initial_Pitch_Range))


        data_dict[f'{Energies[i]}_zpos'].append([Z0 for ptch in pitches]) # z position
        data_dict[f'{Energies[i]}_vpar'].append([V_mag * np.cos(ptch) for ptch in pitches]) # parallel velocity
        data_dict[f'{Energies[i]}_vperp'].append([V_mag * np.sin(ptch) for ptch in pitches]) # perpendicular velocity
        data_dict[f'{Energies[i]}_energies'].append([calcE(V_mag * np.sin(ptch), V_mag * np.cos(ptch), ptcl_mass, ptcl_charge) for ptch in pitches]) # particle energy
        data_dict[f'{Energies[i]}_Bmag'].append([Bmag[np.abs(Alt-zpos).argmin()] for zpos in data_dict[f'{Energies[i]}_zpos'][0]]) # determine B-Field
        data_dict[f'{Energies[i]}_Bgrad'].append([Bgrad[np.abs(Alt-zpos).argmin()] for zpos in data_dict[f'{Energies[i]}_zpos'][0]]) # determine deltaB
        data_dict[f'{Energies[i]}_moment'].append([0.5*m_e*(data_dict[f'{Energies[i]}_vperp'][0][k]**2)/data_dict[f'{Energies[i]}_Bmag'][0][k] for k in range(N)])  # determine magnetic moment
        data_dict[f'{Energies[i]}_force'].append([forceFunc(simTime_Index=0, alt_Index=np.abs(Alt - data_dict[f'{Energies[i]}_zpos'][0][k]).argmin(), mu=data_dict[f'{Energies[i]}_moment'][0][k], deltaB=data_dict[f'{Energies[i]}_Bgrad'][0][k]) for k in range(N)]) # determine the initial force function
        data_dict[f'{Energies[i]}_observed'].append([1 if data_dict[f'{Energies[i]}_zpos'][0][k] <= obsHeight else 0 for k in range(N)])
    Done(start_time)

    # --- --- --- --- --- --- --- --- -
    # --- IMPLEMENT 1 STEP OF EULER ---
    # --- --- --- --- --- --- --- --- -
    prgMsg('1 Step of Euler')
    for engy in Energies: # number of energies

        # Get the "n" data
        vPars0 = data_dict[f"{engy}_vpar"][0]
        vPerps0 = data_dict[f"{engy}_vperp"][0]
        zPos0 = data_dict[f"{engy}_zpos"][0]
        mus0 = data_dict[f"{engy}_moment"][0]
        deltaBs0 = data_dict[f"{engy}_Bgrad"][0]
        Bmags0 = data_dict[f"{engy}_Bmag"][0]

        # add new list to each of the variables
        for vNam in varNames:
            data_dict[f"{engy}_{vNam}"].append([])

        # Determine the "n+1" data
        for k in range(N): # numer of particles

            # determine new parallel velocity
            newVpar = Euler(yn=vPars0[k], funcVal_n=forceFunc(simTime_Index=0, alt_Index=np.abs(Alt - zPos0[k]).argmin(), deltaB=deltaBs0[k], mu=mus0[k]))
            data_dict[f"{engy}_vpar"][-1].append(newVpar)

            # determine new z-position
            newPos = Euler(yn=zPos0[k], funcVal_n=newVpar)
            data_dict[f"{engy}_zpos"][-1].append(newPos)

            # determine new B
            newB = Bmag[np.abs(Alt-newPos).argmin()]
            data_dict[f"{engy}_Bmag"][-1].append(newB)

            # determine new Vperp
            newVperp = vPerps0[k] * np.sqrt(newB/Bmags0[k])
            data_dict[f"{engy}_vperp"][-1].append(newVperp)

            # determine new deltaB
            newdB = Bgrad[np.abs(Alt-newPos).argmin()]
            data_dict[f"{engy}_Bgrad"][-1].append(newdB)

            # determine new moments
            newMu = mus0[k]
            data_dict[f"{engy}_moment"][-1].append(newMu)

            # determine new force
            data_dict[f"{engy}_force"][-1].append(forceFunc(simTime_Index=1, alt_Index= np.abs(Alt - newPos).argmin(), deltaB=newdB, mu=newMu))

            # determine if new particle is observed
            data_dict[f"{engy}_observed"][-1].append(1) if newPos <= obsHeight else data_dict[f"{engy}_observed"][-1].append(0)

            # determine new particle energies
            data_dict[f"{engy}_energies"][-1].append(calcE(newVperp,newVpar,ptcl_mass,ptcl_charge))
    Done(start_time)

    # --- --- --- --- --- --- --- --- -
    # --- IMPLEMENT ADAMS BASHFORTH ---
    # --- --- --- --- --- --- --- --- -
    prgMsg('Adams Bashforth')
    print('\n')
    for j, tme in tqdm(enumerate(simTime)):

        if tme > deltaT: # only start after 2 timesteps (we have inital T0 and used Euler to get T1)

            # loop through all the energies
            for engy in Energies:

                vPerp_initial = data_dict[f"{engy}_vperp"][0]
                Bmag_initial = data_dict[f"{engy}_Bmag"][0]

                vPars_n = data_dict[f"{engy}_vpar"][j-2]
                zPos_n = data_dict[f"{engy}_zpos"][j-2]
                force_n = data_dict[f"{engy}_force"][j-2]
                mu_n = data_dict[f"{engy}_moment"][j-2]
                Bmag_n = data_dict[f"{engy}_Bmag"][j-2]
                deltaB_n = data_dict[f"{engy}_Bgrad"][j - 2]

                vPars_n1 = data_dict[f"{engy}_vpar"][j - 1]
                zPos_n1 = data_dict[f"{engy}_zpos"][j - 1]
                force_n1 = data_dict[f"{engy}_force"][j - 1]
                mu_n1 = data_dict[f"{engy}_moment"][j - 1]
                Bmag_n1 = data_dict[f"{engy}_Bmag"][j - 1]
                deltaB_n1 = data_dict[f"{engy}_Bgrad"][j - 1]

                # add new list to each of the variables
                for vNam in varNames:
                    data_dict[f"{engy}_{vNam}"].append([])

                # for each particle in a particluar energy
                for i in range(N):

                    # new parallel velocity
                    newVpar = AdamsBashforth(yn1=vPars_n1[i], funcVal_n=force_n[i], funcVal_n1=force_n1[i])
                    data_dict[f"{engy}_vpar"][-1].append(newVpar)

                    # new zPosition
                    newPos = AdamsBashforth(yn1=zPos_n1[i], funcVal_n=vPars_n[i], funcVal_n1=vPars_n1[i])
                    data_dict[f"{engy}_zpos"][-1].append(newPos)

                    # new Brad
                    newdB = Bgrad[np.abs(newPos - Alt).argmin()]
                    data_dict[f"{engy}_Bgrad"][-1].append(newdB)

                    # new Force
                    data_dict[f"{engy}_force"][-1].append(forceFunc(simTime_Index=j, alt_Index=np.abs(Alt - newPos).argmin(), deltaB=newdB, mu=mu_n[i],))

                    # new mu
                    data_dict[f"{engy}_moment"][-1].append(mu_n[i])

                    # new Bmag
                    newBmag = Bmag[np.abs(newPos-Alt).argmin()]
                    data_dict[f"{engy}_Bmag"][-1].append(newBmag)

                    # new Vperp
                    newVperp = vPerp_initial[i]*np.sqrt(newBmag/Bmag_initial[i])
                    data_dict[f"{engy}_vperp"][-1].append(newVperp)

                    # determine if new particle is observed
                    data_dict[f"{engy}_observed"][-1].append(1) if newPos <= obsHeight else data_dict[f"{engy}_observed"][-1].append(0)

                    # determine new partcle energy
                    data_dict[f"{engy}_energies"][-1].append(calcE(newVperp,newVpar,ptcl_mass,ptcl_charge))

    # --- CHECK THE OBSERVED PARTICLES ---
    # description: Look through the observed data and find when a particle if/was observed. If it was
    # observed (obs==1), find the first time and then set the rest of the observed variable==0. Also, find the highest ptcl energy
    # that is observed (used later)
    obsPtclData = []
    TEmaxE = 0
    for tme in range(simLen):
        for engy in Energies:
            observed = np.array(data_dict[f'{engy}_observed'][tme])
            for ptclN in range(N):
                if observed[ptclN] == 1:

                    # get the observed particle information
                    obs_ptcl_engy = data_dict[f'{engy}_energies'][tme][ptclN]
                    obs_ptcl_vperp = data_dict[f'{engy}_vperp'][tme][ptclN]
                    obs_ptcl_vpar = data_dict[f'{engy}_vpar'][tme][ptclN]
                    obs_ptcl_pitch = np.degrees(np.arccos(-1 * obs_ptcl_vpar / np.sqrt(obs_ptcl_vpar**2 + obs_ptcl_vperp**2))) if obs_ptcl_vperp > 0 else -1*np.degrees(np.arccos(-1 * obs_ptcl_vpar / np.sqrt(obs_ptcl_vpar**2 + obs_ptcl_vperp**2)))

                    # record the particle information
                    obsPtclData.append([simTime[tme], obs_ptcl_engy, obs_ptcl_pitch])

                    # find the largest observed Parallel energy
                    # TEmaxE = obs_ptcl_engy*np.cos(np.radians(obs_ptcl_pitch)) if obs_ptcl_engy*np.cos(np.radians(obs_ptcl_pitch)) > TEmaxE else TEmaxE # Epar attempt
                    TEmaxE = obs_ptcl_engy  if obs_ptcl_engy  > TEmaxE else TEmaxE

                    # record that we've now observed this particle and don't need to plot it again
                    for s in range(tme+1, simLen):
                        data_dict[f'{engy}_observed'][s][ptclN] = 0
    Done(start_time)
    # --- --- --- --- ----
    # --- DIAGNOSTICS ---
    # --- --- --- --- ----

    if plotInputEnergies:
        # Show the Energies
        fig, ax = plt.subplots()
        for i in range(len(Energies)):
            ax.scatter(data_dict[f'{Energies[i]}_vperp'],data_dict[f'{Energies[i]}_vpar'],color=colors[i])
        plt.show()

    if plotEField:
        fig, ax = plt.subplots()
        wTimes = [i for i in range(0, len(simTime), int(len(simTime)/len(colors)))]
        for i,index in enumerate(wTimes):
            ax.plot(E_Field[index],Alt, color= f'{colors[i]}', label='$\Delta t = $' + f'{wTimes[i]}')
        fig.legend()
        plt.show()

    if plotMirror_vs_electrostatic:

        chosenPitch = [0, 45, 90]
        fig, ax = plt.subplots(len(chosenPitch))

        for k in range(len(chosenPitch)):
            scaling = 1E8
            ax[k].plot([ptcl_charge*(Amplitude)/ (scaling*ptcl_mass) for i in range(len(Alt))], Alt/m_to_km,label='E-Force')
            ax[k].plot([(0.5*((np.sin(np.radians(chosenPitch[k]))*np.sqrt(2*10*ptcl_charge/ptcl_mass) )**2)*Bgrad[i])/(Bmag[i]*scaling) for i in range(len(Alt))]  ,Alt/m_to_km,label=r'$\nabla$ B Force')
            ax[k].legend()
            ax[k].set_title(f'{chosenPitch[k]}$^\circ$')
            ax[k].set_xlabel(f'Force [{scaling} N]')
            ax[k].set_ylabel('Altitude [m]')
        plt.show()

    if plotEuler_1step:
        fig, ax = plt.subplots(2)
        fig.suptitle('1 step of Euler')
        for i in range(2):
            for j,engy in enumerate(Energies):
                ax[i].scatter(np.array(data_dict[f"{engy}_vperp"][i])/m_to_km, np.array(data_dict[f"{engy}_vpar"][i])/m_to_km,color=colors[j])
            if i == 1:
                ax[i].set_title(f'$\Delta t$ = {simTime[1]}')
            elif i ==0:
                ax[i].set_title(f'$\Delta t$ = {simTime[0]}')
            ax[i].set_ylabel('vpar [km/s]')
            ax[i].set_ylabel('vperp [km/s]')
        plt.show()


    # --- --- --- --- --- --- --- --- --- --- --- ----
    # --- CREATE OBSERVED PARTICLE DATAST VARIABLE ---
    # --- --- --- --- --- --- --- --- --- --- --- ----

    # get the ACESII information
    rocketAttrs, b, c = ACES_mission_dicts()
    instr_engy_bins = rocketAttrs.Instr_Energy[0][0:41]
    instr_ptch_bins = rocketAttrs.Instr_sector_to_pitch[0]
    VperpGrid = np.zeros(shape=(len(instr_ptch_bins), len(instr_engy_bins)))
    VparaGrid = np.zeros(shape=(len(instr_ptch_bins), len(instr_engy_bins)))

    # create the meshgrids for plotting
    timeGrid, EngyGrid = np.meshgrid(simTime, instr_engy_bins)

    for ptch in range(len(instr_ptch_bins)):
        for engy in range(len(instr_engy_bins)):
            Emag = np.sqrt(2 * q0 * instr_engy_bins[engy] / (m_e))
            VperpGrid[ptch][engy] = np.sin(np.radians(instr_ptch_bins[ptch])) * Emag
            VparaGrid[ptch][engy] = np.cos(np.radians(instr_ptch_bins[ptch])) * Emag
    VparaGrid, VperpGrid = np.array(VparaGrid) / 1000, np.array(VperpGrid) / 1000

    # create the empty dataset variable
    simCounts = np.zeros(shape=(len(simTime), len(instr_ptch_bins), len(instr_engy_bins)))

    # fill in the counts data
    for obsPtcl in obsPtclData:

        # see if ptcl's ptch is within the detector. If so, then count it
        if -15 <= obsPtcl[2] <= 195:
            # determine where the particle should go in the counts variable
            obsTime_index = np.abs(np.array(simTime) - obsPtcl[0]).argmin()
            obsPtch_index = np.abs(np.array(instr_ptch_bins) - obsPtcl[2]).argmin()
            obsEngy_index = np.abs(np.array(instr_engy_bins) - obsPtcl[1]).argmin()

            # increment the counts value
            simCounts[obsTime_index][obsPtch_index][obsEngy_index] += 1

    # --- create a differential Energy Flux variable and fill it in for all time ---
    geo_factor = rocketAttrs.geometric_factor[0]
    countInterval = 0.8992 * 1E-3
    deadtime = 674E-9
    diffEFlux = np.zeros(shape=(len(simTime), len(instr_ptch_bins), len(instr_engy_bins)))

    for tme, ptch, engy in itertools.product(*[range(len(simTime)), range(len(instr_ptch_bins)), range(len(instr_engy_bins))]):
        measuredT = (countInterval) - (simCounts[tme][ptch][engy] * deadtime)
        diffEFlux[tme][ptch][engy] = int((simCounts[tme][ptch][engy]) / (geo_factor[ptch] * measuredT))


    # --- --- --- --- ---
    # --- OUTPUT DATA ---
    # --- --- --- --- ---

    if outputObservedParticleData:
        prgMsg('outputting Data')

        # get some model data
        globalAttrsMod = rocketAttrs.globalAttributes[0]
        globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L1'
        data_dict_counts = loadDictFromFile('C:\Data\ACESII\L1\high\ACESII_36359_l1_eepaa_fullCal.cdf',wKeys=['eepaa','Energy','Pitch_Angle','Epoch'])
        data_dict_diffFlux = loadDictFromFile('C:\Data\ACESII\L2\high\ACESII_36359_l2_eepaa_fullCal.cdf', wKeys=['Differential_Energy_Flux'])

        # create the output data_dict
        data_dict_output = {'Energy':deepcopy(data_dict_counts['Energy']),
                     'Pitch_Angle':deepcopy(data_dict_counts['Pitch_Angle']),
                     'simCounts':[simCounts,deepcopy(data_dict_counts['eepaa'][1])],
                     'simDiffEFLux':[diffEFlux,deepcopy(data_dict_diffFlux['Differential_Energy_Flux'][1])]}

        for key,val in data_dict_output.items():
            data_dict_output[key][1]['DEPEND_0'] = 'Time'

        # handle the Epoch variable
        exampleEpoch = {'DEPEND_0': 'Time', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                        'FORMAT': None, 'UNITS': 's', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data',
                        'MONOTON': 'INCREASE', 'TIME_BASE': None, 'TIME_SCALE': None,
                        'REFERENCE_POSITION': None, 'SCALETYP': 'linear', 'LABLAXIS': 'Time'}
        data_dict_output = {**data_dict_output, **{'Time':[np.array(simTime),exampleEpoch]}}


        # output the data
        outputCDFdata(outputPath='C:\Data\ACESII\science\simulations\TestParticle\TestParticleData.cdf',
                      data_dict=data_dict_output,
                      ModelData=L1_ACES_Quick(0),
                      globalAttrsMod=globalAttrsMod,
                      instrNam='EEPAA')
        Done(start_time)
    # --- --- --- --- ----
    # --- ANIMATE PLOT ---
    # --- --- --- --- ----

    if outputAnimation:
        prgMsg('Creating Animation')
        print('\n')
        #################
        # --- ANIMATE ---
        #################
        fig = plt.figure()
        figure_height = 17
        figure_width = 15
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)

        # increase all the plot label sizes
        plt.rcParams['axes.labelsize'] = 20
        plt.rcParams['axes.titlesize'] = 30
        # increase the tick label sizes
        plt.rcParams['xtick.labelsize'] = 15
        plt.rcParams['ytick.labelsize'] = 15

        # --- --- --- --- --
        # --- INITIALIZE ---
        # --- --- --- --- --
        title = fig.suptitle(f'$\Delta t$= {round(simTime[0],2)}\n' +
                             r'$E_{\parallel}$ = ' + f'{round(Amplitude*1000,3)} [mV/m]',fontsize=30)

        # gridspec
        gs0 = gridspec.GridSpec(nrows=5, ncols=2, figure=fig,hspace=0.5,wspace=0.3)

        # Altitude Plot
        axAlt = fig.add_subplot(gs0[0:3, 0])
        axAlt.axhline(y=obsHeight/m_to_km, linestyle='--')
        axAlt.text(x=-180,y=(1.2)*obsHeight/m_to_km, s=f'{obsHeight/m_to_km} km',fontsize=10)
        axAlt.set_xticks([-180 + 60*i for i in range(7)])
        axAlt.set_ylim(0, simAltHigh)
        axAlt.set_xlabel('PA $[^{\circ}$]')
        axAlt.set_ylabel('Alt [km]')

        # Show the E-Field
        axE = axAlt.twiny()
        E_Field_Wave, = axE.plot(E_Field[0], Alt/m_to_km, color='black')

        TopArrowX = -1*Amplitude / 2 if flipEField else Amplitude/2
        TopArrowY = (Z0_wave - (0.45*lamda))/m_to_km if flipEField else (Z0_wave - (0.05*lamda))/m_to_km
        BotArrowX = Amplitude / 2 if flipEField else -1*Amplitude/2
        BotArrowY = (Z0_wave - ((1-0.45)*lamda))/m_to_km if flipEField else (Z0_wave - ((1-0.05)*lamda))/m_to_km
        dx = 0
        dyTop = 0.35*lamda/m_to_km if flipEField else -0.35*lamda/m_to_km
        dyBot = -0.35*lamda/m_to_km if flipEField else 0.35*lamda/m_to_km

        topAr = axE.annotate("E$_{\parallel}$", xy=(TopArrowX, TopArrowY), xytext=(TopArrowX + dx, TopArrowY + dyTop),horizontalalignment="center", arrowprops=dict(arrowstyle="->",lw=5),fontsize=17)
        botAr = axE.annotate("E$_{\parallel}$", xy=(BotArrowX, BotArrowY), xytext=(BotArrowX + dx, BotArrowY + dyBot),horizontalalignment="center", arrowprops=dict(arrowstyle="->",lw=5),fontsize=17)

        axE.set_ylim(100, 8000)
        axE.set_xlim(-1.5*Amplitude, 1.5*Amplitude)
        axE.axis('off')

        # Vspace Plot
        axVspace = fig.add_subplot(gs0[0:3, 1])
        axVspace.set_ylabel('Vpar [km/s]')
        axVspace.set_xlabel('Vperp [km/s]')
        axVspace.set_xlim(-5E3, 5E3)
        axVspace.set_ylim(-5E3, 5E3)

        # Enegy Plot vs Time
        axTspaceE = fig.add_subplot(gs0[3, 0:2])
        TEplot = axTspaceE.scatter([],[])
        axTspaceE.set_xlim(0, simTime[-1])
        axTspaceE.set_ylim(20, 1.2*TEmaxE) if TEmaxE != 0 else axTspaceE.set_ylim(20, 500)
        axTspaceE.set_ylabel('Energy [eV]')
        axTspaceE.set_xlabel('Time [s]')

        patches = [mpatches.Patch(color=simColors[i][0], label=f'<{simEnergyRanges[i][1]} eV') for i in range(len(simEnergyRanges))]
        axTspaceE.legend(handles=patches)

        # Time vs Pitch Angle
        axTspacePitch = fig.add_subplot(gs0[4, 0:2])
        TPplot = axTspacePitch.scatter([], [])
        axTspacePitch.set_xlim(0, simTime[-1])
        axTspacePitch.set_yticks([-90 + 45 * i for i in range(5)])
        axTspacePitch.set_ylabel('Pitch Angle [$^\circ$]')
        axTspacePitch.set_xlabel('Time [s]')

        # --- INITIALIZE THE PLOTS ---
        AltitudeArtists = []
        VspaceArtists = []
        ptcls_to_plot = []

        for i in range(len(Energies)):

            vperps = np.array(data_dict[f'{Energies[i]}_vperp'][0]) / m_to_km
            vpars = np.array(data_dict[f'{Energies[i]}_vpar'][0]) / m_to_km

            Pitches = np.array([np.degrees(np.arccos(-1*vpars[k] /np.sqrt(vpars[k]**2 + vperps[k]**2))) if vperps[k] > 0 else -1*np.degrees(np.arccos(-1*vpars[k] /np.sqrt(vpars[k]**2 + vperps[k]**2))) for k in range(len(vpars))])
            zpos = np.array(data_dict[f'{Energies[i]}_zpos'][0]) / m_to_km

            # Plot altitudes
            AltitudeArtists.append(axAlt.scatter(Pitches, zpos, color=colors[i]))

            # plot energies
            VspaceArtists.append(axVspace.scatter(vperps, vpars, color=colors[i]))

        # --- Animation Function ---
        def animate_func(j):
            print('', end='\r' + color.RED + f'{round(100 * j / simLen, 1)} %' + color.END)

            # update the title
            title.set_text(f'$\Delta t$= {round(simTime[j],2)}\n'+
                           r'$|E_{\parallel}|$ = ' + f'{round(Amplitude*1000,3)} [mV/m]\n' +
                           f'$\omega$/k = {WaveSpeed_Static/m_to_km} km/s')
            axAlt.set_xticks([-180 + 60 * i for i in range(7)])

            # update the electric field
            E_Field_Wave.set_xdata(E_Field[j])
            E_Field_Wave.set_ydata(Alt/m_to_km)

            topAr.xy = (TopArrowX, (TopArrowY - WaveSpeed * simTime[j] / m_to_km))
            topAr.set_position((TopArrowX + dx, (TopArrowY - WaveSpeed * simTime[j] / m_to_km) + dyTop))
            botAr.xy = (BotArrowX, (BotArrowY - WaveSpeed * simTime[j] / m_to_km))
            botAr.set_position((BotArrowX + dx, (BotArrowY - WaveSpeed * simTime[j] / m_to_km) + dyBot))

            # Show the particles updated
            perpMax, perpMin = 0, 0
            parMax, parMin = 0, 0

            for i in range(len(Energies)):
                vperps = np.array(data_dict[f'{Energies[i]}_vperp'][j]) / m_to_km
                vpars = np.array(data_dict[f'{Energies[i]}_vpar'][j]) // m_to_km
                Pitches = np.array([np.degrees(np.arccos(-1*vpars[k] /np.sqrt(vpars[k]**2 + vperps[k]**2))) if vperps[k] > 0 else -1*np.degrees(np.arccos(-1*vpars[k] /np.sqrt(vpars[k]**2 + vperps[k]**2))) for k in range(len(vpars))])
                zpos = np.array(data_dict[f'{Energies[i]}_zpos'][j]) / m_to_km
                observed = np.array(data_dict[f'{Energies[i]}_observed'][j])

                # update axes limits of Vspace plot
                perpMin = vperps.min() if vperps.min() < perpMin else perpMin
                perpMax = vperps.max() if vperps.max() > perpMax else perpMax
                parMin = vpars.min() if vpars.min() < parMin else parMin
                parMax = vpars.max() if vperps.max() > parMax else parMax

                # update altitudes
                AltitudeArtists[i].set_offsets(np.stack([Pitches, zpos]).T)

                # update energies
                VspaceArtists[i].set_offsets(np.stack([vperps, vpars]).T)

                # update the Energy vs Time plot and Pitch vs Time
                for k in range(N):
                    if observed[k] == 1:
                        energy = 0.5 * (ptcl_mass / ptcl_charge) * ((vperps[k] * m_to_km) ** 2 + (vpars[k] * m_to_km) ** 2)
                        ptcls_to_plot.append([simTime[j], energy, Pitches[k], colors[i]])

            # update axes limits of Vspace plot
            axVspace.set_xlim(perpMin*1.2, perpMax*1.2)
            axVspace.set_ylim(parMin*1.2, parMax*1.2)

            # update the Energy vs Time plot
            if len(ptcls_to_plot) != 0:

                plotData = [[ptcls_to_plot[d][0] for d in range(len(ptcls_to_plot))],
                            [ptcls_to_plot[d][1] for d in range(len(ptcls_to_plot))],
                            [ptcls_to_plot[d][2] for d in range(len(ptcls_to_plot))],
                            [ptcls_to_plot[d][3] for d in range(len(ptcls_to_plot))]]

                # update the Energy_parallel vs Time Plot (time energy ptch color)
                # TEplot.set_offsets(np.stack([plotData[0], plotData[1]*np.cos(np.radians(np.array(plotData[1])))]).T) # Epar attempt
                TEplot.set_offsets(np.stack([plotData[0], plotData[1]]).T)
                TEplot.set_facecolor(plotData[3])

                # update the Pitch Angle vs Time Plot
                TPplot.set_offsets(np.stack([plotData[0], plotData[2]]).T)
                TPplot.set_facecolor(plotData[3])

        anim = animation.FuncAnimation(fig=fig, func=animate_func, frames=simLen, interval=1000 / fps)
        print('', end='\r' + color.RED + f'{round(100 * simLen / simLen, 1)} %' + color.END)
        anim.save(f'C:\Data\ACESII\science\simulations\{outputTitle}', fps=fps)
        print('\n')
        Done(start_time)

    # --- --- --- --- --- --- --- --
    # --- SIMULATE EEPAA RESULTS ---
    # --- --- --- --- --- --- --- --

    if simulated_EEPAA_plots:

        prgMsg('Simulating EEPAA Response')

        # --- --- --- ----
        # --- PLOTTING ---
        # --- --- --- ----

        for plotIndex in pitches_to_plot:

            # determine the data to plot
            plotThisData = np.array(simCounts[:, plotIndex, :])
            # plotThisData = plotThisData/plotThisData.max()

            # figure
            fig, ax = plt.subplots()

            # colormesh
            cmapPlot = ax.pcolormesh(timeGrid, EngyGrid, plotThisData.T, cmap='turbo', shading='auto',vmin=0,vmax=10)
            ax.set_ylabel('Energy [eV]')
            ax.set_xlabel('Time [s]')
            fig.suptitle('Simulated ACES II EEPAA\n'
                         r'$\alpha = $' + f'{instr_ptch_bins[plotIndex]}' + '$^{\circ}$')

            ax.set_ylim(min(instr_engy_bins), TEmaxE*1.5)

            # colorbar
            cbar = plt.colorbar(cmapPlot)
            cbar.set_label('Counts')

            plt.tight_layout()
            plt.savefig(f'C:\Data\ACESII\science\simulations\TestParticle\TestParticle_{instr_ptch_bins[plotIndex]}deg')
            plt.close()


        Done(start_time)

    if simulated_pitchAnglePlot_EEPAA_animation:

        # --- --- --- --- --- --- --
        # --- START THE ANIMATION ---
        # --- --- --- --- --- --- --
        prgMsg('Animating EEPAA Polar Response')

        # --- determine the plot's limits ---
        wEngyLim = 25 # corresponds to
        EnergyLimit = instr_engy_bins[wEngyLim]
        VelLimit = np.sqrt(2 * EnergyLimit * q0 / (m_e)) / 1000

        # --- INITIALIZE THE PLOT ---
        fig, ax = plt.subplots()
        figure_height = 15
        figure_width = 10
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)
        title = ax.set_title(f'$\Delta t$= {round(simTime[0],2)}\n'+
                           r'$|E_{\parallel}|$ = ' + f'{round(Amplitude*1000,3)} [mV/m]\n' +
                           f'$\omega$/k = {WaveSpeed_Static/m_to_km} km/s')
        ax.set_xlabel('V$_{\perp}$ [km/s]', fontsize=14)
        ax.set_ylabel('V$_{\parallel}$ [km/s]', fontsize=14)
        ax.set_xlim(-5000, VelLimit)
        ax.set_ylim(-VelLimit, VelLimit)


        # plot the data
        cbarLow, cbarHigh = 1E6, 5E8
        cmapPlot = ax.pcolormesh(VperpGrid, VparaGrid, diffEFlux[0], cmap='turbo',norm='log', shading='nearest', vmin=cbarLow, vmax=cbarHigh,edgecolor='k',linewidth=0.1)
        ax.invert_yaxis()  # flip the x-axis

        # create the grid
        # ax.pcolormesh(Vperp, Vpara, diffEFlux[0],facecolor='none',edgecolor='k',linewidth=0.1,alpha=1, cmap='turbo', shading='nearest', vmin=0, vmax=1) # creates the grid

        # colorbar
        cbar = fig.colorbar(mappable = cmapPlot)
        cbar.set_label('Differential Energy Flux [eV/cm$^{2}$-s-sr-eV]', fontsize=16)

        # --- ANIMATION FUNCTION ---
        def animate_func(j):
            print('', end='\r' + color.RED + f'{round(100 * j / simLen, 1)} %' + color.END)

            # update the title
            title.set_text(f'$\Delta t$= {round(simTime[j],2)}\n'+
                           r'$|E_{\parallel}|$ = ' + f'{round(Amplitude*1000,3)} [mV/m]\n' +
                           f'$\omega$/k = {WaveSpeed_Static/m_to_km} km/s')


            # update the data
            cmapPlot.set_array(diffEFlux[j])
            # cmapPlot.set_array(diffEFlux[j]/diffEFlux[j].max())

        anim = animation.FuncAnimation(fig=fig, func=animate_func, frames=simLen, interval=1000 / fps_ptichAnglePlot)
        print('', end='\r' + color.RED + f'{round(100 * simLen / simLen, 1)} %' + color.END)
        anim.save(f'C:\Data\ACESII\science\simulations\TestParticle\EEPAA_Polar_Animation.mp4', fps=fps)

        Done(start_time)





# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

testParticle_Sim()