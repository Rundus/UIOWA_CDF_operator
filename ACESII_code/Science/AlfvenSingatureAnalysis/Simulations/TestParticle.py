# --- TestParticle.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:
import numpy as np

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.myImports import *
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from scipy.stats import linregress
from ACESII_code.class_var_func import CHAOS
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import matplotlib.patches as mpatches


start_time = time.time()



colors = ['tab:purple','tab:orange','tab:red','tab:green','tab:blue']
m_to_km = 1E3

# --- --- --- --- --- --- --
# --- SIMULATION TOGGLES ---
# --- --- --- --- --- --- --
# simToggles
simLen = 235 # how many delta T steps to simulate
deltaT = 0.01  # in seconds
simAltLow = 1000*m_to_km # low altitude in meters
simAltHigh = 8000*m_to_km # high altitude in meters

# particle toggles
Z0 = 6500*m_to_km # initial altitude of particles (in meters)
N = 64 # number of particles in each energy
ptcl_mass = m_e # mass of the particle
ptcl_charge = q0 # charge of the particle
Energies = [10, 20, 30, 40, 50] # choose some energies

# E-Field Toggles
Z0_wave = 7500*m_to_km # initial altitude of the wave (in meters)
Amplitude = 0.0012/4  # V/m
lamda = (7500 - 3500)*m_to_km  # in meters
WaveSpeed = 2000*m_to_km  # in meters/s
duty = 49 # duty cycle percent (Max is 50)

# B-Field Toggles
Lshell = 9

# observation togles
obsHeight = 1000*m_to_km # in meters

# output


# --- toggles ---
inverse = False
plotEField = False
plotMirror_vs_electrostatic = False
plotInputEnergies = False
plotEuler_1step = False
# -------------------------
outputAnimation = True
fps = 30
# -------------------------

outputTitle = 'TestParticle_above_invTrue.mp4' if inverse else 'TestParticle_above_invFalse.mp4'
def testParticle_Sim():

    def E_Field(simTime, Alt):

        if duty > 50:
            raise Exception('Duty cannot be more than 50%')

        sign = -1 if inverse else 1
        deltaPos = simTime * WaveSpeed
        riseTime = lamda*(50 - duty)/(2*100)
        dutyTime = lamda*(duty/100)

        if Alt > Z0_wave - deltaPos:  # if we're above the wave
            E = 0
        elif (Z0_wave - deltaPos - riseTime) <= Alt < (Z0_wave - deltaPos):  # if we're in the first risetime
            slope, intercept, r, p, se = linregress(x=[Z0_wave - deltaPos, Z0_wave - deltaPos - riseTime], y=[0, Amplitude])
            E = sign * (slope * Alt + intercept)
        elif (Z0_wave - deltaPos - (dutyTime+riseTime)) <= Alt < (Z0_wave - deltaPos - riseTime):  # if we're in the first max amplitude part
            E = sign * Amplitude
        elif (Z0_wave - deltaPos - (dutyTime+3*riseTime)) <= Alt < (Z0_wave - deltaPos - (dutyTime+riseTime)):  # if we're in the middle transition region
            slope, intercept, r, p, se = linregress(x=[Z0_wave - deltaPos - (dutyTime+riseTime), Z0_wave - deltaPos - (dutyTime+3*riseTime)], y=[Amplitude, -1 * Amplitude])
            E = sign * (slope * Alt + intercept)
        elif (Z0_wave - deltaPos - (2*dutyTime+3*riseTime)) <= Alt < (Z0_wave - deltaPos - (dutyTime+3*riseTime)):  # if we're in the lowest max amplitude part
            E = -1 * sign * Amplitude
        elif (Z0_wave - deltaPos - lamda) <= Alt < (Z0_wave - deltaPos - (2*dutyTime+3*riseTime)):  # if we're in the last risetime
            slope, intercept, r, p, se = linregress(x=[Z0_wave - deltaPos - (2*dutyTime+3*riseTime), Z0_wave - deltaPos - lamda], y=[-1 * Amplitude, 0])
            E = sign * (slope * Alt + intercept)
        elif Alt < (Z0_wave - deltaPos - lamda):  # if we're outside the wave on the bottom edge
            E = 0
        else:
            E = 0

        return E

    def forceFunc(simTime, Alt, mu, deltaB): # gives meters/s^2
        # return (-1 * ptcl_charge / ptcl_mass) * E_Field(simTime, Alt)
        # return  - (mu / ptcl_mass) * deltaB
        return  (-1*ptcl_charge / ptcl_mass) * E_Field(simTime, Alt) - (mu/ptcl_mass)*deltaB

    def AdamsBashforth(yn1, funcVal_n, funcVal_n1):
        yn2 = yn1 + (3 / 2) * deltaT * funcVal_n1 - (1 / 2) * deltaT * funcVal_n
        return yn2

    def Euler(yn, funcVal_n):
        yn1 = yn + deltaT * funcVal_n
        return yn1

    def calcE(Vperp,Vpar,mass,charge):
        return 0.5*(mass/charge)*(Vperp**2 + Vpar**2)


    ##########################
    # --- SIMULATION START ---
    ##########################
    Alt = np.linspace(simAltLow, simAltHigh, 2000)  # in METERS
    simTime = [deltaT * i for i in range(simLen)]

    # --- --- --- --- --- --- --- --- -
    # --- GET THE GEOMAGNETIC FIELD ---
    # --- --- --- --- --- --- --- --- -
    prgMsg('Getting Geomagnetic Field')

    R_REF = 6371.2*m_to_km # in meters
    geomagAlts = [((alt+R_REF)/R_REF) for alt in Alt]
    geomagLats = np.array([np.degrees(np.arccos(radi/Lshell)) for radi in geomagAlts])
    geomagLongs = np.array([111.83 for i in range(len(Alt))])
    times = [dt.datetime(2022, 11, 20, 17, 20, 00, 000) for i in range(len(Alt))]
    Pos = np.array([geomagAlts, geomagLats, geomagLongs]).transpose()
    ISOtime = [times[i].isoformat() for i in range(len(times))]
    cvals_MAG = coord.Coords(Pos, 'MAG', 'sph')
    cvals_MAG.ticks = Ticktock(ISOtime, 'ISO')
    cvals_GDZ = cvals_MAG.convert('GEO','sph')
    Alt_geo = cvals_GDZ.radi
    Lat_geo = cvals_GDZ.lati
    Long_geo = cvals_GDZ.long

    # Get the Chaos model
    B = CHAOS(Lat_geo, [15.25 for i in range(len(Alt))], np.array(Alt)/m_to_km, times)
    Bmag = (1E-9)*np.array([np.linalg.norm(Bvec) for Bvec in B])
    Bgrad = [ (Bmag[i+1] - Bmag[i])/(Alt[i+1] - Alt[i]) for i in range(len(Bmag)-1)] + [0]

    Done(start_time)


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
        pitches = np.radians(np.linspace(0, 360, N))
        data_dict[f'{Energies[i]}_zpos'].append([Z0 for ptch in pitches]) # z position
        data_dict[f'{Energies[i]}_vpar'].append([V_mag * np.cos(ptch) for ptch in pitches]) # parallel velocity
        data_dict[f'{Energies[i]}_vperp'].append([V_mag * np.sin(ptch) for ptch in pitches]) # perpendicular velocity
        data_dict[f'{Energies[i]}_energies'].append([calcE(V_mag * np.sin(ptch), V_mag * np.cos(ptch), ptcl_mass, ptcl_charge) for ptch in pitches]) # particle energy
        data_dict[f'{Energies[i]}_Bmag'].append([Bmag[np.abs(Alt-zpos).argmin()] for zpos in data_dict[f'{Energies[i]}_zpos'][0]]) # determine B-Field
        data_dict[f'{Energies[i]}_Bgrad'].append([Bgrad[np.abs(Alt-zpos).argmin()] for zpos in data_dict[f'{Energies[i]}_zpos'][0]]) # determine deltaB
        data_dict[f'{Energies[i]}_moment'].append([0.5*m_e*(data_dict[f'{Energies[i]}_vperp'][0][k]**2)/data_dict[f'{Energies[i]}_Bmag'][0][k] for k in range(N)])  # determine magnetic moment
        data_dict[f'{Energies[i]}_force'].append([forceFunc(simTime=0, Alt=data_dict[f'{Energies[i]}_zpos'][0][k], mu=data_dict[f'{Energies[i]}_moment'][0][k], deltaB=data_dict[f'{Energies[i]}_Bgrad'][0][k]) for k in range(N)]) # determine the initial force function
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
            newVpar = Euler(yn=vPars0[k], funcVal_n=forceFunc(simTime[0], zPos0[k], deltaB=deltaBs0[k], mu=mus0[k]))
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
            data_dict[f"{engy}_force"][-1].append(forceFunc(simTime=simTime[1], Alt=newPos, deltaB=newdB, mu=newMu))

            # determine if new particle is observed
            data_dict[f"{engy}_observed"][-1].append(1) if newPos <= obsHeight else data_dict[f"{engy}_observed"][-1].append(0)

            # determine new particle energies
            data_dict[f"{engy}_energies"][-1].append(calcE(newVperp,newVpar,ptcl_mass,ptcl_charge))
    Done(start_time)

    # --- IMPLEMENT ADAMS BASHFORTH---
    prgMsg('Adams Bashforth')
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
                    data_dict[f"{engy}_force"][-1].append(forceFunc(simTime=tme, Alt=newPos, deltaB=newdB, mu=mu_n[i]))

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
    TEmaxE = 0
    for tme in range(simLen):
        for engy in Energies:
            observed = np.array(data_dict[f'{engy}_observed'][tme])
            for ptclN in range(N):
                if observed[ptclN] == 1:
                    observed_particle_enery = data_dict[f'{engy}_energies'][tme][ptclN]
                    TEmaxE = observed_particle_enery if observed_particle_enery > TEmaxE else TEmaxE  # find the largest observed energy

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
        wTimes = [simTime[i] for i in range(0,len(simTime),int(len(simTime)/1))]
        for i in range(len(wTimes)):
            ax.plot(Alt,[E_Field(wTimes[i], alt) for alt in Alt],color= f'{colors[i]}',label='$\Delta t = $' + f'{wTimes[i]}')
        fig.legend()
        plt.show()

    if plotMirror_vs_electrostatic:
        chosenPitch = [0, 45, 90]
        fig, ax = plt.subplots(len(chosenPitch))

        # ax.scatter(Bgrad, Alt)
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


    # --- --- --- --- ----
    # --- ANIMATE PLOT ---
    # --- --- --- --- ----

    if outputAnimation:
        prgMsg('Creating Animation')

        #################
        # --- ANIMATE ---
        #################
        fig = plt.figure()
        figure_height = 15
        figure_width = 15
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)

        # --- --- --- --- --
        # --- INITIALIZE ---
        # --- --- --- --- --
        title = fig.suptitle(f'$\Delta t$= {simTime[0]}\n' +
                             r'$E_{\parallel}$ = ' + f'{Amplitude*1000} [mV/m]')

        # gridspec
        gs0 = gridspec.GridSpec(nrows=3, ncols=2, figure=fig)

        # Altitude Plot
        axAlt = fig.add_subplot(gs0[0, 0])
        axAlt.axhline(y=obsHeight/m_to_km, linestyle='--')
        axAlt.set_xticks([-180 + 60*i for i in range(7)])
        axAlt.set_ylim(500, 8000)
        axAlt.set_xlabel('PA $[^{\circ}$]')
        axAlt.set_ylabel('Alt [km]')

        # Show the E-Field
        axE = axAlt.twiny()
        E_Field_Wave, = axE.plot([E_Field(simTime[0], alt) for alt in Alt], Alt/m_to_km, color='black')

        TopArrowX = -1*Amplitude / 2 if inverse else Amplitude/2
        TopArrowY = (Z0_wave - (0.45*lamda))/m_to_km if inverse else (Z0_wave - (0.05*lamda))/m_to_km
        BotArrowX = Amplitude / 2 if inverse else -1*Amplitude/2
        BotArrowY = (Z0_wave - ((1-0.45)*lamda))/m_to_km if inverse else (Z0_wave - ((1-0.05)*lamda))/m_to_km
        dx = 0
        dyTop = 0.35*lamda/m_to_km if inverse else -0.35*lamda/m_to_km
        dyBot = -0.35*lamda/m_to_km if inverse else 0.35*lamda/m_to_km

        topAr = axE.annotate("E$_{\parallel}$", xy=(TopArrowX, TopArrowY), xytext=(TopArrowX + dx, TopArrowY + dyTop),horizontalalignment="center", arrowprops=dict(arrowstyle="->",lw=5),fontsize=17)
        botAr = axE.annotate("E$_{\parallel}$", xy=(BotArrowX, BotArrowY), xytext=(BotArrowX + dx, BotArrowY + dyBot),horizontalalignment="center", arrowprops=dict(arrowstyle="->",lw=5),fontsize=17)

        axE.set_ylim(500, 8000)
        axE.set_xlim(-1.5*Amplitude, 1.5*Amplitude)
        # axE.axis('off')

        # Vspace Plot
        axVspace = fig.add_subplot(gs0[0, 1])
        axVspace.set_ylabel('Vpar [km/s]')
        axVspace.set_xlabel('Vperp [km/s]')
        axVspace.set_xlim(-5E3, 5E3)
        axVspace.set_ylim(-5E3, 5E3)

        # Enegy Plot vs Time
        axTspaceE = fig.add_subplot(gs0[1, 0:2])
        TEplot = axTspaceE.scatter([],[])
        axTspaceE.set_xlim(0, simTime[-1])
        axTspaceE.set_ylim(0, 1.2*TEmaxE) if TEmaxE != 0 else axTspaceE.set_ylim(0, 500)
        axTspaceE.set_ylabel('Energy [eV]')
        axTspaceE.set_xlabel('Time [s]')

        patches = [mpatches.Patch(color=colors[i], label=f'{Energies[i]} eV') for i in range(len(colors))]
        axTspaceE.legend(handles=patches)

        # Time vs Pitch Angle
        axTspacePitch = fig.add_subplot(gs0[2, 0:2])
        TPplot = axTspacePitch.scatter([], [])
        axTspacePitch.set_xlim(0, simTime[-1])
        axTspacePitch.set_yticks([-180 + 60 * i for i in range(7)])
        axTspacePitch.set_ylabel('Pitch Angle [$^\circ$]')
        axTspacePitch.set_xlabel('Time [s]')

        # --- INITIALIZE THE PLOTS ---
        AltitudeArtists = []
        VspaceArtists = []

        for i in range(len(Energies)):

            vperps = np.array(data_dict[f'{Energies[i]}_vperp'][0]) / m_to_km
            vpars = np.array(data_dict[f'{Energies[i]}_vpar'][0]) / m_to_km

            Pitches = np.array([np.degrees(np.arccos(-1*vpars[k] /np.sqrt(vpars[k]**2 + vperps[k]**2))) if vperps[k] > 0 else -1*np.degrees(np.arccos(-1*vpars[k] /np.sqrt(vpars[k]**2 + vperps[k]**2))) for k in range(len(vpars))])
            zpos = np.array(data_dict[f'{Energies[i]}_zpos'][0]) / m_to_km

            # Plot altitudes
            AltitudeArtists.append(axAlt.scatter(Pitches, zpos, color=colors[i]))

            # plot energies
            VspaceArtists.append(axVspace.scatter(vperps, vpars, color=colors[i]))

        # animation function.  This is called sequentially
        ptcls_to_plot = []
        def animate_func(j):
            print('',end='\r' + color.RED + f'{round(100 * j / simLen, 1)} %' + color.END)

            # update the title
            title.set_text(f'$\Delta t$= {simTime[j]}')

            axAlt.set_xticks([-180 + 60 * i for i in range(7)])

            # update the electric field
            E_Field_Wave.set_xdata([E_Field(simTime[j], alt) for alt in Alt])
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


                # update the Energy vs Time Plot
                TEplot.set_offsets(np.stack([plotData[0], plotData[1]]).T)
                TEplot.set_facecolor(plotData[3])


                # update the Pitch Angle vs Time Plot
                TPplot.set_offsets(np.stack([plotData[0], plotData[2]]).T)
                TPplot.set_facecolor(plotData[3])



        anim = animation.FuncAnimation(fig=fig, func=animate_func, frames=[i for i in range(simLen)], interval=1000 / fps)
        anim.save(f'C:\Data\ACESII\science\Movies\{outputTitle}', fps=fps)
        print('\n')
        Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

testParticle_Sim()