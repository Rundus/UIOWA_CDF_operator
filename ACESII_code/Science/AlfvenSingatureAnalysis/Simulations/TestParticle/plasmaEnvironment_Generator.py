# --- imports ---
from simToggles import m_to_km, R_REF, GenToggles, EToggles
from ACESII_code.class_var_func import lightSpeed, u0,q0,m_e,ep0,cm_to_m, IonMasses
from numpy import exp, sqrt, array, pi, abs, tanh
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.geomagneticField_Generator import geomagneticFieldProfile



##################
# --- PLOTTING ---
##################
plot_Temperature = False
plot_Density = False
plot_ionMass = True
plot_PlasmaFreq = False
plot_PlasmaBeta = False
plot_skinDepth = False
plot_ionCyclotron = False
plot_ionLarmorRadius = False
plot_MHDalfvenSpeed = False
plot_lambdaPerp = False
plot_lambdaPara = False
plot_kineticAlfSpeed = False

simulationAlt = GenToggles.simAlt


# --- Temperature ---
def temperatureProfile(altRange, **kwargs):
    plotBool = kwargs.get('showPlot', False)

    # --- Ionosphere Temperature Profile ---
    # ASSUMES Ions and electrons have same temperature profile
    T0 = 1 # Temperature at the Ionospher (in eV)
    T1 = 0.0135 # (in eV)
    h0 = 2000*m_to_km # scale height (in meters)
    T_iono = T1*exp(altRange/h0 ) + T0
    deltaZ = 0.3*R_REF
    T_ps = 2000 # temperature of plasma sheet (in eV)
    z_ps = 3.75*R_REF # height of plasma sheet (in meters)
    w = 0.5*(1 - tanh((altRange - z_ps)/deltaZ)) # models the transition to the plasma sheet

    # determine the overall temperature profile
    T_e = [T_iono[i]*w[i] + T_ps*(1 - w[i]) for i in range(len(altRange))]

    if plotBool:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3)
        fig.set_figwidth(15)
        fig.set_figheight(10)

        ax[0].plot(altRange / R_REF, T_iono)
        ax[0].set_title('Ionospheric Temperature Profile vs Altitude')
        ax[0].set_ylabel('Temperature [eV]')
        ax[0].set_xlabel('Altitude [$R_{E}$]')
        ax[0].set_yscale('log')
        ax[0].axvline(x=400000 / R_REF, label='Observation Height', color='red')

        ax[1].plot(altRange / R_REF, w)
        ax[1].set_title('Weighting Function vs Altitude')
        ax[1].set_ylabel('Weighting Function')
        ax[1].set_xlabel('Altitude [$R_{E}$]')
        ax[1].axvline(x=400000 / R_REF, label='Observation Height', color='red')

        ax[2].plot(altRange / R_REF, T_e)
        ax[2].set_yscale('log')
        ax[2].set_title('Total Electron Temperature vs Altitude')
        ax[2].set_ylabel('Electron Temperature [eV]')
        ax[2].set_xlabel('Altitude [$R_{E}$]')
        ax[2].axvline(x=400000 / R_REF, label='Observation Height', color='red')
        plt.legend()
        plt.tight_layout()
        plt.show()

    return T_e



# --- Kperp ---

def lambdaPerpProfile(altRange, **kwargs):
    plotBool = kwargs.get('showPlot', False)

    Bmag,Bgrad = geomagneticFieldProfile(altRange)
    initindex = abs(altRange - EToggles.Z0_wave).argmin() # the index of the startpoint of the Wave
    initBmag = Bmag[initindex] # <--- This determines where the scaling begins
    LambdaPerp = array([EToggles.lambdaPerp0*sqrt(initBmag/Bmag[i]) for i in range(len(altRange))]) if not EToggles.static_Kperp else array([EToggles.lambdaPerp0 for i in range(len(altRange))])
    kperp = array([2*pi/LambdaPerp[i] for i in range(len(altRange))])

    if plotBool:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2)
        ax[0].plot(altRange / R_REF, LambdaPerp)
        ax[0].set_title('$\lambda_{\perp}$ vs Altitude \n'
                        '$\lambda_{\perp 0}$=' +f'{EToggles.lambdaPerp0}')
        ax[0].set_ylabel('Perpendicular Wavelength [m]')
        ax[0].set_xlabel('Altitude [$R_{E}$]')
        ax[0].axvline(x=400000 / R_REF, label='Observation Height', color='red')

        ax[1].plot(altRange / R_REF, kperp)
        ax[1].set_title('$k_{\perp}$ vs Altitude')
        ax[1].set_ylabel('Perpendicular Wavenumber [m$^{-1}$]')
        ax[1].set_xlabel('Altitude [$R_{E}$]')
        ax[1].axvline(x=400000 / R_REF, label='Observation Height', color='red')
        plt.legend()
        plt.tight_layout()
        plt.show()

    return LambdaPerp, kperp



# --- PLASMA DENSITY ---
# uses the Klezting Model to return an array of plasma density (in m^-3) from [Alt_low, ..., Alt_High]
def plasmaDensityProfile(altRange, **kwargs):
    plotBool = kwargs.get('showPlot', False)

    # --- determine the density over all altitudes ---
    # Description: returns density for altitude "z [km]" in m^-3
    h = 0.06 * (R_REF / m_to_km)  # in km from E's surface
    n0 = 6E4
    n1 = 1.34E7
    z0 = 0.05 * (R_REF / m_to_km)  # in km from E's surface
    n_density = (cm_to_m**3)*array([(n0 * exp(-1 * ((alt / m_to_km) - z0) / h) + n1 * ((alt / m_to_km) ** (-1.55))) for alt in altRange])  # calculated density (in m^-3)

    if plotBool:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(altRange/R_REF, n_density/(cm_to_m**3))
        ax.set_title('$n_{i}$ vs Altitude')
        ax.set_ylabel('Plasma Density [$cm^{-3}$]')
        ax.set_xlabel('Altitude [$R_{E}$]')
        ax.axvline(x=400000/R_REF,label='Observation Height',color='red')
        plt.legend()
        plt.show()

    return n_density


# --- Ion Mass ---
def ionMassProfile(altRange, **kwargs):
    plotBool = kwargs.get('showPlot', False)

    plasmaDensity_total = plasmaDensityProfile(altRange)
    z_i = 2370*m_to_km  #
    h_i = 1800*m_to_km  # height of plasma sheet (in meters)
    n_Op = array([plasmaDensity_total[i]*0.5 * (1 - tanh((altRange[i] - z_i) / h_i)) for i in range(len(altRange))])
    n_Hp = plasmaDensity_total - n_Op
    m_0p = IonMasses[1]
    m_Hp = IonMasses[2]
    m_eff_i = array([ m_Hp*0.5*(1 + tanh( (altRange[i] - z_i)/h_i )) + m_0p*0.5*(1 - tanh( (altRange[i] - z_i)/h_i )) for i in range(len(altRange))])


    if plotBool:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3)
        ax[0].plot(altRange / R_REF, n_Op)
        ax[0].set_title('$n_{0^{+}}$ vs Altitude')
        ax[0].set_ylabel('Monatomic-Oxygen number density [cm^{3}]')
        ax[0].set_xlabel('Altitude [$R_{E}$]')
        ax[0].axvline(x=400000 / R_REF, label='Observation Height', color='red')


        ax[1].plot(altRange / R_REF, n_Hp)
        ax[1].set_title('$n_{H^{+}}$ vs Altitude')
        ax[1].set_ylabel('Monatomic-Hydrogen number density [cm^{3}]')
        ax[1].set_xlabel('Altitude [$R_{E}$]')
        ax[1].axvline(x=400000 / R_REF, label='Observation Height', color='red')

        ax[2].plot(altRange / R_REF, m_eff_i)
        ax[2].set_title('$m_{eff_{i}}$ vs Altitude')
        ax[2].set_ylabel('Effective Ion Mass [kg]')
        ax[2].set_xlabel('Altitude [$R_{E}$]')
        ax[2].axvline(x=400000 / R_REF, label='Observation Height', color='red')
        plt.legend()
        plt.tight_layout()
        plt.show()

    return n_Op,n_Hp, m_eff_i


# --- PLASMA BETA ---
def plasmaBetaProfile(altRange, **kwargs):
    plotBool = kwargs.get('showPlot', False)

    plasmaDensity = plasmaDensityProfile(altRange)
    Bgeo, Bgrad = geomagneticFieldProfile(altRange)
    Te = 50
    plasmaBeta = array([(plasmaDensity[i]*q0*Te)/(Bgeo[i]**2 /(2*u0)) for i in range(len(altRange))])
    IonMass = 5E-26
    ratio = m_e/IonMass

    if plotBool:
        import matplotlib.pyplot as plt

        Temps = [1, 10, 50]
        colors= ['tab:red','tab:blue','tab:green']

        fig, ax = plt.subplots()

        for k,temp in enumerate(Temps):
            plasmaBeta = array([(plasmaDensity[i] * q0 * temp) / (Bgeo[i] ** 2 / (2 * u0)) for i in range(len(altRange))])
            ax.plot(altRange/R_REF, plasmaBeta/ratio, color=colors[k], label=f'T_e = {temp} eV')

        ax.set_title(r'$\beta$ vs Altitude')
        ax.set_ylabel('Plasma Beta / (m_e/m_i)')
        ax.set_xlabel('Altitude [$R_{E}$]')
        ax.axvline(x=400000/R_REF, label='Observation Height', color='red',linestyle='--')
        ax.axhline(y=1, color='black')
        ax.set_yscale('log')
        ax.grid(True)
        plt.legend()
        plt.show()

    return plasmaBeta


# --- PLASMA FREQ ---
def plasmaFreqProfile(altRange,**kwargs):
    plotBool = kwargs.get('showPlot', False)

    plasmaDensity = plasmaDensityProfile(altRange)
    plasmaFreq = array([ sqrt( plasmaDensity[i]* (q0*q0) / (ep0*m_e)) for i in range(len(plasmaDensity))])

    if plotBool:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(altRange/R_REF, plasmaFreq)
        ax.set_title('$\omega_{pe}$ vs Altitude')
        ax.set_ylabel('Plasma Freq [rad/s]')
        ax.set_xlabel('Altitude [$R_{E}$]')
        ax.axvline(x=400000/R_REF,label='Observation Height',color='red')
        plt.legend()
        plt.show()

    return plasmaFreq

# --- SKIN DEPTH ---
def skinDepthProfile(altRange,**kwargs):
    plotBool = kwargs.get('showPlot', False)

    plasmaFreq = plasmaFreqProfile(altRange)
    skinDepth = array([lightSpeed/plasmaFreq[i] for i in range(len(plasmaFreq))])

    if plotBool:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(altRange/R_REF, skinDepth)
        ax.set_title('$\lambda_{e}$ vs Altitude')
        ax.set_ylabel('Skin Depth [m]')
        ax.set_xlabel('Altitude [$R_{E}$]')
        ax.axvline(x=400000/R_REF,label='Observation Height',color='red')
        plt.legend()
        plt.show()
    return skinDepth

# --- ION CYCLOTRON FREQ ---
def ionCyclotronProfile(altRange,**kwargs):
    plotBool = kwargs.get('showPlot', False)

    Bgeo,Bgrad = geomagneticFieldProfile(altRange)
    IonMass = 5E-26

    ionCyclotron = array([q0*Bgeo[i]/IonMass for i in range(len(Bgeo))])

    if plotBool:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2)
        ax[0].plot(altRange/R_REF, ionCyclotron)
        ax[0].set_title('$\omega_{ci}$ vs Altitude')
        ax[0].set_ylabel('Ion Cyclotron [m]')
        ax[0].set_xlabel('Altitude [$R_{E}$]')
        ax[0].axvline(x=400000/R_REF,label='Observation Height',color='red')

        ax[1].plot(altRange / R_REF, ionCyclotron / (2*pi))
        ax[1].set_title('$\Omega_{ci}$ vs Altitude')
        ax[1].set_ylabel('Ion Cyclotron [Hz]')
        ax[1].set_xlabel('Altitude [$R_{E}$]')
        ax[1].axvline(x=400000 / R_REF, label='Observation Height', color='red')
        plt.legend()
        plt.tight_layout()
        plt.show()

    return ionCyclotron


# --- Ion Larmor Radius ---
def ionLarmorRadiusProfile(altRange,**kwargs):
    plotBool = kwargs.get('showPlot', False)

    ionCyclo = ionCyclotronProfile(altRange)
    IonMass = 5E-26
    vth = sqrt(8*q0*Ti/IonMass)
    ionLarmorRadius = array([vth/ionCyclo[i] for i in range(len(ionCyclo))])

    if plotBool:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(altRange / R_REF, ionLarmorRadius)
        ax.set_title(r'$\rho_{i}$ vs Altitude')
        ax.set_ylabel('Ion Larmor Radius [rad/s]')
        ax.set_xlabel('Altitude [$R_{E}$]')
        ax.axvline(x=400000 / R_REF, label='Observation Height', color='red',linestyle='--')
        ax.axvline(x=10000000 / R_REF, label='Magnetosheath Proton Limit', color='tab:green',linestyle='--')
        ax.set_yscale('log')
        ax.set_ylim(0,1E4)
        ax.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

    return ionLarmorRadius


# --- MHD Alfven Speed ---
def MHD_alfvenSpeedProfile(altRange,**kwargs):
    plotBool = kwargs.get('showPlot', False)

    plasmaDensity = plasmaDensityProfile(altRange)
    Bmag,Bgrad = geomagneticFieldProfile(altRange)
    IonMass_Avg = 3*1.67E-27
    VA_MHD = array([ Bmag[i]/sqrt(u0*IonMass_Avg*plasmaDensity[i]) for i in range(len(plasmaDensity))])

    if plotBool:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(altRange / R_REF, VA_MHD/(m_to_km))
        ax.set_title(r'$V_{A}$ vs Altitude')
        ax.set_ylabel('MHD Alfven Speed  [km/s]')
        ax.set_xlabel('Altitude [$R_{E}$]')
        ax.axvline(x=400000 / R_REF, label='Observation Height', color='red', linestyle='--')
        ax.grid(True)

        # plot some thermal velocity comparisons
        Vth_low = sqrt(8*q0*1/(9.11E-31))/m_to_km
        ax.axhline(y=Vth_low)
        ax.text(x=R_REF/R_REF,y=Vth_low*1.3,s='$V_{th_{e}}$ (1 eV)', color='black')

        Vth_high = sqrt(8 * q0 * 50 / (9.11E-31))/m_to_km
        ax.axhline(y=Vth_high)
        ax.text(x=R_REF/R_REF, y=Vth_high * 1.1, s='$V_{th_{e}}$ (50 eV)', color='black')

        plt.legend()
        plt.tight_layout()
        plt.show()

    return VA_MHD


# --- Lambda Parallel Wavelength ---
def lambdaParallelProfile(altRange, **kwargs):
    plotBool = kwargs.get('showPlot', False)

    # collect profiles
    lambdaPerp, kperp = lambdaPerpProfile(altRange)
    ionCyclo = ionCyclotronProfile(altRange)
    ionLarmorRadi = ionLarmorRadiusProfile(altRange)
    skinDepth = skinDepthProfile(altRange)
    alfSpdMHD = MHD_alfvenSpeedProfile(altRange)

    LambdaPara = array([2*pi*alfSpdMHD[i]*sqrt(1 - (EToggles.waveFreqionCyclo[i])**2)*sqrt(1 + (kperp[i]*ionLarmorRadi[i])**2)/(EToggles.waveFreqsqrt(1 + (kperp[i]*skinDepth[i])**2)) for i in range(len(altRange))])
    kpara = array([2*pi/LambdaPara[i] for i in range(len(altRange))])

    if plotBool:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2)
        ax[0].plot(altRange / R_REF, LambdaPara/m_to_km)
        ax[0].set_title('$\lambda_{\parallel}$ vs Altitude')
        ax[0].set_ylabel('Parallel Wavelength [km]')
        ax[0].set_xlabel('Altitude [$R_{E}$]')
        ax[0].axvline(x=400000 / R_REF, label='Observation Height', color='red')

        ax[1].plot(altRange / R_REF, kpara)
        ax[1].set_title('$k_{\parallel}$ vs Altitude')
        ax[1].set_ylabel('Parallel Wavenumber [m$^{-1}$]')
        ax[1].set_xlabel('Altitude [$R_{E}$]')
        ax[1].axvline(x=400000 / R_REF, label='Observation Height', color='red')
        plt.legend()
        plt.tight_layout()
        plt.show()

    return LambdaPara,kpara


# --- Kinetic Alfven Speed ---
def kinetic_alfvenSpeedProfile(altRange, **kwargs):
    plotBool = kwargs.get('showPlot', False)

    # collect profiles
    lambdaPerp, kperp = lambdaPerpProfile(altRange)
    ionCyclo = ionCyclotronProfile(altRange)
    ionLarmorRadi = ionLarmorRadiusProfile(altRange)
    skinDepth = skinDepthProfile(altRange)
    alfSpdMHD = MHD_alfvenSpeedProfile(altRange)
    kineticAlfSpeed = array([ alfSpdMHD[i] * sqrt(1 - (EToggles.waveFreq/ionCyclo[i])**2) *sqrt(1 + (kperp[i]*ionLarmorRadi[i])**2)/(sqrt(1 + (kperp[i]*skinDepth[i])**2)) for i in range(len(altRange))])

    if plotBool:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(altRange / R_REF, kineticAlfSpeed/(m_to_km), label='kinetic Alf speed', color='blue')
        ax.set_title(r'$\omega_{wave}/k_{\parallel}$ vs Altitude' +
                     '\n' + r'$\omega_{wave0}$=' + f'{EToggles.waveFreq} rad/s')
        ax.set_ylabel('Kinetic Alfven Speed  [km/s]')
        ax.set_xlabel('Altitude [$R_{E}$]')
        ax.axvline(x=400000 / R_REF, label='Observation Height', color='red', linestyle='--')
        ax.grid(True)

        # plot the MHD alfven speed
        ax.plot(altRange / R_REF, alfSpdMHD/(m_to_km), color='red')

        # plot some thermal velocity comparisons
        Vth_low = sqrt(8*q0*1/(9.11E-31))/m_to_km
        ax.axhline(y=Vth_low)
        ax.text(x=R_REF/R_REF,y=Vth_low*1.3,s='$V_{th_{e}}$ (1 eV)', color='black')

        Vth_high = sqrt(8 * q0 * 50 / (9.11E-31))/m_to_km
        ax.axhline(y=Vth_high)
        ax.text(x=R_REF/R_REF, y=Vth_high * 1.1, s='$V_{th_{e}}$ (50 eV)', color='black')

        plt.legend()
        plt.tight_layout()
        plt.show()

    return kineticAlfSpeed


#################
# --- EXECUTE ---
#################
if plot_Temperature:
    temperatureProfile(simulationAlt, showPlot=plot_Temperature)

if plot_Density:
    plasmaDensityProfile(simulationAlt, showPlot=plot_Density)

if plot_ionMass:
    ionMassProfile(simulationAlt, showPlot=plot_ionMass)

if plot_PlasmaBeta:
    plasmaBetaProfile(simulationAlt, showPlot=plot_PlasmaBeta)

if plot_PlasmaFreq:
    plasmaFreqProfile(simulationAlt, showPlot=plot_PlasmaFreq)

if plot_skinDepth:
    skinDepthProfile(simulationAlt, showPlot=plot_skinDepth)

if plot_ionCyclotron:
    ionCyclotronProfile(simulationAlt,showPlot=plot_ionCyclotron)

if plot_ionLarmorRadius:
    ionLarmorRadiusProfile(simulationAlt,showPlot=plot_ionLarmorRadius)

if plot_MHDalfvenSpeed:
    MHD_alfvenSpeedProfile(simulationAlt,showPlot=plot_MHDalfvenSpeed)

if plot_kineticAlfSpeed:
    kinetic_alfvenSpeedProfile(simulationAlt, showPlot=plot_kineticAlfSpeed)
