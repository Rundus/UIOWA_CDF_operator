# --- imports ---
from simulation_Toggles import m_to_km, R_REF, simAlt, simTime
from scipy.stats import linregress
from ACESII_code.class_var_func import lightSpeed, u0, IonMasses,q0,m_e,ep0,cm_to_m
from numpy import exp, sqrt, array,abs

##############################
# --- WAVE E-Field Toggles ---
##############################
Z0_wave = (1.3*6371+3000)*m_to_km # initial altitude of the wave (in meters)

# static
Amplitude_static = 0.18/1000  # V/m
lamda_static = (7500 - 4500)*m_to_km  # in meters
WaveSpeed_Static = 20000*m_to_km  # in meters/s
duty = 40 # duty cycle percent (Max is 50)

# dynamic
Eperp = 12E-3 # V/m
alpha = 1/10 # defines the ratio of the wave's wavelength to the total simulated altitude size
kperp_obs = 1/alpha

flipEField = True
dynamicWaveSpeed = False


def generateAlfvenSpeed(Bmag):
# --- determine the density over all altitudes ---
    # Description: returns density for altitude "z [km]" in m^-3
    h = 0.06 * (R_REF / m_to_km)  # in km from E's surface
    n0 = 6E4
    n1 = 1.34E7
    z0 = 0.05 * (R_REF / m_to_km)  # in km from E's surface
    n = [(cm_to_m ** 3) * (n0 * exp(-1 * ((alt / m_to_km) - z0) / h) + n1 * ((alt / m_to_km) ** (-1.55))) for alt in simAlt]  # calculated density (in cm^-3)

    # --- determine the electron plasma density and skin depth ---
    plasmaFreq = [sqrt((n[j] * q0 * q0) / (m_e * ep0)) for j in range(len(n))]
    skinDepth = [lightSpeed / freq for freq in plasmaFreq]

    # --- determine MHD alfven speed over all altitudes ---
    VA_MHD = array([Bmag[j] / sqrt(u0 * IonMasses[0] * n[j]) for j in range(len(simAlt))])

    # --- determine the Wave Speed at all times ---
    # Note: Use the Speed of the wave at the CENTER of the wave (lambda/2), not at the edges
    WaveCenterPosition = []
    WaveSpeed = []

    for i, tme in enumerate(simTime):

        # get the altitude of the center of the wave
        if i ==0:
            WaveSpeed.append(VA_MHD[abs(simAlt - Z0_wave).argmin()])
            WaveCenterPosition.append(Z0_wave - lamda_static/2)
        else:
            Z_center =  WaveCenterPosition[i-1] - WaveSpeed[i-1]*tme
            WaveCenterPosition.append(Z_center)
            WaveSpeed.append(VA_MHD[abs(simAlt - Z_center).argmin()])

    return WaveSpeed, VA_MHD, WaveCenterPosition



#
# import matplotlib.pyplot as plt
# from geomagnetic_Field_Toggles import generateBField
# Bmag, Bgrad = generateBField()
# WaveSpeed, VA_MHD, WavePos = generateAlfvenSpeed(Bmag)
#
# fig, ax = plt.subplots(3)
# ax[0].set_ylabel('Wave Speed [m/s]')
# ax[0].scatter(simTime, WaveSpeed)
# ax[0].set_yscale('log')
# ax[0].set_ylim(1E6, 1E8)
# ax[0].set_ylabel('Wave Position [m/s]')
# ax[1].scatter(simTime,WavePos)
# ax[2].set_ylabel('Alfven Speed [m/s]')
# ax[2].set_xlabel('Altitude [Re]')
# ax[2].plot(simAlt/R_REF, VA_MHD)
# plt.show()


def generateEField(Bmag):

    # Define the Wave Parameters
    E_Field_WaveSpeed = generateAlfvenSpeed(Bmag) if dynamicWaveSpeed else [WaveSpeed_Static for i in range(len(simTime))]
    sign = -1 if flipEField else 1  # flip the sign of the amplitude
    E_Field_lambda = lamda_static if dynamicWaveSpeed else lamda_static
    E_Field_Amplitude = Amplitude_static if dynamicWaveSpeed else Amplitude_static

    # prepare data for output
    E_Field = []
    E_Field_Position = []

    # Loop through all times to generate an E-Field list for all altitudes at that time
    for k, tme in enumerate(simTime):

        # make sure we have a reasonable duty
        if duty > 50:
            raise Exception('Duty cannot be more than 50%')

        # prepare the E-Field for storage
        E_temp = []

        # determine how much the wave has moved THIS time tick
        deltaAlt = tme * E_Field_WaveSpeed[k] if dynamicWaveSpeed else tme * WaveSpeed_Static  # calc the change in position for this tick
        deltaPos = E_Field_Position[k] - deltaAlt if dynamicWaveSpeed else E_Field_WaveSpeed[k]*tme  # calc the total distance moved thus far

        # Determine the wave amplitude throughout all altitudes at this time tick
        for j, alt in enumerate(simAlt):

            # determine the wave parameters
            riseTime = E_Field_lambda * (50 - duty) / (2 * 100)
            dutyTime = E_Field_lambda * (duty / 100)

            if alt > Z0_wave - deltaPos:  # if we're above the wave
                E = 0
            elif (Z0_wave - deltaPos - riseTime) <= alt < (Z0_wave - deltaPos):  # if we're in the first risetime
                slope, intercept, r, p, se = linregress(x=[Z0_wave - deltaPos, Z0_wave - deltaPos - riseTime],
                                                        y=[0, E_Field_Amplitude])
                E = sign * (slope * alt + intercept)
            elif (Z0_wave - deltaPos - (dutyTime + riseTime)) <= alt < (
                    Z0_wave - deltaPos - riseTime):  # if we're in the first max amplitude part
                E = sign * E_Field_Amplitude
            elif (Z0_wave - deltaPos - (dutyTime + 3 * riseTime)) <= alt < (
                    Z0_wave - deltaPos - (dutyTime + riseTime)):  # if we're in the middle transition region
                slope, intercept, r, p, se = linregress(
                    x=[Z0_wave - deltaPos - (dutyTime + riseTime), Z0_wave - deltaPos - (dutyTime + 3 * riseTime)],
                    y=[E_Field_Amplitude, -1 * E_Field_Amplitude])
                E = sign * (slope * alt + intercept)
            elif (Z0_wave - deltaPos - (2 * dutyTime + 3 * riseTime)) <= alt < (
                    Z0_wave - deltaPos - (dutyTime + 3 * riseTime)):  # if we're in the lowest max amplitude part
                E = -1 * sign * E_Field_Amplitude
            elif (Z0_wave - deltaPos - E_Field_lambda) <= alt < (
                    Z0_wave - deltaPos - (2 * dutyTime + 3 * riseTime)):  # if we're in the last risetime
                slope, intercept, r, p, se = linregress(
                    x=[Z0_wave - deltaPos - (2 * dutyTime + 3 * riseTime), Z0_wave - deltaPos - E_Field_lambda],
                    y=[-1 * E_Field_Amplitude, 0])
                E = sign * (slope * alt + intercept)
            elif alt < (Z0_wave - deltaPos - E_Field_lambda):  # if we're outside the wave on the bottom edge
                E = 0
            else:
                E = 0

            E_temp.append(E)

        E_Field.append(E_temp)

    return E_Field, E_Field_Amplitude, E_Field_WaveSpeed,E_Field_lambda