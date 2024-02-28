# --- imports ---
import numpy as np
from simToggles import m_to_km, R_REF, GenToggles,EToggles
from scipy.stats import linregress
from ACESII_code.class_var_func import lightSpeed, u0, IonMasses,q0,m_e,ep0,cm_to_m
from numpy import exp, sqrt, array,abs
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.geomagneticField_Generator import geomagneticFieldProfile



# TODO:
#  (1) Get Average Ion Mass vs Altitude up to 2Re to fix: Cyclotron Freq, Larmor Radius, MHD_AlfvenSpeed
#  (2) Implement Temperature vs Altitude profile for Ions for Larmor Radius. Note: Tanaka 2005 states <10,000 km Ionosphere is cold (<1eV) but > 10,000km magnetosheath protons are imporant (400eV)

##############################
# --- WAVE E-Field Profile ---
##############################

# Generates the entire E-Field for all time in the simulation and returns a variable that looks like:
# [
#   [Ez(x=0,t=0),Ez(x=1,t=0),Ez(x=2,t=0)...., Ez(x=len(Alt),t=0)],
#   [Ez(x=0,t=1),Ez(x=1,t=1),Ez(x=2,t=1)...., Ez(x=len(Alt),t=1)]
#   ,...]

def alfvenWaveEFieldProfile(Bmag):

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


#################
# --- EXECUTE ---
#################
