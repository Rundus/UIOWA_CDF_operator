# --- imports ---
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.simToggles import GenToggles,BgeoToggles,m_to_km, R_REF
from numpy import array,degrees,arccos
from numpy.linalg import norm
from datetime import datetime
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from ACESII_code.class_var_func import CHAOS


########################################
# --- GENERATE THE B-FIELD & TOGGLES ---
########################################
plot_BField = False

# --- OUTPUT DATA ------
outputData = True


def generateGeomagneticField(outputData, **kwargs):
    plotBool = kwargs.get('showPlot', False)

    def geomagneticFieldProfile(altRange, plotBool):

        geomagAlts = [((alt + R_REF) / R_REF) for alt in altRange]
        geomagLats = array([degrees(arccos(radi / BgeoToggles.Lshell)) for radi in geomagAlts])
        geomagLongs = array([111.83 for i in range(len(altRange))])
        times = [datetime(2022, 11, 20, 17, 20, 00, 000) for i in range(len(altRange))]
        Pos = array([geomagAlts, geomagLats, geomagLongs]).transpose()
        ISOtime = [times[i].isoformat() for i in range(len(times))]
        cvals_MAG = coord.Coords(Pos, 'MAG', 'sph')
        cvals_MAG.ticks = Ticktock(ISOtime, 'ISO')
        cvals_GDZ = cvals_MAG.convert('GEO', 'sph')
        Lat_geo = cvals_GDZ.lati

        # Get the Chaos model
        B = CHAOS(Lat_geo, [15.25 for i in range(len(altRange))], array(altRange) / m_to_km, times)
        Bgeo = (1E-9) * array([norm(Bvec) for Bvec in B])
        Bgrad = [(Bgeo[i + 1] - Bgeo[i]) / (altRange[i + 1] - altRange[i]) for i in range(len(Bgeo) - 1)]
        Bgrad = array(Bgrad + [Bgrad[-1]]) # add the high altitude value Bgrad again to model the starting point (MAYBE it SHOULD BE 0?)


        if plotBool:

            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(2)
            ax[0].plot(altRange/R_REF, Bgeo/(1E-9))
            ax[0].set_title('|B| vs Altitude')
            ax[0].set_ylabel('$B_{geo}$ [nT]')
            ax[0].set_xlabel('Altitude [$R_{E}$]')
            ax[0].axvline(x=400000/R_REF,label='Observation Height',color='red')
            ax[0].legend()

            ax[1].plot(altRange / R_REF, Bgrad/(1E-9))
            ax[1].set_title(r'$\nabla B$ vs Altitude')
            ax[1].set_ylabel(r'$\nabla B$ [nT/m]')
            ax[1].set_xlabel('Altitude [$R_{E}$]')
            ax[1].axvline(x=400000 / R_REF, label='Observation Height', color='red')
            ax[1].legend()
            plt.tight_layout()
            plt.show()

        return Bgeo, Bgrad


    # get all the variables and plot them if required
    Bgeo, Bgrad = geomagneticFieldProfile(altRange=GenToggles.simAlt, plotBool=plotBool)

    if outputData:

        from copy import deepcopy
        from ACESII_code.class_var_func import outputCDFdata

        # --- Construct the Data Dict ---
        exampleVar = {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808,
                      'FORMAT': 'I5', 'UNITS': 'm', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                      'SCALETYP': 'linear', 'LABLAXIS': 'simAlt'}

        data_dict = {'Bgeo': [Bgeo, {'DEPEND_0':'simAlt', 'UNITS':'T', 'LABLAXIS': 'Bgeo'}],
                      'Bgrad': [Bgrad, {'DEPEND_0':'simAlt', 'UNITS':'T', 'LABLAXIS': 'Bgrad'}],
                      'simAlt': [GenToggles.simAlt, {'DEPEND_0':'simAlt', 'UNITS':'m', 'LABLAXIS': 'simAlt'}]}

        # update the data dict attrs
        for key, val in data_dict.items():
            newAttrs = deepcopy(exampleVar)

            for subKey, subVal in data_dict[key][1].items():
                newAttrs[subKey] = subVal

            data_dict[key][1] = newAttrs

        outputPath = rf'{GenToggles.simOutputPth}\geomagneticField\geomagneticfield.cdf'
        outputCDFdata(outputPath, data_dict)



#################
# --- EXECUTE ---
#################
generateGeomagneticField(outputData=outputData,showPlot=plot_BField)