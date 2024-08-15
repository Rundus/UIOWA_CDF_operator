# --- imports ---
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.simToggles import GenToggles,BgeoToggles,m_to_km, R_REF, runFullSimulation
from numpy import array,degrees,arccos
from numpy.linalg import norm
from datetime import datetime
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from myspaceToolsLib.models import CHAOS


########################################
# --- GENERATE THE B-FIELD & TOGGLES ---
########################################
plot_BField = True

# --- OUTPUT DATA ------
outputData = False if not runFullSimulation else True


def generateGeomagneticField(outputData, **kwargs):
    plotting = kwargs.get('plotting', False)

    def geomagneticFieldProfile(altRange, **kwargs):
        plotBool = kwargs.get('showPlot', False)

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
            figure_width = 12  # in inches
            figure_height = 8  # in inches
            Title_FontSize = 25
            Label_FontSize = 25
            Tick_FontSize = 25
            Tick_FontSize_minor = 20
            Tick_Length = 10
            Tick_Width = 2
            Tick_Length_minor = 5
            Tick_Width_minor = 1
            Plot_LineWidth = 2.5
            Legend_fontSize = 20
            dpi = 100

            fig, ax = plt.subplots(2,sharex=True)
            fig.set_size_inches(figure_width, figure_height)
            ax[0].plot(altRange/R_REF, Bgeo/(1E-9),linewidth=Plot_LineWidth)
            ax[0].set_title('|B| vs Altitude', fontsize=Title_FontSize)
            ax[0].set_ylabel('$B_{geo}$ [nT]', fontsize=Label_FontSize)
            ax[0].set_yscale('log')
            ax[0].axvline(x=400000/R_REF,label='Observation Height',color='red',linewidth=Plot_LineWidth)
            ax[0].legend(fontsize=Legend_fontSize)

            ax[1].plot(altRange / R_REF, Bgrad/(1E-9),linewidth=Plot_LineWidth)
            ax[1].set_title(r'$\nabla B$ vs Altitude', fontsize=Title_FontSize)
            ax[1].set_ylabel(r'$\nabla B$ [nT/m]', fontsize=Label_FontSize)
            ax[1].set_xlabel('Altitude [$R_{E}$]', fontsize=Label_FontSize)
            ax[1].axvline(x=400000 / R_REF, label='Observation Height', color='red',linewidth=Plot_LineWidth)
            ax[1].legend(fontsize=Legend_fontSize)

            for i in range(2):
                ax[i].tick_params(axis='y', which='major', labelsize=Tick_FontSize, width=Tick_Width,
                                           length=Tick_Length)
                ax[i].tick_params(axis='y', which='minor', labelsize=Tick_FontSize_minor,
                                           width=Tick_Width_minor, length=Tick_Length_minor)
                ax[i].tick_params(axis='x', which='major', labelsize=Tick_FontSize, width=Tick_Width,
                                           length=Tick_Length)
                ax[i].tick_params(axis='x', which='minor', labelsize=Tick_FontSize_minor,
                                           width=Tick_Width_minor, length=Tick_Length_minor)

            plt.tight_layout()
            plt.savefig('C:\Data\ACESII\science\simulations\TestParticle\geomagneticField\MODEL_Bgeo.png',dpi=dpi)
            # plt.show()

        return Bgeo, Bgrad


    if plotting:
        # get all the variables and plot them if required
        geomagneticFieldProfile(altRange=GenToggles.simAlt, showPlot=plotting)

    if outputData:

        # get all the variables and plot them if required
        Bgeo, Bgrad = geomagneticFieldProfile(altRange=GenToggles.simAlt)

        from copy import deepcopy
        from myspaceToolsLib.CDF_load import outputCDFdata

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

        outputPath = rf'{GenToggles.simOutputPath}\geomagneticField\geomagneticfield.cdf'
        outputCDFdata(outputPath, data_dict)



#################
# --- EXECUTE ---
#################

generateGeomagneticField(outputData=outputData,plotting=plot_BField)