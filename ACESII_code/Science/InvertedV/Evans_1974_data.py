# # --- Evans_1974_data.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: Provides needed functions and data needed to run the Evans 1974 model

# --- Imports ---
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import simpson

# --- TOGGLES ---
plot_Evans1974Curves = False
plot_Example_BackScatterCurve = False


########################
# --- SECONDARIES ---
########################

# EVANs - Secondary Electron Backscatter Curve
Energy_secondaries = [10.17496576561469, 10.718446060311226, 11.290955526950686, 11.997646044427388, 12.85961219779013, 13.546488791572546, 14.270053852147703, 15.032267038074046, 15.97312280060254, 17.12070551529857, 18.350735858020062, 19.499291530622624, 20.540816480342805, 21.826447283974872, 22.79373225515666, 23.803884496641345, 24.64414455946939, 25.960473995484463, 27.110967065377004, 28.31244665824677, 29.311855734157486, 30.610872342779636, 31.691414198105562, 32.81009840651562, 33.96827136573758, 35.47364683805397, 36.7258411479822, 38.022237019621116, 39.70727251629559, 41.46698390917224, 43.3046806178358, 46.82018711552433, 45.22381872576847, 48.89512229647873, 51.06201259915312, 53.78941220859147, 56.173203274592055, 58.66263706121415, 62.33427753861964, 65.6637680063783, 69.17109813495638, 73.50045352199562, 77.4263682681127, 81.56197976921872, 85.17657064716352, 89.72614740082682, 95.34202440844852, 102.19183459028312, 107.65025574452392, 114.3879862011722, 119.45733067563364, 126.93405508837688, 133.71404425547448, 139.63986368735334, 148.37979427108883, 156.30527492380824, 163.23227231294564, 178.02082518094767, 171.95107545235493, 189.16298482974503, 201.00252199891452, 211.73876548190316, 223.04846905374944, 234.96226321142643, 249.66833507399298, 267.60565824583045, 286.8316813342009, 302.15235860140325, 321.0637967119472, 344.13049869757356, 365.66930997210375, 385.2009790629936, 409.3103539141118, 434.9287123538119, 454.2034940680203, 470.2365518874762, 486.8355651572963, 499.66822499934193, 512.8391451707628, 530.9419994508138, 549.6838715128936, 564.1731706697959, 584.0880792596014, 594.3076210570079, 609.9731724020412, 620.6456147134156, 642.5539601049684, 653.7964546628194, 682.7708212456762, 700.7682033726139, 719.239984453556, 738.1986693275925, 751.1146188630573, 764.2565532984336, 784.4018448146026, 798.1261917473589, 826.2994743104199, 840.7568863253938, 862.9186754473221, 885.6646345059752, 901.1607335913945, 932.9710616381673, 949.2948612477545, 982.8042895037572]
NFlux_up_PeriE_secondaries = [0.01, 0.00962040327106476, 0.008903890941296734, 0.008082822192522706, 0.007337467970651373, 0.0067909850299557074, 0.006285203269240097, 0.005817091329374358, 0.005383843622033489, 0.004887374631624418, 0.004436687330978606, 0.004027559979850294, 0.003586095482011829, 0.003255405118685032, 0.002898577214651068, 0.0025808615404180717, 0.0022979709690469606, 0.0020068781676649667, 0.0017190722018585747, 0.0014725404276241208, 0.0012613637220374424, 0.0010597662486760708, 0.0009620403271064761, 0.0008240743309893649, 0.0007337467970651374, 0.000653320145959377, 0.0005817091329374357, 0.0005280670574582637, 0.0004701851489299847, 0.00041864772882899473, 0.0003727593720314938, 0.0003071814301268697, 0.00033190087959146774, 0.0002788548171726292, 0.0002580861540418071, 0.00022979709690469603, 0.0002126821747296208, 0.00019306977288832496, 0.00017526592085367583, 0.0001622123939129173, 0.00015013107289081742, 0.00013628679593557727, 0.00012613637220374423, 0.00011674193588234575, 0.00010804718223181286, 0.00009808365445406671, 0.00008733261623828437, 0.00007626985859023451, 0.0000692366640267872, 0.00006164757056337162, 0.00005489032450921557, 0.00004982863514647296, 0.00004436687330978606, 0.000039503780135666406, 0.00003586095482011829, 0.00003131831005243845, 0.00002843030459302665, 0.000023886430789846003, 0.000025808615404180768, 0.00002168375310987433, 0.000019306977288832496, 0.000017526592085367586, 0.000015910384392721502, 0.000014443214647272789, 0.000013628679593557728, 0.00001213482767749143, 0.000011015819387358994, 0.000010195378685327302, 0.00000943604310147889, 0.000008733261623828437, 0.000008082822192522723, 0.000007337467970651388, 0.000006408002744415996, 0.000005705615758781021, 0.000005080218046913012, 0.000004436687330978607, 0.0000039503780135666406, 0.0000033838551534282333, 0.000002955209235202888, 0.0000025314033152415694, 0.000002168375310987433, 0.0000018574090746385719, 0.0000016221239391291696, 0.0000014443214647272788, 0.0000012860080806105652, 0.0000011450475699382812, 9.80836544540667e-7, 8.565902153685499e-7, 7.480826455225093e-7, 6.53320145959377e-7, 5.705615758781021e-7, 5.080218046913023e-7, 4.4366873309786064e-7, 3.9503780135666407e-7, 3.51737350096938e-7, 3.0718143012686963e-7, 2.843030459302665e-7, 2.4353120734334173e-7, 2.1268217472962083e-7, 1.8936988889503602e-7, 1.6861288382868804e-7, 1.4725404276241208e-7, 1.286008080610565e-7, 1.1231045018329516e-7]
Energy_secondaries, NFlux_up_PeriE_secondaries = zip(*sorted(zip(Energy_secondaries,NFlux_up_PeriE_secondaries)))

def integrate_SecondariesCurve():
    return simpson(x= Energy_secondaries,y= NFlux_up_PeriE_secondaries)

def generate_SecondariesCurve():
    return CubicSpline(np.array(Energy_secondaries), np.array(NFlux_up_PeriE_secondaries))

R = integrate_SecondariesCurve()

def calcSecondaries(SampledEnergies,InputOmniFlux):

    # InputOmniFlux is the number flux (not diffNFlux) array
    # This function returns the Upgoing number flux as a function of Energy
    # by multiplying the Evans Secondaries Curve by Omniflux and sampling it at "SampledEnergies"

    curve = 5



########################
# --- BACKSCATTER ---
########################

# EVANS - Primary Electron Backscatter Curve
Energy_backscatter = [1, 0.9582558549602653, 0.9261187281287937, 0.8874586936177703, 0.8576958985908941, 0.8289312615337154, 0.8011313071180073, 0.767688765637801, 0.7356422544596414, 0.6989473207273487, 0.6697703623957719, 0.6363612169654032, 0.5994842503189411, 0.5695810810737686, 0.5365740058950261, 0.5098088935466992, 0.48437886495410437, 0.46021732414665023, 0.4372609970601657, 0.4119218362588828, 0.3913745601980384, 0.37185221293765525, 0.3533036694992736, 0.33568035509467253, 0.3189361179185748, 0.3030271082866398, 0.29037750072735075, 0.27355020025126336, 0.2576980374514879, 0.24276450335386565, 0.23065504755640526, 0.2172886476740597, 0.2046968271807521, 0.192834699402793, 0.1801173528334133, 0.1711328304161781, 0.16121572750178864, 0.15187331811625285, 0.14429765031006772, 0.13709986812211702, 0.13026112205370022, 0.12271252398511903, 0.11659144011798317, 0.10983499655437803, 0.10347008713411984, 0.09830884473994818, 0.09340505282049072, 0.08950593874810804, 0.08504124891946163, 0.08149127469020741, 0.07808949110006652, 0.07419427072122865, 0.07170600970409613, 0.06930119775694155, 0.06584435056818477, 0.06309573444801933, 0.060461856957832634, 0.05793792843161315, 0.05551935914386209, 0.05320175096324757, 0.05098088935466992, 0.048437886495410414, 0.046415888336127795, 0.04447829676127633, 0.04262158829015325, 0.04084238652674522, 0.03913745601980384, 0.037503696379226896, 0.03593813663804628, 0.034145488738336005, 0.03300034791125285, 0.03189361179185748, 0.030823992397451434, 0.02953727118810853, 0.02830426305555415, 0.027122725793320285, 0.02599051079393097, 0.024694065319965788, 0.023865897868585808, 0.023065504755640526, 0.022102654979706378, 0.021000141557086554, 0.0199526231496888, 0.019119717955005014, 0.018321581675572463, 0.01740767383330236, 0.016681005372000592, 0.015984671064343196, 0.015187331811625286, 0.014677992676220698, 0.01406527242105237, 0.013478129649084574, 0.012805820423628654, 0.012271252398511897, 0.011659144011798317, 0.011077568505097092, 0.010615144878732719, 0.01008564560714737]
NFlux_up_PeriE_backscatter =[0.00011659144011798335, 0.0001081230103192746, 0.00010026967115829288, 0.00009125021956617401, 0.00008462240814424688, 0.0000770104582390694, 0.00007141692874235842, 0.0000662296761714834, 0.00006141919126211213, 0.000058042042927286296, 0.00005382625345418286, 0.00004898447745345474, 0.0000454265695305298, 0.00004212708446819444, 0.0000383376714314434, 0.000035553076963491, 0.000032970736989915416, 0.00003115783344567483, 0.000028894734843002547, 0.000026796012729932755, 0.000024849727887245227, 0.00002348335989328803, 0.00002177768420652337, 0.000020580233239763535, 0.000019448624389373654, 0.000018379237311466208, 0.000017368650727593748, 0.00001641363147908973, 0.000015511124183259427, 0.000014658241458327077, 0.000013852254685866867, 0.000013090585281163766, 0.0000123707964435753, 0.000011690585360501962, 0.000011258019285837975, 0.000010638994510302386, 0.000010245338593872244, 0.000009866248431789506, 0.00000950118507318144, 0.00000932375049624682, 0.000008978760230238898, 0.000008646535029500377, 0.000008326602570874948, 0.000007868762416533264, 0.000007436096708208832, 0.000006895987907029672, 0.000006395109030978996, 0.000005930610678191438, 0.000005499850408475976, 0.000005005128186391172, 0.0000045549071886755406, 0.000004224068906483475, 0.000003917260525325609, 0.0000036327366723874875, 0.0000034329898122952065, 0.000003124185713602668, 0.000002897265560913933, 0.0000026868273847837817, 0.0000024451420451043524, 0.0000022675431258708044, 0.000002102843815547939, 0.000001950107172003647, 0.0000017746912457903893, 0.0000016150543227062977, 0.0000014697770507853085, 0.0000013630221830031345, 0.0000012640212815719117, 0.0000011288378916846906, 0.0000010468465857176266, 9.526807029019838e-7, 8.834842880970382e-7, 8.193138424413288e-7, 7.59804312832608e-7, 6.914584397494283e-7, 6.292604106421371e-7, 5.835551032264557e-7, 5.211455867173361e-7, 4.7426752939905745e-7, 4.3981987805811335e-7, 4.078742758969022e-7, 3.7118519290075845e-7, 3.4422475956860925e-7, 3.1326107414037545e-7, 2.850826323317469e-7, 2.643761185749101e-7, 2.4059497342849544e-7, 2.2311974848625545e-7, 2.06913808111479e-7, 1.9188495988212056e-7, 1.7794770762290308e-7, 1.6502276503431532e-7, 1.5303660464837528e-7, 1.4192103954525943e-7, 1.3161283545126592e-7, 1.2205335101141175e-7, 1.1318820419024433e-7, 1.0496696290308788e-7, 9.734285811778847e-8]
Energy_backscatter, NFlux_up_PeriE_backscatter = zip(*sorted(zip(Energy_backscatter,NFlux_up_PeriE_backscatter)))

def generate_BackscatterCurve(E_Incident):
    return CubicSpline(np.array(Energy_backscatter)*E_Incident, (E_Incident/10000)*np.array(NFlux_up_PeriE_backscatter))










#############################
# --- Diagnostic Plotting ---
#############################
if plot_Evans1974Curves:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,2)
    fig.set_size_inches(12, 8)
    xData, yData = zip(*sorted(zip(Energy_backscatter,NFlux_up_PeriE_backscatter)))
    ax[0].set_title('primary electron backscatter')
    ax[0].plot(xData, yData,color='black')
    ax[0].set_yscale('log')
    ax[0].set_ylim(1E-10,1E-3)
    ax[0].set_xscale('log')
    ax[0].set_xlim(1E-2, 1)
    ax[0].set_ylabel(r'Units of $( \frac{10000}{E(Incident)} ) \cdot (cm^{-2}sec^{-2}eV^{-1}$) Upgoing Flux per Incident Electron')
    ax[0].set_xlabel('E in Units of E(Backscatter)/E(Incident)')

    xData, yData = zip(*sorted(zip(Energy_secondaries, NFlux_up_PeriE_secondaries)))
    ax[1].set_title('secondary electron backscatter ' )
    ax[1].plot(xData, yData, color='black')
    ax[1].set_yscale('log')
    ax[1].set_ylim(1E-8, 1E-1)
    ax[1].set_xscale('log')
    ax[1].set_xlim(1E1, 1E3)
    ax[1].set_ylabel(r'Upgoing Flux $(cm^{-2}sec^{-2}eV^{-1}$) per Incident Electron')
    ax[1].set_xlabel('Energy (eV)')

    plt.show()

if plot_Example_BackScatterCurve:
    import matplotlib.pyplot as plt
    Eincdient_test = 100
    Erange = np.linspace(1, Eincdient_test, 100)
    test_backscatterCurve = generate_BackscatterCurve(Eincdient_test)

    fig, ax = plt.subplots()
    ax.plot(Erange, test_backscatterCurve(Erange))
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.show()

