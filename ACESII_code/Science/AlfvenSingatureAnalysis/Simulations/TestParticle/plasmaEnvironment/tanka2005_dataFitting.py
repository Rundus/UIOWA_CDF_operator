# --- File_template_ACESII.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
plotDigitizedData = False
fitCurveToData = True
plotFit = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

Altitude_Total = np.array([499.99999999999994, 599.4468440705969,  893.4456945885175, 739.8373622075442, 1094.718645036434, 1322.009011920112, 1596.4904183575156, 1942.0009769047372, 2311.4201281269657, 2791.327482702818, 3370.875342340215,  3785.840592147318, 4221.149441025492, 4810.085685325583, 5285.8978140377285,  5725.089226678763,  6155.942138747528, 6764.885954904924,  7380.320174706805,  7935.740765263604,   8532.960630799225,   9309.245109713238,  8977.560293719682,  10156.151922220071,   10763.140371047517,  6429.865535454216, 5024.1219741141385, 4538.820390234139, 3983.0972682257147, 3572.3377082941183, 3067.4446993938636, 2540.0650636771625, 2118.674148343684])
Density_Total = np.array([7442.01675466667,    5661.078983820722, 2491.8673895899983, 3829.985877007741, 1621.260048131898, 975.5286453470695,   586.985498710856,  314.1261601199634,  168.1050804294565, 76.94435356339415, 27.85809704090042, 13.259207513855252,  6.56226693804075, 2.569028915108545, 1.4296045766291532, 0.9671946556202097, 0.7650550092510098, 0.629276686421801, 0.5596687134504675, 0.5382212313898729, 0.49776048529048134, 0.46034137315021834, 0.4786854346003639, 0.46034137315021834,  0.46034137315021834, 0.6804278022455501, 1.8793486616648172, 3.948582698900836,  9.699647111541799, 19.598354513371504, 50.061575059253656,  118.2630298443682, 238.95310396987702])

Altitude_Op = np.array([503.6411782616502, 561.5515567417896, 630.6803017582541, 688.0563228832785, 778.3856991899237, 880.5736922297984, 1018.0995499443744, 960.6837172720645, 1094.718645036434, 1229.481918628943, 1370.8519382290901, 1495.5649339383822, 1631.6236708360916, 1872.8082304368338, 2043.1866048574295, 2229.065119654618, 2414.272305493529, 2595.963142468009, 2751.1124207702, 2936.7662479340643, 3089.782925310358, 3227.270220497889, 3445.056692327947, 3650.9525507015987, 4012.1036025999188, 4408.979600379052, 4845.114444193743, 5247.682320459124, 5561.313169626089, 5851.0787114576415, 6111.436479434778, 6245.928135864294, 3869.1538972809585, 4190.63176644437, 3572.3377082941183, 2328.2527135751193, 1741.730832104719, 1160.1451547363683, 830.9136449667715, 739.8373622075442, 595.1130188953433])
Density_Op = np.array([5886.666312725545, 5235.507739034487, 4141.308472338776, 3406.3287494518627, 2591.1653697116644, 1895.5423744283125, 1333.5265289162521, 1621.260048131898, 1054.8250499402072, 742.0763610652137, 522.055601242125, 381.9048083651314, 279.3788292001547, 155.46779193061766, 105.18126475988069, 73.99571162030468, 46.298199406870935, 32.57108782320815, 23.827069422816194, 16.762500068740252, 12.262447450503846, 8.97047677015572, 6.310789413841213, 3.948582698900836, 2.1972940297983485, 1.087488068184512, 0.5382212313898729, 0.27699208196975655, 0.16666895036806434, 0.10428269598473719, 0.07336356257638207, 0.05580703940486837, 2.7778539016792294, 1.6714629197207898, 4.991860764875156, 56.28790981744149, 212.52110424471556, 938.1446847685236, 2304.5414204839553, 3029.5348134305505, 4656.377621784846])

Altitude_Hp = np.array([503.6411782616502, 549.4598391629211, 590.8105260072398, 639.8994260671726, 813.0218360905542, 698.1141552819029, 750.6521166743614, 867.8871386854019, 939.9976091980864, 1003.4316413208621, 1078.9468739412332, 1160.1451547363683, 1238.435444299352, 1322.009011920112, 1452.7817181231128, 1667.5300851298155, 1832.4816247495214, 2013.750117607759, 2181.0673435498607, 2467.4021191683287, 2731.222683444829, 3001.3943527274355, 3250.772352840586, 3546.510752580111, 3897.3304554044767, 4282.853130382592, 4638.704185897042, 5024.1219741141385, 5561.313169626089, 6245.928135864294, 7065.905765899767, 7764.8628913432785, 8410.025023149248, 9108.791987666478, 9865.617670125288, 10531.38096789753, 5247.682320459124, 5893.68835268037, 6667.423256841222, 4070.751513038073, 1550.8200548944383])
Density_Hp = np.array([302.0882973235955, 326.6437175684845, 353.19513921013396, 381.9048083651314, 367.2695494531804, 381.9048083651314, 367.2695494531804, 367.2695494531804, 367.2695494531804, 353.19513921013396, 339.6600849359795, 339.6600849359795, 314.1261601199634, 290.51174644295224, 258.37654386736995, 212.52110424471556, 174.8038698612576, 138.27059061659116, 113.73098410005477, 76.94435356339415, 52.05646985416136, 33.86900730508808, 22.913974529853924, 13.787571314288687, 7.378439235465839, 3.797266117208541, 2.375902333335134, 1.374819633726121, 0.7955415501552152, 0.581970852406614, 0.49776048529048134, 0.4786854346003639, 0.4786854346003639, 0.46034137315021834, 0.44270028815634965, 0.4257352403338333, 1.0458136270857745, 0.6543526211045144, 0.5382212313898729, 5.190780536076251, 238.95310396987702])



# Sort the lists
def sortListBasedOnAnother(List1,List2):
    combined_lists = list(zip(List1, List2))
    sorted_lists = sorted(combined_lists, key=lambda x: x[0])
    sorted_List2 = list([item[1] for item in sorted_lists])
    sorted_List1 = list(sorted(List1))
    return sorted_List1, sorted_List2


Altitude_Total, Density_Total = sortListBasedOnAnother(Altitude_Total, Density_Total)
Altitude_Op, Density_Op = sortListBasedOnAnother(Altitude_Op, Density_Op)
Altitude_Hp, Density_Hp = sortListBasedOnAnother(Altitude_Hp, Density_Hp)

functionalSum = np.array(Density_Hp) + np.array(Density_Op)

if plotDigitizedData:
    fig, ax = plt.subplots()
    ax.plot(Altitude_Total, Density_Total,label='Total', color='black')
    ax.plot(Altitude_Op, Density_Op,label='Op', color='red')
    ax.plot(Altitude_Hp, Density_Hp,label='Hp', color='blue')
    ax.plot(Altitude_Hp, functionalSum,label='FunctionalSum', color='orange')

    ax.set_ylabel('n_i [cm^-3]')
    ax.set_xlabel('Altitude [km]')
    ax.legend()
    ax.set_yscale('log')
    ax.set_ylim(1E-2, 1E5)
    ax.set_xscale('log')
    ax.set_xlim(500,11000)
    ax.grid(True)
    plt.show()

if fitCurveToData:

    # def fitFunc(x,n0,n1,z0,h,H,a):
    #     return n0 * np.exp(-1 * (x - z0) / h) + np.tanh(x/a)*n1 * (x**(H))
    def fitFunc(x,n0,n1,z0,h,H,a):
        return a*(n0 * np.exp(-1 * (x - z0) / h) +  n1*(x**(H)))
    xData = np.array(Altitude_Total)
    xData = [xData[0]] + [xData[i]*1.3 for i in range(len(xData)) if i != 0]
    yData = 1.4*np.array(Density_Total)
    guess = np.array([          2.2E7,        5E6,       500,       681,           -0.68,        0.0005])
    boundVals =      [    [1E7,2.5E7],  [2E6,1E7], [490,600], [680,700], [-0.74, -0.65],  [1E-4, 1E-2]]
    bounds = tuple([[boundVals[i][0] for i in range(len(boundVals))], [boundVals[i][1] for i in range(len(boundVals))]])
    params, cov = curve_fit(fitFunc, xdata=xData, ydata=yData, p0=guess, bounds=bounds, maxfev=1E6)
    newXData = np.linspace(np.array(xData).min(),20000, 2000)
    newYData = np.array([fitFunc(z, *params) for z in newXData])

    print('--- FIT PARAMS ---')
    print(f'n_0  = {params[0]}')
    print(f'n_1  = {params[1]}')
    print(f'z_0  = {params[2]}')
    print(f'h  = {params[3]}')
    print(f'H  = {params[4]}')
    print(f'a  = {params[5]}')

    if plotFit:
        fig, ax = plt.subplots()
        ax.scatter(xData, yData, label='Tanaka', color='blue')
        ax.plot(newXData, newYData, label='fit', color='black')
        ax.set_ylabel('n_i [cm^-3]')
        ax.set_xlabel('Altitude [km]')
        ax.legend()
        ax.set_yscale('log')
        ax.set_ylim(1E-2, 1E5)
        ax.set_xscale('log')
        ax.set_xlim(500, newXData.max())
        ax.grid(True)
        plt.show()