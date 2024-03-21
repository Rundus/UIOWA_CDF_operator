import matplotlib.pyplot as plt
from ACESII_code.myImports import *
from numpy.fft import rfft, fftfreq

a = ['E_East','E_North','E_Up','Emag','ILat','ILong']
coordinatesSets = [['_east','_north','_up'],['_x','_y','_z'],['_e','_p','_r']]



for set in coordinatesSets:

    coords = []

    for key in a:

        for coordStr in set:
            # print(key.lower(),coordStr, coordStr in key.lower())
            if coordStr in key.lower():
                coords.append(key)

    print(coords)
