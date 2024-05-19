from ACESII_code.myImports import *
from ACESII_code.class_var_func import Re

from ACESII_code.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes

for time in dispersionAttributes.keyDispersionDeltaT:

    print( (pycdf.lib.datetime_to_tt2000(time[1]) - pycdf.lib.datetime_to_tt2000(time[0]))/1E9)
