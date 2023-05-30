# # --- scratchPaper.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering


test = [0,11231,2,5,1,23,1,3,12,3]

preparestring = '[-10]'

pads = ['-10','0','10','110']

for thing in test:
    preparestring = preparestring + f'{thing:8}'

print(preparestring)