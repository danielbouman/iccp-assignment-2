"""

"""
## Import modules
import simulation
import numpy as np
import save_data as save
import os
import sys

sys.setrecursionlimit(10000)
np.seterr(divide='ignore', invalid='ignore')

# Simulation parameters
minBeads, maxBeads, minT, maxT, stepT, minBend, maxBend, stepBend, plot_data = simulation.variables()

print(minBend)
print(maxBend)
# Actual simulation
for nBeadsVar in range(minBeads, maxBeads):
    print('beads')
    for tempVar in simulation.drange(minT, maxT,stepT):
        print('temp')
        for bendVar in simulation.drange(minBend, maxBend, stepBend):
            existingPos = simulation.start(nBeadsVar,tempVar,bendVar)
# Plot data
if plot_data == 'y':
    simulation.plot(existingPos)