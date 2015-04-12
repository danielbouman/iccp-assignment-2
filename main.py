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
minBeads, maxBeads, plot_data = simulation.variables()  # User input for simualation variables
# Actual simulation
# print(minBeads)
# print(maxBeads)
for iii in range(minBeads, maxBeads):
    existingPos = simulation.start(iii)
# Plot data
if plot_data == 'y':
    simulation.plot(existingPos)