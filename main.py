"""

"""
## Import modules
import simulation               # simulation module
import numpy as np
import save_data as save
import os

# Simulation
sigma, epsilon, T, minBeads, maxBeads, plot_data, bending_energy, amount_of_polymers = simulation.user_input()  # User input for simualation variables

# Check if multiple polymer lengths are simulated
if (minBeads+1) == maxBeads:
    write_mode="w"
    multi = False
else:             
    # os.remove("R_squared.dat")
    write_mode="a"
    multi = True

for iii in range(minBeads, maxBeads):
    simulation.simulation(iii,multi,write_mode)


if plot_data == 'y':                                            # Plot when chosen
    simulation.plot(beads_pos,end_to_end_distance)