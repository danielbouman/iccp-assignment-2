"""

"""
## Import modules
from datetime import datetime   # timer functions
import simulation               # simulation module
import numpy as np
import save_data as save
import os

# Simulation
sigma, epsilon, T, minBeads, maxBeads, plot_data, bending_energy, amount_of_polymers = simulation.user_input()  # User input for simualation variables
start_time=datetime.now() # =========== Start timer                                             # Startup message
weight_factors = np.zeros((amount_of_polymers),dtype=float)    # initialize all end_to_end distances
end_to_end_distance_squared = np.zeros((amount_of_polymers),dtype=float)    # initialize all end_to_end distances, squared
radius_of_gyration_squared = np.zeros((amount_of_polymers),dtype=float)    # initialize all end_to_end distances, squared

# Check if multiple polymer lengths are simulated
if (minBeads+1) == maxBeads:
    write_mode="w"
    multi = False
else:             
    # os.remove("R_squared.dat")
    write_mode="a"
    multi = True

for iii in range(minBeads, maxBeads):
    start_time=datetime.now() # =========== Start timer
    
    # Simulate polymers
    for ii in range(0, amount_of_polymers):
        beads_pos,weight_factors[ii],end_to_end_distance_squared[ii],radius_of_gyration_squared[ii] = simulation.start(iii,sigma,epsilon,T,bending_energy)            # Start simulation
    
    print(str(iii) + " beads, done in: " + str(datetime.now() - start_time))
    
    # Collect properties
    radius_of_gyration = np.sqrt(radius_of_gyration_squared)
    exp_end_to_end_distance_squared = simulation.calculate_expectation_value(weight_factors,end_to_end_distance_squared)
    exp_radius_of_gyration = simulation.calculate_expectation_value(weight_factors,radius_of_gyration)
    
    # Save data to file
    save.save(exp_end_to_end_distance_squared,"R_squared",header="",write_mode=write_mode)
    save.save(exp_radius_of_gyration,"exp_radius_of_gyration",header="",write_mode=write_mode)

if plot_data == 'y':                                            # Plot when chosen
    simulation.plot(beads_pos,end_to_end_distance)