"""

"""
## Import modules
from datetime import datetime   # timer functions
import simulation               # simulation module
import numpy as np

# Simulation
sigma, epsilon, T, Nbeads, plot_data, bending_energy, amount_of_polymers = simulation.user_input()  # User input for simualation variables
start_time=datetime.now() # =========== Start timer                                             # Startup message
weight_factors = np.zeros((amount_of_polymers,1),dtype=float)    # initialize all end_to_end distances
beads_pos = simulation.start(Nbeads,sigma,epsilon,T,bending_energy)            # Start simulation
end_to_end_distance_squared = np.zeros((amount_of_polymers,1),dtype=float)    # initialize all end_to_end distances, squared
radius_of_gyration_squared = np.zeros((amount_of_polymers,1),dtype=float)    # initialize all end_to_end distances, squared

start_time=datetime.now() # =========== Start timer                                             # Startup message
for ii in range(0, amount_of_polymers):
    beads_pos,weight_factors[ii],end_to_end_distance_squared[ii],radius_of_gyration_squared[ii] = simulation.start(Nbeads,sigma,epsilon,T,bending_energy)            # Start simulation

print(datetime.now() - start_time) # =========== End timer

end_to_end_distance = np.sqrt(end_to_end_distance_squared)
radius_of_gyration = np.sqrt(radius_of_gyration_squared)

exp_end_to_end_distance = simulation.calculate_expectation_value(weight_factors,end_to_end_distance)
exp_radius_of_gyration = simulation.calculate_expectation_value(weight_factors,radius_of_gyration)

#print(weight_factors)
print(exp_end_to_end_distance)
print(exp_radius_of_gyration)

if plot_data == 'y':                                            # Plot when chosen
    simulation.plot(beads_pos,end_to_end_distance)