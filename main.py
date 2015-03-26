"""

"""
## Import modules
from datetime import datetime   # timer functions
import simulation               # simulation module        

# Simulation
sigma, epsilon, T, Nbeads, plot_data, bending_energy, amount_of_polymers = simulation.user_input()  # User input for simualation variables
start_time=datetime.now() # =========== Start timer                                             # Startup message
beads_pos = simulation.start(Nbeads,sigma,epsilon,T,bending_energy)            # Start simulation
end_to_end_distance = np.zeros((amount_of_polymers,1),dtype=float)    # initialize all end_to_end distances
start_time=datetime.now() # =========== Start timer                                             # Startup message
for ii in range(0, number_of_beads):
    beads_pos,end_to_end_distance[ii] = simulation.start(amount_of_beads,sigma,epsilon,T)            # Start simulation

print(datetime.now() - start_time) # =========== End timer
if plot_data == 'y':                                            # Plot when chosen
    simulation.plot(beads_pos)