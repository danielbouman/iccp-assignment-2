"""

"""
## Import modules
from datetime import datetime   # timer functions
import simulation               # simulation module        

# Simulation
sigma, epsilon, T, Nbeads, plot_data = simulation.user_input()  # User input for simualation variables
start_time=datetime.now() # =========== Start timer                                             # Startup message
beads_pos = simulation.start(Nbeads,sigma,epsilon,T)            # Start simulation
print(datetime.now() - start_time) # =========== End timer
if plot_data == 'y':                                            # Plot when chosen
    simulation.plot(beads_pos)
