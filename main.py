"""

"""
## Import modules
from datetime import datetime   # timer functions
import running as start         # start message
import simulation               # simulation module        

# Simulation
sigma, epsilon, T, Nbeads, plot_data = simulation.user_input()  # USer input for simualation variables
start_time=datetime.now() # =========== Start timer
start.message()                                                 # Startup message
beads_pos = simulation.start(Nbeads,sigma,epsilon,T)            # Start simulation
if plot_data == 'y':                                            # Plot when chosen
    simulation.plot(beads_pos)
print(datetime.now() - start_time) # =========== End timer