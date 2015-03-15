"""

"""
## Import modules
from datetime import datetime           # timer functions    
import sys                              # progress messages
import running as start                 # startup message
import simulate

## User input
sigma = input('Sigma value of L-J potential (default: 0.5): ') or 0.8
epsilon = input('Epsilon value of L-J potential (default: 0.5): ') or 7.25
T = input('Temperature, expressed in epsilon (default: 1): ') or 1
number_of_beads = input('Amount of beads per polymer: ') or 250
plot_data = input('Plot data? (y/n, default: y): ') or 'y'

# =========== Start timer
start_time=datetime.now()

# =========== Start simulation
start.message()
beads_pos = simulate.start(number_of_beads,sigma,epsilon,T)

if plot_data == 'y':
    import matplotlib.pyplot as plt 			           # plotting tools
    plt.plot(beads_pos[:,0],beads_pos[:,1], 'b')
    plt.plot(beads_pos[:,0],beads_pos[:,1], '.r')
    #plt.axis([-1*(number_of_beads), number_of_beads, -1*(number_of_beads), number_of_beads])
    plt.show()


# =========== End timer
print(datetime.now() - start_time)