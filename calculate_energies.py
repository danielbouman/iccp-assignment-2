"""

"""
## Import libraries
import numpy as np
def calculate_energies(possible_beads_pos,beads_pos,epsilon,sigma):
    energies  = np.zeros(len(possible_beads_pos),dtype=float) # initialize all possible energies corresponding to each possible new bead position
    for ii in range(0, len(possible_beads_pos)-1):
        for iii in range(0,len(beads_pos)-1):
            abs_distance_squared = np.divide(sum(np.square(np.subtract(possible_beads_pos[ii,:],beads_pos[iii,:]))),np.square(sigma)) # distance alluar
            V = 4*epsilon*(np.power(abs_distance_squared,-6)-np.power(abs_distance_squared,-3))
            energies[ii] = energies[ii]+V

    energies = np.divide(energies,2)
    return energies;