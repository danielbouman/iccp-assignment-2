"""

"""
## Import libraries
import numpy as np
def calculate_energies(possible_beads_pos,beads_pos,epsilon,sigma_squared):
    energies  = np.zeros(len(possible_beads_pos),dtype=float) # initialize all possible energies corresponding to each possible new bead position
    for ii in range(len(possible_beads_pos)):
        for iii in range(len(beads_pos)):
            abs_distance_squared = np.divide(sum(np.square(np.subtract(possible_beads_pos[ii,:],beads_pos[iii,:]))),sigma_squared) # distance
            V = 4*epsilon*(np.power(abs_distance_squared,-6)-np.power(abs_distance_squared,-3))
            #V = 1/(abs_distance_squared)
            energies[ii] = energies[ii]+V
    return energies;