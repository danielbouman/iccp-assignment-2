"""

"""
## Import libraries
import numpy as np
import list_tracking as track
def calculate_energies(possible_beads_pos,beads_pos,epsilon,sigma_squared,cutoff_length):
    energies  = np.zeros(len(possible_beads_pos),dtype=float) # initialize all possible energies corresponding to each possible new bead position
    relevant_beads = track.get(beads_pos[-1,:],2*1+)
    for ii in range(len(possible_beads_pos)):
        for iii in relevant_beads: #range(len(beads_pos)):
            abs_distance_squared = (sum(np.square(np.subtract(possible_beads_pos[ii,:],beads_pos[iii,:]))))/sigma_squared # distance
            V = 4*epsilon*(((abs_distance_squared)**(-6))-(abs_distance_squared)**(-3))
            energies[ii] = energies[ii]+V
    return energies;