"""

"""
## Import libraries
import numpy as np
import list_tracking as track
def calculate_energies2(possible_beads_pos,beads_pos,epsilon,sigma_squared):
    energies  = np.zeros(len(possible_beads_pos),dtype=float) # initialize all possible energies corresponding to each possible new bead position
    for ii in range(len(possible_beads_pos)):
        beads_in_range = track.get(possible_beads_pos[ii],2)
        V = 4*epsilon*(np.power(abs_distance_squared,-6)-np.power(abs_distance_squared,-3))
        energies[ii] = energies[ii]+V
    return energies;