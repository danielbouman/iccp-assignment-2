"""

"""

## Import libraries
import numpy as np
def positions(previous_beads_pos, angles):
    offset = 2*np.pi*(np.random.random())               # create random offset of each angle to avoid preferential angles
    angles = np.add(angles,offset)                      # add offset to all possible angles
    possible_beads_pos = np.zeros((len(angles),2),dtype=float)      # initialize the possible bead possitions, wrt the origin
    possible_rel_prev_beads_pos = np.zeros((len(angles),2),dtype=float) # initialize the possible bead positions, relative to the previous bead

    possible_rel_prev_beads_pos[:,0] = np.cos(angles)         # calcate x-position relative to previous point
    possible_rel_prev_beads_pos[:,1] = np.sin(angles)         # calcate y-position relative to previous point
    return previous_beads_pos[:] + possible_rel_prev_beads_pos[:,:];
"""

"""
## Import libraries
import numpy as np
def roulette(energies,T):
    probabilities = np.exp(np.divide(energies,-T))   # unnormalized probabilities
    probabilities = np.divide(probabilities,sum(probabilities)) # normalize probabilities by dividing by their sum
    
    cumsum_probabilities = np.cumsum(probabilities)
    RNG = np.random.random() # rng is a random number chosen from a uniform distribution between 0 and 1. It is used to select one of the possible bead positions with corresponding probability

    for ii in range(0, len(probabilities)):
        if cumsum_probabilities[ii] > RNG:
            break
    return ii,probabilities[ii];